# Differential Isoform Usage
# tutorial: https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#overview-of-isoform-switch-workflow

# install packages if needed
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

# import libraries
library(IsoformSwitchAnalyzeR)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ggrepel)
library(dplyr)


# set working environment
working_dir   <- "C:/Users/irisg/Desktop/stringtie_quant/gse109843"
setwd(working_dir)

# Set variables for the analysiS
metadata_path <- "C:/Users/irisg/Desktop/stringtie_quant/gse109843/metadados_109843_subtypes.csv"
kallisto_path <- "../gse109843"
gtf_path      <- "C:/Users/irisg/Desktop/stringtie_quant/gse109843/merged_gse109843.gtf"

# 1. Data Import
metadata <- read.csv2(metadata_path, stringsAsFactors = FALSE)
quant_data <- importIsoformExpression(parentDir = kallisto_path,
                                      readLength = 300)


##Bloco para trabalhar os varios subtipos##
# 1. Definir o vetor de subtipos (excluindo 'Control')
# Usando a coluna 'disease' do ficheiro 
todos_subtipos <- unique(metadata$disease)
subtipos_alvo <- "PD1POS"
subtipos_alvo
# 2. Iniciar o Loop
for (sub in subtipos_alvo) {
  
  
  
  # A. Filtrar metadata para o Control + Subtipo atual
  metadata_sub <- metadata[metadata$disease %in% c("Control", "PD1POS"), ]
  
  # B. Identificar colunas a manter
  amostras_keep <- metadata_sub$id
  
  # Preservar as primeiras colunas
  cols_indices <- c(1, 2, which(colnames(quant_data$counts) %in% amostras_keep))
  
  quant_sub <- list(
    counts    = quant_data$counts[, cols_indices],
    abundance = quant_data$abundance[, cols_indices]
  )
  
  # C. Criar design matrix para iteração
  # IMPORTANTE: Usar factor e levels para garantir que Control é o grupo de referência
  design_sub <- data.frame(
    sampleID  = colnames(quant_sub$abundance)[-(1:2)],
    condition = factor(
      metadata_sub$disease[match(colnames(quant_sub$abundance)[-(1:2)], metadata_sub$id)],
      levels = c("Control", "PD1POS")
    )
  )
  design_sub <- na.omit(design_sub)
  
  # 3. Create switchAnalyzeRlist
  message("-> A criar switchAnalyzeRlist...")
  switch_list <- importRdata(
    isoformCountMatrix   = quant_sub$counts,       # Usar o subset
    isoformRepExpression = quant_sub$abundance,    # Usar o subset
    designMatrix         = design_sub,             # Usar o subset
    isoformExonAnnoation = gtf_path,
    showProgress         = FALSE
  )
  
  # 4. Perform Differential Transcript Usage (DTU) test
  message("-> A filtrar dados...")
  switch_list_filtered <- preFilter(
    switch_list,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0.1,
    removeSingleIsoformGenes = TRUE
  )

  # Using DEXSeq (Robust but slower)
  message("-> A correr teste estatístico (DEXSeq) - Isto pode demorar...")
  switch_list_stats <- isoformSwitchTestSatuRn(
    switch_list_filtered,
    reduceToSwitchingGenes = FALSE 
  )
  
  # 5. Guardar os Resultados
  nome_ficheiro <- paste0("IsoformSwitchResults_Control_vs_PD1POS.rds")
  saveRDS(switch_list_stats, file = nome_ficheiro)
  
  
  # Limpar a memória no final de cada loop
  gc()
}


load("C:/Users/irisg/Desktop/stringtie_quant/gse109843/IsoformSwitchResults_Control_PD1POS.RData")
# Summarize switching features

extractSwitchSummary(switch_list_stats)

IsoformSwitchResults_Control_vs_Tfh <- analyzeORF(
  switch_list_stats,
  genomeObject = Hsapiens,
  orfMethod = "longest"
)


# 5. Extract Fasta Sequences (Optional)
# This step allows for external functional annotation (Pfam, CPAT)
switch_list_stats <- extractSequence(IsoformSwitchResults_Control_vs_Tfh, 
                                     genomeObject = Hsapiens,
                                     dIFcutoff = 0.1,
                                     alpha = 0.05,
                                     pathToOutput = '.',
                                     writeToFile = T
)

save(switch_list_stats, file = "PD1POS_ORF_fasta.RData")
# 6.EXTERNAL SEQUENCE ANALYSIS

# This step requires CPC2 analysis!
# needs the output from CPC2 (input: isoforms_nt.fa)
# ran directly on the website: https://cpc2.gao-lab.org/batch.php

switch_list_stats <- analyzeCPC2(
  switchAnalyzeRlist   = switch_list_stats,
  pathToCPC2resultFile = "C:/Users/irisg/Desktop/stringtie_quant/gse109843/PD1POS/Real/result_cpc2_PD1POS.txt",
  codingCutoff         = 0.325, # the coding potential cutoff we suggested for human
  removeNoncodinORFs   = TRUE   # if ORF was predicted de novo
)

# 7. Run Alternative Splicing Analysis
switch_list_stats = analyzeIntronRetention(switch_list_stats, showProgress = T)
switch_list_stats = analyzeAlternativeSplicing(switch_list_stats, showProgress = T)


# 8. The consequences of Splicing Analysis and Coding Potential
consequencesOfInterest <- c('intron_retention', 'NMD_status', 'coding_potential')

# other things that can be analyzed: 'domains_identified','ORF_seq_similarity'

switch_list_stats <- analyzeSwitchConsequences(switch_list_stats,
                                               consequencesToAnalyze = consequencesOfInterest, 
                                               dIFcutoff = 0.1, # higher threshold, faster running times
                                               showProgress=T)

# 9. save results as R object
saveRDS(switch_list_stats, file = "SwitchAnalysisResults_PD1POS_Real.rds")

switch_list_stats <- readRDS ("C:/Users/irisg/Desktop/stringtie_quant/gse109843/SwitchAnalysisResults_PD1POS_Real.rds")
# Quick Check
switchPlot(switch_list, gene = 'XIST')
switchPlot(switch_list_stats, gene = 'XIST')

# 10. Consequeces 
extractConsequenceSummary(
  SwitchAnalysisResults_PD1POS_Real,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE # enables analysis of fraction of significant features
)

# 11. Consequeces Enrichment
extractConsequenceEnrichment(
  SwitchAnalysisResults_PD1POS_Real,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE
)

extractSplicingEnrichment(switch_list_stats)


