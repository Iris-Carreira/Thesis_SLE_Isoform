# Differential Isoform Usage
# tutorial: https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#overview-of-isoform-switch-workflow

library(IsoformSwitchAnalyzeR)
setwd("C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_228066")

# 1. Import metadata
metadata <- read.csv2("metadados_228066.csv", stringsAsFactors = FALSE)

# 2. Import quantifications (Ensure rhdf5 is installed)
kallistoQuant <- importIsoformExpression(parentDir = "../kalisto_228066")

# 3. Create myDesign with explicit matching

myDesign <- data.frame(
  sampleID  = colnames(kallistoQuant$abundance)[-1],
  condition = metadata$disease[match(colnames(kallistoQuant$abundance)[-1], metadata$id)]
)


# 4. Create switchAnalyzeRlist
Tcell_SwitchList <- importRdata(
  isoformCountMatrix   = kallistoQuant$counts,
  isoformRepExpression = kallistoQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/gencode.v39.chr_patch_hapl_scaff.annotation.gtf/gencode.v39.chr_patch_hapl_scaff.annotation.gtf",
  showProgress         = TRUE,
  removeNonConvensionalChr = TRUE
)

# 5. Filtering the noise
Tcell_SwitchList <- preFilter(
  switchAnalyzeRlist = Tcell_SwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

# 6. Identify Isoform Switches
Tcell_SwitchList <- isoformSwitchTestDRIMSeq(
  switchAnalyzeRlist = switch_list_filtered,
  reduceToSwitchingGenes = FALSE,
  showProgress = TRUE
)
save(Tcell_SwitchList, file = "SwitchList_AposDEXSeq228066.RData")

# 7. ORF analyzing
options(timeout = 3600)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update = FALSE)

library(BSgenome.Hsapiens.UCSC.hg38)
Tcell_SwitchList <- analyzeORF(
  Tcell_SwitchList,
  genomeObject = Hsapiens,
  orfMethod = "longest"
)

# Look at the isoform switch in a given gene
switchPlot(Tcell_SwitchList, gene='XIST')


# 8. summary + overlap
extractSwitchSummary(
  Tcell_SwitchList,
  filterForConsequences=TRUE
)

extractSwitchOverlap(
  Tcell_SwitchList,
  filterForConsequences=TRUE
)

# 9. Consequences
extractConsequenceSummary(
  Tcell_SwitchList,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)

extractConsequenceEnrichment(
  Tcell_SwitchList,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE
)
