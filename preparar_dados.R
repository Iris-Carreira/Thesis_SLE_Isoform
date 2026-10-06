
# SCRIPT SETUP DADOS PARA ANALISE

# Objetivo: Ler os ficheiros raw do Kallisto, criar os dados e GRAVAR em disco.


if (!require("tximport", quietly = TRUE)) BiocManager::install("tximport")
library(tximport)
library(openxlsx)

# 1. Carregar a Anotação dos Genes (GeneInfo)
load('C:/Users/irisg/Desktop/Tese Lupus/gencode.v38_geneInfo.RData') 

# 2. Definir onde estão as pastas do Kallisto
dir_PBMC <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_228066"
dir_bcell <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_118256"
dir_tcell <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto109843"

# Função para listar ficheiros
get_files <- function(path) {
  folders <- list.dirs(path, full.names = FALSE, recursive = FALSE)
  folders <- folders[grep("RNA-seq", folders)] # Filtra pastas RNA-seq
  files <- file.path(path, folders, "abundance.tsv")
  names(files) <- folders
  return(files)
}

# CORREÇÃO: Usar os nomes corretos das variáveis
files_PBMC <- get_files(dir_PBMC)
files_b <- get_files(dir_bcell)
files_t <- get_files(dir_tcell)
all_files <- c(files_PBMC, files_b, files_t)


# 3. Criar a Tabela de Metadados (metaSuper)
# CORREÇÃO: Atualizado para os comprimentos dos novos ficheiros
batch_vec <- c(rep("Batch_PBMC", length(files_PBMC)), 
               rep("Batch_B", length(files_b)), 
               rep("Batch_T", length(files_t)))

names(all_files) <- paste(batch_vec, names(all_files), sep="_")
sample_names <- names(all_files)

# Definir Grupos (MÁXIMA GRANULARIDADE)
group_vec <- rep(NA, length(sample_names))

# A. Controlos
group_vec[grep("Control|Healthy|Contro", sample_names, ignore.case=T)] <- "Control"

# B. SLE (apanha todos os Lúpus de B e T)
group_vec[grep("SLE", sample_names, ignore.case=T)] <- "SLE"

# C. Especificar os graus de severidade no PBMC 
group_vec[grep("mSLE", sample_names, ignore.case=T)] <- "mSLE"
group_vec[grep("sSLE", sample_names, ignore.case=T)] <- "sSLE"

metaSuper <- data.frame(SampleID = sample_names, 
                        Batch = factor(batch_vec), 
                        Group = factor(group_vec), 
                        row.names = sample_names)

print("Resumo dos Metadados criados:")
print(table(metaSuper$Group, metaSuper$Batch))

# 4. Importar 

txi <- tximport(all_files, type = 'kallisto', txOut = F, 
                tx2gene = geneInfo[,1:2], 
                ignoreTxVersion = T, 
                countsFromAbundance = "lengthScaledTPM")

# Atribuir nomes de genes (
gene_symbols <- geneInfo[match(rownames(txi$abundance), geneInfo$Ensembl_GeneID),"GeneSymbol"]
rownames(txi$abundance) <- gene_symbols

# Criar logTPMs
logTPMs <- log2(txi$abundance + 1)

# Remover duplicados e linhas vazias
valid_genes <- !is.na(rownames(logTPMs)) & !duplicated(rownames(logTPMs))
logTPMs <- logTPMs[valid_genes, ]

#guardar também as contagens brutas limpas 
counts_matrix <- round(txi$counts)
counts_matrix <- counts_matrix[valid_genes, ]
rownames(counts_matrix) <- rownames(logTPMs)

print("Dados processados com sucesso!")

# 5. SALVAR TUDO NUM FICHEIRO .RData

save(logTPMs, counts_matrix, metaSuper, geneInfo, file = "Dados_SuperDataset_Prontos.RData")

