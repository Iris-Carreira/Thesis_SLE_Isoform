# CORRELAÇÃO CLÍNICA (SLEDAI vs INTERFERÃO)


library(tximport)
library(limma)
library(ggplot2)
library(ggpubr)

# Caminho metadados
file_sledai <- "C:/Users/irisg/Desktop/Tese Lupus/sledai/SLEDAI.csv"

# Caminho resultados 
dir_novo_dataset <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_228066"

# Caminho do geneInfo

load('C:/Users/irisg/Desktop/Tese Lupus/gencode.v38_geneInfo.RData')


# 2. CARREGAR METADADOS CLÍNICOS
# header=FALSE porque o ficheiro começa com dados
meta_clinico <- read.csv(file_sledai, header = FALSE, sep=";", stringsAsFactors = FALSE)
head(meta_clinico)
colnames(meta_clinico) <- c("SampleID", "SLEDAI") # Dar nomes às colunas

# Limpar espaços em branco que possam existir nos nomes
meta_clinico$SampleID <- trimws(meta_clinico$SampleID)

head(meta_clinico)


# 3. IMPORTAR DADOS DE EXPRESSÃO (Kallisto)

folders <- list.dirs(dir_novo_dataset, recursive = FALSE, full.names = FALSE)
files <- file.path(dir_novo_dataset, folders, "abundance.tsv")
names(files) <- folders


txi <- tximport(files_to_use, type = 'kallisto', txOut = F, 
                tx2gene = geneInfo[,1:2], ignoreTxVersion = T, 
                countsFromAbundance = "lengthScaledTPM")

# Converter para Log2TPM
logTPMs_new <- log2(txi$abundance + 1)
# Colocar nomes de genes (Symbols)
rownames(logTPMs_new) <- geneInfo[match(rownames(logTPMs_new), geneInfo$Ensembl_GeneID),"GeneSymbol"]


# 4. LOOP DE CORRELAÇÃO 
# Lista de genes (Interferão + XIST)
target_genes <- c("IFI27", "IFI44", "IFI44L", 
                  "RSAD2", "STAT1", 
                  "ISG15", "IFIT2", 
                  "IRF7", "MX1", "OAS3", "XIST")

pdf("Gene_sledai_correlation.pdf", width=7, height=6)

for (gene in target_genes) {
  
  # Verificar se o gene existe nos dados
  if (gene %in% rownames(logTPMs_new)) {
    
    # Criar tabela para o gráfico
    df_plot <- data.frame(
      SampleID = colnames(logTPMs_new),
      Expression = logTPMs_new[gene, ]
    )
    
    # Juntar com o SLEDAI (pelo SampleID)
    df_plot <- merge(df_plot, meta_clinico, by="SampleID")
    
    # Fazer o Plot
    p <- ggplot(df_plot, aes(x=SLEDAI, y=Expression)) +
      geom_point(color = "darkblue", size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.2) + # Linha de tendência
      theme_classic() +
      labs(title = paste("Gene:", gene, "- Correlation with Disease Severity"),
           x = "SLEDAI Score (Disease Activity)",
           y = paste(gene, "Log2 TPM")) +
      stat_cor(method = "spearman", label.x.npc = "middle") # Correlação de Spearman (melhor para scores)
    
    print(p)
    print(paste("Gráfico gerado para:", gene))
    
  } else {
    print(paste("Aviso: Gene", gene, "não encontrado nos dados."))
  }
}

dev.off()

