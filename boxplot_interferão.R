# Interferão
library(limma)
library(openxlsx)
library(ggplot2)
library(ggpubr)

setwd("C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!")
# 1. Carregar Dados (se não tiveres já carregados)
load("Dados_SuperDataset_Prontos.RData")

# 2. Preparar os dados para o gráfico
#  criar uma coluna nova "CellType" baseada no Batch

metaSuper$CellType <- "Unknown"
metaSuper$CellType[metaSuper$Batch =="Batch_PBMC"] <- "PBMC-Cells"
metaSuper$CellType[metaSuper$Batch =="Batch_B"] <- "B-Cells"
metaSuper$CellType[metaSuper$Batch == "Batch_T"] <- "T-Cells"

#VALIDAÇÃO DE GENES ALVO (Interferão & SLEDAI)


# 1. Definir as Listas 
list_interferon <- c("IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", 
                     "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", 
                     "RSAD2", "STAT1", "IRF7")


# BOXPLOTS DO INTERFERÃO DIVIDIDOS POR TIPO CELULAR

pdf("Boxplots_Interferao_CellType.pdf", width = 10, height = 6)


for(gene in list_interferon) {
  
  if (gene %in% rownames(logTPMs)) {
    df_gene <- data.frame(
      Expression = logTPMs[gene, ],
      Group = metaSuper$Group,
      CellType = metaSuper$CellType
    )
    
    p <- ggplot(df_gene, aes(x=Group, y=Expression, fill=Group)) +
      geom_boxplot(outlier.shape = NA, alpha=0.7) +
      geom_jitter(width = 0.2, size=1.5, alpha=0.5) + 
      scale_fill_manual(values = c("orange","#87CEEB", "blue","#000080")) +
      theme_bw() +
      labs(title = paste(gene, "Expression by Cell Type (Interferon List)"),
           y = "Log2 TPM", x = "") +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      facet_wrap(~CellType) # separação
    
    print(p)
  }
}
dev.off()
print("PDF 'Boxplots_Interferao_CellType.pdf' criado!")
