# XIST POR TIPO CELULAR (B vs T)

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

df_xist_split <- data.frame(
  Expression = logTPMs["XIST", ],
  Group = metaSuper$Group,
  Batch = metaSuper$Batch,
  CellType = metaSuper$CellType,
  StudyID = metaSuper$Batch # Para separar os gráficos por cell type
)

# 3. Gerar o Gráfico Dividido (Faceted Plot)

p_split <- ggplot(df_xist_split, aes(x=Group, y=Expression, fill=Group)) +
  geom_boxplot(alpha=0.6, outlier.shape = NA) +
  geom_jitter(width=0.2, alpha=0.5) +
  scale_fill_manual(values = c("orange", "#87CEEB","blue","#000080" )) +
  theme_bw() +
  labs(title = "XIST expression by cell type",
       y="XIST Log2 TPM") +
  stat_compare_means(label = "p.format", method = "wilcox.test", 
                     label.y = max(df_xist_split$Expression) + 0.5) +
  facet_wrap(~Batch) # <--- Divide o gráfico em 3

print(p_split)

# 4. Guardar
pdf("XIST_Cell_Type.pdf", width=10, height=5)
print(p_split)
dev.off()

print("Gráfico 'XIST_por_Tipo_Celular.pdf' gerado com sucesso!")
