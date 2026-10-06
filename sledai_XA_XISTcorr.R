
# X-Ratio vs SLEDAI and XIST 


library(ggplot2)
library(ggpubr)
library(AnnotationDbi)
library(org.Hs.eg.db)


setwd("C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!")


# 1. Carregar logTPMs e metaSuper (para o XIST)
load("Dados_SuperDataset_Prontos.RData")

# 2. Carregar as Raw Counts das PBMCs (GSE228066)
get_counts <- function(caminho) {
  temp_env <- new.env()
  load(caminho, envir = temp_env)
  obj <- ls(temp_env)[1]
  return(temp_env[[obj]]$counts)
}
counts_pbmc <- get_counts("C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_228066/readCounts_gse228066.RData")


# 3. Calcular X/A Ratio COM RAW COUNTS E SEM XIST

print("A calcular Rácios com as Raw Counts...")
chr_data <- mapIds(org.Hs.eg.db, keys = rownames(counts_pbmc), column = "CHR", 
                   keytype = "SYMBOL", multiVals = "first")

genes_X_limpos <- names(chr_data)[!is.na(chr_data) & chr_data == "X" & names(chr_data) != "XIST"]
genes_Auto <- names(chr_data)[!is.na(chr_data) & chr_data %in% as.character(1:22)]

sum_X <- colSums(counts_pbmc[rownames(counts_pbmc) %in% genes_X_limpos, , drop=FALSE])
sum_Auto <- colSums(counts_pbmc[rownames(counts_pbmc) %in% genes_Auto, , drop=FALSE])
racio_pbmc <- sum_X / sum_Auto


# 2. Criar Tabela de Correlação
# Juntar: SLEDAI + X/A Ratio + XIST 
df_final <- data.frame(
  SampleID = colnames(logTPMs_new),
  XA_Ratio = XA_Ratio_New,
  XIST_Expr = logTPMs_new["XIST", ]
)

# Juntar com os dados clínicos 
df_final <- merge(df_final, meta_clinico, by="SampleID")


# 3. Gerar os Gráficos de Correlação

pdf("Analise_Final_Integrada_X_SLEDAI.pdf", width=6, height=12)
par(mfrow=c(3,1)) # Preparar layout vertical

# GRÁFICO 1: XIST vs SLEDAI
p1 <- ggplot(df_final, aes(x=SLEDAI, y=XIST_Expr)) +
  geom_point(color="purple", size=3, alpha=0.6) +
  geom_smooth(method="lm", color="black", alpha=0.2) +
  theme_classic() +
  labs(title="1. XIST vs Disease Severity", 
       y="XIST Expression") +
  stat_cor(method = "spearman")

# GRÁFICO 2: X/A Ratio vs SLEDAI
p2 <- ggplot(df_final, aes(x=SLEDAI, y=XA_Ratio)) +
  geom_point(color="orange", size=3, alpha=0.6) +
  geom_smooth(method="lm", color="black", alpha=0.2) +
  theme_classic() +
  labs(title="2. X-Activity vs Disease Severity", 
       y="X-to-Autosome Ratio") +
  stat_cor(method = "spearman")


# Imprimir num layout bonito
print(ggarrange(p1, p2, ncol=1))

dev.off()

