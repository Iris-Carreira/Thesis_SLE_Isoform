
#  META-ANÁLISE ESTATÍSTICA (Random Effects Model)


# 1. Instalar pacotes necessários
if (!require("meta", quietly = TRUE)) install.packages("meta")
if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")

library(meta)
library(openxlsx)
library(dplyr)

# PASSO 1: Carregar os Dados COMPLETOS


# Caminhos para os ficheiros com TODOS os genes (edgeR sem filtros)
path_PBMC <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_228066/DEGs_228066.xlsx"
path_bcell <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_118256/DEG/DEGs_118256.xlsx"
path_tcell   <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto109843/DEG/DEGs_109843.xlsx"


df_pbmc <- read.xlsx(path_PBMC, rowNames = TRUE)
df_b    <- read.xlsx(path_bcell, rowNames = TRUE)
df_t    <- read.xlsx(path_tcell,  rowNames = TRUE)


# PASSO 2: Preparar e Cruzar os Dados


prepare_for_meta <- function(df, prefix) {
  df$Gene <- rownames(df)
  
  # Calcular Standard Error (SE) a partir do P-value e LogFC
  df$SE <- abs(df$logFC) / qnorm(1 - df$PValue / 2)
  
  # Limpeza de valores infinitos ou problemáticos
  df <- df[is.finite(df$SE) & is.finite(df$logFC), ]
  
  # Selecionar colunas
  res <- df[, c("Gene", "logFC", "SE")]
  
  # Renomear com o prefixo do estudo (Ex: logFC_PBMC)
  colnames(res) <- c("Gene", paste0("logFC_", prefix), paste0("SE_", prefix))
  return(res)
}

dat_pbmc <- prepare_for_meta(df_pbmc, "PBMC")
dat_b    <- prepare_for_meta(df_b, "B")
dat_t    <- prepare_for_meta(df_t, "T")

# Juntar os 3 datasets (apenas genes presentes em TODOS)
meta_data <- inner_join(dat_pbmc, dat_b, by="Gene") %>%
  inner_join(dat_t, by="Gene")

cat(paste("Total de genes comuns prontos para análise estatística:", nrow(meta_data), "\n"))


# PASSO 3: Random Effects Model

cat("A calcular o modelo de Efeitos Aleatórios... aguarde.\n")

# Vetores para guardar resultados
res_logFC <- numeric(nrow(meta_data))
res_pval  <- numeric(nrow(meta_data))
res_I2    <- numeric(nrow(meta_data))

for(i in 1:nrow(meta_data)) {
  # Extrair os dados do gene i para os 3 estudos
  TE   <- c(meta_data$logFC_PBMC[i], meta_data$logFC_B[i], meta_data$logFC_T[i])
  seTE <- c(meta_data$SE_PBMC[i],    meta_data$SE_B[i],    meta_data$SE_T[i])
  
  # Meta-análise (DerSimonian-Laird)
  m <- metagen(TE, seTE, sm="MD", random=TRUE, fixed=FALSE, verbose=FALSE)
  
  res_logFC[i] <- m$TE.random
  res_pval[i]  <- m$pval.random
  res_I2[i]    <- m$I2
}

# Guardar resultados no dataframe
meta_data$Meta_LogFC <- res_logFC
meta_data$Meta_PValue <- res_pval
meta_data$Heterogeneity_I2 <- res_I2
meta_data$Meta_FDR <- p.adjust(res_pval, method="fdr")

# Filtro Final: |LogFC Combinado| > 1 e FDR Combinado < 0.05
final_genes <- meta_data %>%
  filter(abs(Meta_LogFC) > 1 & Meta_FDR < 0.05) %>%
  arrange(Meta_FDR)

cat(paste("Número de genes significativos após Meta-Análise Estatística:", nrow(final_genes), "\n"))

# Guardar em Excel
write.xlsx(final_genes, "Resultados_Estatisticos_MetaAnalise.xlsx")


