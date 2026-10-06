#META-ANÁLISE (3 DATASETS)

# 1. Instalação de Pacotes
if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!require("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!require("futile.logger", quietly = TRUE)) install.packages("futile.logger")

library(openxlsx)
library(VennDiagram)



# PASSO 1: Caminhos dos Ficheiros


path_PBMC <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_228066/DEG/DEGs_228066.xlsx"
path_bcell <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto_118256/DEG/DEGs_118256.xlsx"
path_tcell   <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!/kalisto109843/DEG/DEGs_109843.xlsx"

# Verificação de segurança
if(!file.exists(path_PBMC)) stop("ERRO: Ficheiro B-Cell 1 não encontrado.")
if(!file.exists(path_bcell)) stop("ERRO: Ficheiro B-Cell 2 não encontrado.")
if(!file.exists(path_tcell))   stop("ERRO: Ficheiro T-Cell não encontrado.")

# PASSO 2: Função de Leitura e Filtragem

get_filtered_degs <- function(filepath) {
  # Lê o ficheiro (que já tem FDR < 0.05)
  degs_table <- read.xlsx(filepath, rowNames = TRUE)
  # Aplica filtro extra de magnitude: |logFC| > 1
  filtered_table <- degs_table[abs(degs_table$logFC) > 1, ]
  return(rownames(filtered_table))
}

# Carregar os dados
degs_PBMC <- get_filtered_degs(path_PBMC)
degs_b <- get_filtered_degs(path_bcell)
degs_t  <- get_filtered_degs(path_tcell)

# Lista organizada para análise
deg_list <- list(
  "PBMC" = degs_PBMC,
  "B-Cells" = degs_b,
  "T-Cells" = degs_t
)


# PASSO 3: diagrama venn

library(VennDiagram)
venn.diagram(
  x = deg_list,
  filename = 'Venn_DEGs_Overlap.png',
  output = TRUE,
  imagetype = "png",
  fill = c("orange", "blue", "gray"), # Cores para os 3 grupos
  col = "transparent",
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  cat.pos = 0,
  margin = 0.1
)


# PASSO 4: Exportar Excel

library(openxlsx)

# Prepara as 3 tabelas
df_PBMC_B <- if(length(common_PBMC_B)>0) data.frame(Genes = common_PBMC_B) else data.frame(Info="Zero genes em comum")
df_PBMC_T <- if(length(common_PBMC_T)>0) data.frame(Genes = common_PBMC_T) else data.frame(Info="Zero genes em comum")
df_All    <- if(length(common_All)>0) data.frame(Genes = common_All) else data.frame(Info="Zero genes em comum")

# A lista dita os separadores (sheets) que o ficheiro Excel vai ter
list_to_save <- list(
  "PBMC_vs_BCells"  = df_PBMC_B,
  "PBMC_vs_TCells"  = df_PBMC_T,
  "Comuns_Aos_Tres" = df_All
)

# Guardar o ficheiro
write.xlsx(list_to_save, file = "Genes_Em_Comum_PBMC_B_T.xlsx")


