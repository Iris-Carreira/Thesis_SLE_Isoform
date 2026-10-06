
#X/A RATIO E CORRELAÇÃO XIST x X/A, COM CONTAGENS BRUTAS para dataset e subtipos


library(tximport)
library(ggplot2)
library(ggpubr)
library(AnnotationDbi)
library(org.Hs.eg.db)

#paleta
paleta_cores <- c(
  "Control" = "orange",  # laranja
  "SLE"     = "blue",  # azul
  "mSLE"    = "skyblue",  # azul claro
  "sSLE"    = "navy",  # azul escuro
  # subtipos de B-cells e T-cells: tons de azul/roxo diferentes
  "rN"      = "slateblue",
  "T3"      = "deepskyblue",
  "aN"      = "cyan4",
  "DN2"     = "royalblue4",
  "SM"      = "lightskyblue",
  "ASC"     = "midnightblue",
  "Tfh"     = "turquoise",
  "PD1POS"  = "cornflowerblue"
)


#CAMINHOS:


base_dir <- "C:/Users/irisg/Desktop/Tese Lupus/Coisas organizadas!!"

dir_B      <- file.path("C:/Users/irisg/Desktop/kalisto_118256")
meta_B_csv <- file.path(base_dir, "kalisto_118256/metadados_118256_subtypes.csv")

dir_T      <- file.path(base_dir, "kalisto109843")
meta_T_csv <- file.path(dir_T, "metadados_109846_linhascelulares.csv")

dir_PBMC      <- file.path(base_dir, "kalisto_228066")
meta_PBMC_csv <- file.path(dir_PBMC, "metadados_228066.csv")

load(file.path(base_dir, "..", "gencode.v38_geneInfo.RData"))


dir_output <- getwd()  # onde os gráficos e a tabela final vão ser guardados

# 1. IMPORTAÇÃO — devolve as contagens brutas (txi$counts)
# -------------------------------------------------------------------------
importar_dataset <- function(dir_dados, meta_csv, nome_dataset, formato = "padrao") {
  
  if (formato == "padrao") {
    # T-cells, PBMC: csv com colunas SampleID;Group, pasta = SampleID
    meta <- read.csv(meta_csv, header = TRUE, sep = ";", stringsAsFactors = FALSE)
    colnames(meta) <- c("SampleID", "Group")
    meta$SampleID <- trimws(meta$SampleID)
    meta$Group    <- trimws(meta$Group)
    
  } else if (formato == "b_cells_srr") {
  
    meta <- read.csv(meta_csv, header = TRUE, sep = ";", stringsAsFactors = FALSE)
    colnames(meta) <- c("SampleID", "Group")
    meta$SampleID <- trimws(meta$SampleID)
    meta$Group    <- trimws(meta$Group)
    
    # extrai o código curto entre parêntesis, ex: "Activated Naive B cells (aN)" -> "aN"
    # se não houver parêntesis (caso do "Control"), mantém o valor como está
    tem_codigo <- grepl("\\(([^)]+)\\)", meta$Group)
    meta$Group[tem_codigo] <- sub(".*\\(([^)]+)\\).*", "\\1", meta$Group[tem_codigo])
  }
  
  meta$Disease <- ifelse(meta$Group == "Control", "Control", "SLE")
  
  folders <- list.dirs(dir_dados, full.names = FALSE, recursive = FALSE)
  
  # algumas pastas (B-cells) têm sufixo "_quant" que não está no csv de
  # metadados (ex: "SRR14530639_quant" vs "SRR14530639") — remover esse
  # sufixo só para efeitos de correspondência.
  folders_clean <- sub("_quant$", "", folders)
  match_idx <- folders_clean %in% meta$SampleID
  folders <- folders[match_idx]
  folders_clean <- folders_clean[match_idx]
  
  files <- file.path(dir_dados, folders, "abundance.tsv")
  names(files) <- folders_clean   # as colunas da matriz ficam com "SRR..." limpo
  files <- files[file.exists(files)]
  
  #debug
  if (length(files) == 0) stop(paste("ERRO:", nome_dataset, "- nenhum ficheiro encontrado. Verifica os caminhos!"))
  cat("A importar", length(files), "amostras de", nome_dataset, "...\n")
  
  txi <- tximport(files,
                  type = "kallisto",
                  txOut = FALSE,
                  tx2gene = geneInfo[, 1:2],
                  ignoreTxVersion = TRUE,
                  countsFromAbundance = "no")
  
  raw_counts <- txi$counts
  gene_symbols <- geneInfo$GeneSymbol[match(rownames(raw_counts), geneInfo$Ensembl_GeneID)]
  valid <- !is.na(gene_symbols) & !duplicated(gene_symbols)
  raw_counts <- raw_counts[valid, ]
  rownames(raw_counts) <- gene_symbols[valid]
  
  logTPM <- log2(txi$abundance + 1)
  logTPM <- logTPM[valid, ]
  rownames(logTPM) <- gene_symbols[valid]
  
  meta_aligned <- meta[match(colnames(raw_counts), meta$SampleID), ]
  rownames(meta_aligned) <- meta_aligned$SampleID
  
  list(raw_counts = raw_counts, logTPM = logTPM, meta = meta_aligned)
}



dados_B    <- importar_dataset(dir_B, meta_B_csv, "B-cells", formato = "b_cells_srr")
dados_T    <- importar_dataset(dir_T, meta_T_csv, "T-cells", formato = "padrao")
dados_PBMC <- importar_dataset(dir_PBMC, meta_PBMC_csv, "PBMC", formato = "padrao")

# 2. X/A RATIO COM CONTAGENS BRUTAS — soma de contagens de genes do X dividida pela soma de contagens de genes do autossoma
calc_XA_ratio_raw <- function(raw_counts_matrix) {
  chr_data <- mapIds(org.Hs.eg.db,
                     keys = rownames(raw_counts_matrix),
                     column = "CHR",
                     keytype = "SYMBOL",
                     multiVals = "first")
  genes_X    <- names(chr_data)[!is.na(chr_data) & chr_data == "X"]
  genes_Auto <- names(chr_data)[!is.na(chr_data) & chr_data %in% as.character(1:22)]
  
  sum_X    <- colSums(raw_counts_matrix[rownames(raw_counts_matrix) %in% genes_X, ])
  sum_Auto <- colSums(raw_counts_matrix[rownames(raw_counts_matrix) %in% genes_Auto, ])
  sum_X / sum_Auto
}

XA_Ratio_B    <- calc_XA_ratio_raw(dados_B$raw_counts)
XA_Ratio_T    <- calc_XA_ratio_raw(dados_T$raw_counts)
XA_Ratio_PBMC <- calc_XA_ratio_raw(dados_PBMC$raw_counts)



# 3. GRÁFICOS DE X/A RATIO POR CONDIÇÃO

plot_XA_por_grupo <- function(XA_Ratio, meta_aligned, nome_dataset, arquivo) {
  df <- data.frame(XA_Ratio = XA_Ratio, Group = meta_aligned$Group)
  p <- ggplot(df, aes(x = Group, y = XA_Ratio, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.6) +
    stat_compare_means() +
    scale_fill_manual(values = paleta_cores) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none") +
    labs(title = paste("X/A Ratio —", nome_dataset),
         x = "Group", y = "X-to-Autosome Ratio")
  ggsave(file.path(dir_output, arquivo), p, width = 7, height = 5)
  print(p)
}

plot_XA_por_grupo(XA_Ratio_B, dados_B$meta, "B-cells", "XA_ratio_rawcounts_Bcells.pdf")
plot_XA_por_grupo(XA_Ratio_T, dados_T$meta, "T-cells", "XA_ratio_rawcounts_Tcells.pdf")
plot_XA_por_grupo(XA_Ratio_PBMC, dados_PBMC$meta, "PBMC", "XA_ratio_rawcounts_PBMC.pdf")


# 4. CORRELAÇÃO XIST x X/A RATIO 

resultados_familia1 <- data.frame(
  Dataset = character(), Comparison = character(),
  R = numeric(), p_raw = numeric(), n = integer()
)

add_result <- function(df, dataset, comparison, xist, xa) {
  ok <- complete.cases(xist, xa) & is.finite(xist) & is.finite(xa)
  teste <- cor.test(xist[ok], xa[ok], method = "spearman")
  rbind(df, data.frame(
    Dataset = dataset, Comparison = comparison,
    R = unname(teste$estimate), p_raw = teste$p.value, n = sum(ok)
  ))
}

plot_correlacao <- function(xist, xa, group_labels, titulo, arquivo) {
  df <- data.frame(XIST = xist, XA_Ratio = xa, Group = group_labels)
  p <- ggplot(df, aes(x = XIST, y = XA_Ratio, color = Group)) +
    geom_point(size = 2.5, alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
    stat_cor(method = "spearman", size = 3.2) +
    scale_color_manual(values = paleta_cores) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(title = titulo, x = "XIST Log2 TPM", y = "X/A Ratio (raw counts)")
  ggsave(file.path(dir_output, arquivo), p, width = 6, height = 5)
  print(p)
}

#B-cells: overall, por grupo, por subtipo
xist_B <- as.numeric(dados_B$logTPM["XIST", ])

resultados_familia1 <- add_result(resultados_familia1, "B-cells", "Overall", xist_B, XA_Ratio_B)
plot_correlacao(xist_B, XA_Ratio_B, dados_B$meta$Disease,
                "XIST vs X/A - B-cells (overall)", "Corr_XIST_XA_Bcells_overall.pdf")

for (grp in c("Control", "SLE")) {
  idx <- dados_B$meta$Disease == grp
  resultados_familia1 <- add_result(resultados_familia1, "B-cells", paste0("Group: ", grp),
                                    xist_B[idx], XA_Ratio_B[idx])
}

subtypes_order_B <- c("rN", "T3", "aN", "DN2", "SM", "ASC")
for (sub in subtypes_order_B) {
  idx <- dados_B$meta$Group %in% c("Control", sub)
  resultados_familia1 <- add_result(resultados_familia1, "B-cells", paste0("Subtype: ", sub),
                                    xist_B[idx], XA_Ratio_B[idx])
  plot_correlacao(xist_B[idx], XA_Ratio_B[idx],
                  factor(dados_B$meta$Group[idx], levels = c("Control", sub)),
                  paste("XIST vs X/A - B-cells,", sub, "vs Control"),
                  paste0("Corr_XIST_XA_Bcells_", sub, ".pdf"))
}

# T-cells: overall, por grupo, por subtipo
xist_T <- as.numeric(dados_T$logTPM["XIST", ])

resultados_familia1 <- add_result(resultados_familia1, "T-cells", "Overall", xist_T, XA_Ratio_T)
plot_correlacao(xist_T, XA_Ratio_T, dados_T$meta$Disease,
                "XIST vs X/A - T-cells (overall)", "Corr_XIST_XA_Tcells_overall.pdf")

for (grp in c("Control", "SLE")) {
  idx <- dados_T$meta$Disease == grp
  resultados_familia1 <- add_result(resultados_familia1, "T-cells", paste0("Group: ", grp),
                                    xist_T[idx], XA_Ratio_T[idx])
}

subtypes_order_T <- c("PD1POS", "Tfh")
for (sub in subtypes_order_T) {
  idx <- dados_T$meta$Group %in% c("Control", sub)
  resultados_familia1 <- add_result(resultados_familia1, "T-cells", paste0("Subtype: ", sub),
                                    xist_T[idx], XA_Ratio_T[idx])
  plot_correlacao(xist_T[idx], XA_Ratio_T[idx],
                  factor(dados_T$meta$Group[idx], levels = c("Control", sub)),
                  paste("XIST vs X/A - T-cells,", sub, "vs Control"),
                  paste0("Corr_XIST_XA_Tcells_", sub, ".pdf"))
}

# PBMC: overall, por grupo 
xist_PBMC <- as.numeric(dados_PBMC$logTPM["XIST", ])

resultados_familia1 <- add_result(resultados_familia1, "PBMC", "Overall", xist_PBMC, XA_Ratio_PBMC)
plot_correlacao(xist_PBMC, XA_Ratio_PBMC, dados_PBMC$meta$Group,
                "XIST vs X/A - PBMC (overall)", "Corr_XIST_XA_PBMC_overall.pdf")

grupos_pbmc <- unique(dados_PBMC$meta$Group)


for (grp in grupos_pbmc[grupos_pbmc != "Control"]) {
 
  idx <- dados_PBMC$meta$Group == grp
  resultados_familia1 <- add_result(resultados_familia1, "PBMC", paste0("Group: ", grp),
                                    xist_PBMC[idx], XA_Ratio_PBMC[idx])
}


# 5. AJUSTE CONJUNTO


resultados_familia1$p_adj <- p.adjust(resultados_familia1$p_raw, method = "BH")

write.csv(resultados_familia1, file.path(dir_output, "Familia1_XIST_XA_pvalues_ajustados_rawcounts.csv"), row.names = FALSE)


# 6. GRÁFICO COMBINADO — os 3 datasets lado a lado no mesmo painel,


df_combinado <- rbind(
  data.frame(XIST = xist_B, XA_Ratio = XA_Ratio_B,
             Group = dados_B$meta$Disease, Dataset = "B-Cells"),
  data.frame(XIST = xist_T, XA_Ratio = XA_Ratio_T,
             Group = dados_T$meta$Disease, Dataset = "T-Cells"),
  data.frame(XIST = xist_PBMC, XA_Ratio = XA_Ratio_PBMC,
             Group = dados_PBMC$meta$Group, Dataset = "PBMC-Cells")
)

# ordem dos painéis B-Cells, PBMC-Cells, T-Cells
df_combinado$Dataset <- factor(df_combinado$Dataset,
                               levels = c("B-Cells", "PBMC-Cells", "T-Cells"))

p_combinado <- ggplot(df_combinado, aes(x = XIST, y = XA_Ratio)) +
  geom_point(aes(color = Group), size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, alpha = 0.15) +
  stat_cor(method = "spearman", size = 3.5, label.x.npc = "left", label.y.npc = "top") +
  scale_color_manual(values = paleta_cores) +
  facet_wrap(~ Dataset, scales = "free", nrow = 1) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.title = element_blank(),
        strip.background = element_rect(fill = "grey85")) +
  labs(title = "Correlation: XIST vs X-Chromosome Activity (X/A)",
       x = "XIST Expression (Log2 TPM)", y = "X-to-Autosome Ratio")

ggsave(file.path(dir_output, "Correlacao_XIST_XA_combinado_3datasets.pdf"),
       p_combinado, width = 14, height = 5)
print(p_combinado)



# 6b. GRÁFICO COMBINADO, GRUPOS SEPARADOS DENTRO DE CADA PAINEL — mesmo facet_wrap dos 3 datasets, com uma reta e um R/p (Spearman) por grupo dentro de cada painel


df_combinado_grupos <- rbind(
  data.frame(XIST = xist_B, XA_Ratio = XA_Ratio_B,
             Group = dados_B$meta$Disease, Dataset = "B-Cells"),
  data.frame(XIST = xist_T, XA_Ratio = XA_Ratio_T,
             Group = dados_T$meta$Disease, Dataset = "T-Cells"),
  data.frame(XIST = xist_PBMC, XA_Ratio = XA_Ratio_PBMC,
             Group = dados_PBMC$meta$Group, Dataset = "PBMC-Cells")
)


df_combinado_grupos$Dataset <- factor(df_combinado_grupos$Dataset,
                                      levels = c("B-Cells", "PBMC-Cells", "T-Cells"))

p_combinado_grupos <- ggplot(df_combinado_grupos, aes(x = XIST, y = XA_Ratio, color = Group)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  stat_cor(method = "spearman", size = 3.2, label.x.npc = "left",
           show.legend = FALSE,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "*`, `*"))) +
  scale_color_manual(values = paleta_cores) +
  facet_wrap(~ Dataset, scales = "free", nrow = 1) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.title = element_blank(),
        strip.background = element_rect(fill = "grey85")) +
  labs(title = "Correlation: XIST vs X-Chromosome Activity (X/A)",
       x = "XIST Expression (Log2 TPM)", y = "X-to-Autosome Ratio")

ggsave(file.path(dir_output, "Correlacao_XIST_XA_combinado_3datasets_porgrupo.pdf"),
       p_combinado_grupos, width = 14, height = 5)
print(p_combinado_grupos)



# 7. FACET_WRAP POR SUBTIPO — uma figura só por dataset Cada painel = Control vs. um subtipo

#  B-cells: 6 subtipos num facet_wrap só 
df_subtipos_B <- do.call(rbind, lapply(subtypes_order_B, function(sub) {
  idx <- dados_B$meta$Group %in% c("Control", sub)
  data.frame(
    XIST = xist_B[idx], XA_Ratio = XA_Ratio_B[idx],
    Group = factor(dados_B$meta$Group[idx], levels = c("Control", sub)),
    Subtype = sub
  )
}))
df_subtipos_B$Subtype <- factor(df_subtipos_B$Subtype, levels = subtypes_order_B)

p_subtipos_B <- ggplot(df_subtipos_B, aes(x = XIST, y = XA_Ratio, color = Group)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  stat_cor(method = "spearman", size = 3, label.x.npc = "left",
           show.legend = FALSE,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "*`, `*"))) +
  scale_color_manual(values = paleta_cores) +
  facet_wrap(~ Subtype, scales = "free", nrow = 2) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "grey85")) +
  labs(title = "XIST vs X/A Ratio — B-cell subtypes (Control vs. each subtype)",
       x = "XIST Expression (Log2 TPM)", y = "X-to-Autosome Ratio")

ggsave(file.path(dir_output, "Correlacao_XIST_XA_Bcells_subtipos_facet.pdf"),
       p_subtipos_B, width = 12, height = 8)
print(p_subtipos_B)


#  T-cells: 2 subtipos num facet_wrap só 
df_subtipos_T <- do.call(rbind, lapply(subtypes_order_T, function(sub) {
  idx <- dados_T$meta$Group %in% c("Control", sub)
  data.frame(
    XIST = xist_T[idx], XA_Ratio = XA_Ratio_T[idx],
    Group = factor(dados_T$meta$Group[idx], levels = c("Control", sub)),
    Subtype = sub
  )
}))
df_subtipos_T$Subtype <- factor(df_subtipos_T$Subtype, levels = subtypes_order_T)

p_subtipos_T <- ggplot(df_subtipos_T, aes(x = XIST, y = XA_Ratio, color = Group)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  stat_cor(method = "spearman", size = 3.2, label.x.npc = "left",
           show.legend = FALSE,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "*`, `*"))) +
  scale_color_manual(values = paleta_cores) +
  facet_wrap(~ Subtype, scales = "free", nrow = 1) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "grey85")) +
  labs(title = "XIST vs X/A Ratio — T-cell subtypes (Control vs. each subtype)",
       x = "XIST Expression (Log2 TPM)", y = "X-to-Autosome Ratio")

ggsave(file.path(dir_output, "Correlacao_XIST_XA_Tcells_subtipos_facet.pdf"),
       p_subtipos_T, width = 10, height = 5)
print(p_subtipos_T)


# 8.PAINEL DE 3 GRÁFICOS POR SUBTIPO (XIST expression, X/A ratio, e a correlação entre os dois)

gerar_painel_subtipo <- function(xist, xa, group_labels, sub, dataset_label) {
  df <- data.frame(XIST = xist, XA_Ratio = xa, Group = group_labels)
  
  # Painel 1: XIST expression por grupo
  p_xist <- ggplot(df, aes(x = Group, y = XIST, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
    stat_compare_means(size = 3) +
    scale_fill_manual(values = paleta_cores) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none") +
    labs(title = "XIST Expression", x = NULL, y = "XIST (Log2 TPM)")
  
  # Painel 2: X/A ratio por grupo
  p_xa <- ggplot(df, aes(x = Group, y = XA_Ratio, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
    stat_compare_means(size = 3) +
    scale_fill_manual(values = paleta_cores) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none") +
    labs(title = "X/A Ratio", x = NULL, y = "X-to-Autosome Ratio")
  
  # Painel 3: correlação XIST vs X/A, com R/p (Spearman) por grupo
  p_corr <- ggplot(df, aes(x = XIST, y = XA_Ratio, color = Group)) +
    geom_point(size = 2, alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    stat_cor(method = "spearman", size = 2.8, label.x.npc = "left",
             show.legend = FALSE,
             aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "*`, `*"))) +
    scale_color_manual(values = paleta_cores) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(title = "XIST vs X/A Correlation", x = "XIST (Log2 TPM)", y = "X-to-Autosome Ratio")
  
  painel <- ggarrange(p_xist, p_xa, p_corr, ncol = 3, nrow = 1, widths = c(1, 1, 1.3))
  painel <- annotate_figure(painel,
                            top = text_grob(paste0(dataset_label, " — Control vs. ", sub),
                                            face = "bold", size = 13))
  painel
}

pdf(file.path(dir_output, "Annex_XIST_XA_subtipos.pdf"), width = 13, height = 5)

for (sub in subtypes_order_B) {
  idx <- dados_B$meta$Group %in% c("Control", sub)
  painel <- gerar_painel_subtipo(
    xist_B[idx], XA_Ratio_B[idx],
    factor(dados_B$meta$Group[idx], levels = c("Control", sub)),
    sub, "B-cells"
  )
  print(painel)
}

for (sub in subtypes_order_T) {
  idx <- dados_T$meta$Group %in% c("Control", sub)
  painel <- gerar_painel_subtipo(
    xist_T[idx], XA_Ratio_T[idx],
    factor(dados_T$meta$Group[idx], levels = c("Control", sub)),
    sub, "T-cells"
  )
  print(painel)
}

dev.off()

