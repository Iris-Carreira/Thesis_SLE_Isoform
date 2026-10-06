

#Most steps will be performed using edgeR R package, which complete tutorial can be found here: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


#########################################################################################################################
# Install packages from R and BioConductor repositories ####
#########################################################################################################################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("tximport","rhdf5", "PCAtools", "edgeR"))

install.packages(c("gplots", "RColorBrewer", "openxlsx"))


#########################################################################################################################
# Import read counts to R ####
# We will import the read counts obtained with Kallisto (Turorial 3) using the tximport R package
#########################################################################################################################

library(tximport)

# First, set the R working directory to the folder containing Kallisto results 
setwd("C:/Users/irisg/Desktop/Tese Lupus/kalisto_118256") #Complete with the path in your computer

# Load the file containing the EnsemblIDs and gene names
load('../gencode.v38_geneInfo.RData');
head(geneInfo)

# Define filenames to import
folders <- dir(pattern="RNA-seq")
fileNames <- file.path(".", folders, 'abundance.tsv');
names(fileNames) <- folders 

# Import all data using tximport function
readCounts_gene <- tximport(fileNames, type = 'kallisto',txOut = F,tx2gene=geneInfo[,1:2], ignoreTxVersion = T,countsFromAbundance = "lengthScaledTPM")

# Explore the object created with tximport
names(readCounts_gene)
View(head(readCounts_gene$abundance)) # matrix containing TPMs
View(head(readCounts_gene$counts)) # matrix containing read counts
head(readCounts_gene$length) 
head(readCounts_gene$countsFromAbundance) # indicates the type of normalization applied to get the values in "abundance"

# To facilitate interpretation of downstream results, replace EnsemblIDs with Official gene symbol (i.s gene name)
geneNames <- geneInfo[match(rownames(readCounts_gene$abundance), geneInfo$Ensembl_GeneID),"GeneSymbol"]
rownames(readCounts_gene$abundance) <- rownames(readCounts_gene$counts) <- geneNames
head(readCounts_gene$abundance)
head(readCounts_gene$counts)


# Save readCounts in an RData file for downstream analyses
save(readCounts_gene, file="readCounts_gse118256.RData")

#########################################################################################################################
# Explore transcriptome profiles using PCA ####
#########################################################################################################################

library(PCAtools)

# Prepare gene expression matrix for PCA analysis: Step 1) get log TPMs (we add 1 TPM to each gene to avoid infinite values after log)
logTPMs <- log2(readCounts_gene$abundance+1)

# Prepare gene expression matrix for PCA analysis: Step 2) remove duplicated genes
uniqueGenes <- unique(rownames(logTPMs))
logTPMs <- logTPMs[uniqueGenes,]


#Prepare metadata with sample type
sampleTypes <- gsub("_RNA-seq", "", colnames(logTPMs))
metaData <- data.frame(sampleTypes); rownames(metaData) <- colnames(logTPMs)
metaData


#Run PCA
pca.res <- pca(logTPMs, metadata=metaData)

#Plot variance explained by each component
screeplot(pca.res)

#Plot 2 selected components/eigenvectors
biplot(pca.res)
biplot(pca.res, colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res, lab="", colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res, x="PC1", y="PC9",lab="",colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3
#biplot para todos??

#Plot several components - fazer para mais de 5??
pairsplot(pca.res, colby="sampleTypes")
pairsplot(pca.res)

# Plot the component loadings and label genes most responsible for variation
# NOTE: Loadings are interpreted as the coefficients of the linear combination of the initial variables from which the principal components are constructed.
# From a statistical point of view, the loadings are equal to the coordinates of the variables divided by the square root of the eigenvalue associated with the component.
PC_genes <- pca.res$loadings
PC1_genes <- PC_genes[order(PC_genes$PC1, decreasing=T),]
head(PC1_genes)
tail(PC1_genes)

pdf("Loadings_plot.pdf") 
plotloadings(pca.res, components = c("PC1", "PC2", "PC3"),rangeRetain =0.1) #retaining 1% of the loadings per PC
#Plot nao funcional, ver o que se passa - so aparece em pdf
dev.off()



# install.packages("gtools")
library(gtools) 

# nomes de coluna LONGOS
long_names <- colnames(logTPMs)

# Extrai o nome do GRUPO ('control' ou 'SLE') para a COR
# remove o número (1, 2, ..., 20) e o sufixo "-RNA-seq"
# Ex: "control1-RNA-seq" -> "control"
sampleTypes <- gsub("[0-9]*[-_]RNA-seq", "", long_names) 

# Coloca a primeira letra em maiúscula (para ficar 'Control' e 'SLE')
sampleTypes <- gsub("control", "Control", sampleTypes)
sampleTypes <- gsub("SLE", "SLE", sampleTypes) 

# Se os nomes dos grupos tiverem variações como 'Contro1',correção extra:
sampleTypes[grep("^Contro", sampleTypes)] <- "Control" 
sampleTypes[grep("^SLE", sampleTypes)] <- "SLE"

# Confirma os grupos 
print(unique(sampleTypes))

# 2. Definir o mapa de cores
colorMap <- c("Control"="orange", "SLE"="blue")

# 3. Cria o vetor de cores final
plotCol <- colorMap[sampleTypes]

# 4. Definir a legenda
legenda_nomes <- names(colorMap)
legenda_cores <- unname(colorMap)

# 5. Cria os rótulos curtos para o Eixo X
# Remove APENAS o sufixo '-RNA-seq' e coloca a primeira letra em maiúscula.
short_labels <- gsub("[-_]RNA-seq", "", long_names)
short_labels <- gsub("control", "Control", short_labels)


# ----------------------------------------------------------------------------------
# Ordenar o vetor de dados, rótulos e cores para o BARPLOT
# ----------------------------------------------------------------------------------

# Usa 'mixedorder' para ordenar as amostras
ordenacao_natural <- gtools::mixedorder(short_labels)

# Ordena os dados e os rótulos
dados_RP11_ordenados <- logTPMs["IGHD", ][ordenacao_natural]
dados_MTTL2_ordenados <- logTPMs["IGHA2", ][ordenacao_natural]
rotulos_ordenados <- short_labels[ordenacao_natural]
cores_ordenadas <- plotCol[ordenacao_natural]


# --- Barplot 1: Gene positivamente correlacionado com PC1 ---
barplot(dados_RP11_ordenados, 
        col=cores_ordenadas, 
        las=2, 
        main="IGHD", 
        ylab="Expression levels (logTPMs)",
        names.arg=rotulos_ordenados) 
legend("topright", fill=legenda_cores, legend=legenda_nomes, bty="n")

# --- Barplot 2: Gene negativamente correlacionado com PC1 ---
barplot(dados_MTTL2_ordenados, 
        col=cores_ordenadas, 
        las=2, 
        main="IGHA2", 
        ylab="Expression levels (logTPMs)",
        names.arg=rotulos_ordenados) 
legend("topright", fill=legenda_cores, legend=legenda_nomes, bty="n")

# You can save some plot in a pdf file
pdf("PCA_plots.pdf") # open pdf file
biplot(pca.res, colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
plotloadings(pca.res, components = c("PC1", "PC2", "PC3"),rangeRetain =0.1) 
dev.off() # closes pdf file


#########################################################################################################################
# Differential expression analysis using edge R package ####
#########################################################################################################################


library(edgeR)

# Create DGEList data class (specific for edgeR package)
y <- DGEList(counts=readCounts_gene$counts, group= sampleTypes)


# Filter out lowly expressed genes
keep <- filterByExpr(y, group=sampleTypes)
y <- y[keep, ,keep.lib.sizes=FALSE]

# Normalize for library sizes
y <- calcNormFactors(y)
y$samples


# Define the design and contrast matrix based on the experimental design, meaning define which comparison to be made
# See more detailed information here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
design_matrix <- model.matrix(~0+sampleTypes)
colnames(design_matrix) <- gsub("sampleTypes", "", colnames(design_matrix))
rownames(design_matrix) <- colnames(logTPMs)
design_matrix
contrast_matrix <- makeContrasts(SLE-Control, levels = design_matrix) # NOTE: REPPLACE HERE WITH THE NAMES OF YOUR SAMPLES
contrast_matrix


# Estimate the dispersion (Biological coefficient of variation) and Fit model
# edgeR uses the negative binomial (NB) distribution to model the read counts for each gene in each sample. The dispersion parameter of the NB distribution accounts for variability between biological replicates (McCarthy, Chen, and Smyth 2012). edgeR estimates an empirical Bayes moderated dispersion for each individual gene. It also estimates a common dispersion, which is a global dispersion estimate averaged over all genes, and a trended dispersion where the dispersion of a gene is predicted from its abundance.
# The vertical axis of the plotBCV plot shows square-root dispersion, also known as biological coefficient of variation (BCV) (McCarthy, Chen, and Smyth 2012).
# For RNA-seq studies, the NB dispersions tend to be higher for genes with very low counts. The dispersion trend tends to decrease smoothly with abundance and to asymptotic to a constant value for genes with larger counts. From our past experience, the asymptotic value for the BCV tends to be in range from 0.05 to 0.2 for genetically identical mice or cell lines, whereas somewhat larger values (>0.3) are observed for human subjects.
y <- estimateDisp(y, design_matrix)
plotBCV(y)
fit <- glmQLFit(y,design_matrix) # Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 


# Differential expression analysis
qlf <- glmQLFTest(fit, contrast = contrast_matrix) # test for differential expression between the experimental groups using quasi-likelihood F-test

# Get differentially expressed genes (DEGs)
print(summary(decideTests(qlf,p.value = 0.05,adjust.method = "fdr"))) 

DEGs <- topTags(qlf, nrow(qlf$table), p.value=0.05, adjust.method = "fdr")$table
upReg <- rownames(DEGs)[which(DEGs$logFC > 0)]
downReg <- rownames(DEGs)[which(DEGs$logFC < 0)]

# Write table with DEGs (Se a tabela estiver vazia, significa que não há DEGs com FDR < 0.05)
library(openxlsx)
write.xlsx(DEGs, file = "DEGs.xlsx", rowNames=T)



# Obter todos os genes com os resultados
allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 0.05)$table
plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)


upReg_visual <- rownames(allGenes)[which(allGenes$logFC > 1 & allGenes$FDR < 0.05)] 
downReg_visual <- rownames(allGenes)[which(allGenes$logFC < -1 & allGenes$FDR < 0.05)]


# Produce Volcano Plot (customized)
plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
abline(h=-log10(0.05), v=c(-1,1),lty=2)

# conjuntos VISUAIS para colorir
points(plotData[upReg_visual,], col="tomato", pch=20)
points(plotData[downReg_visual,], col="steelblue", pch=20)

# Adicionar texto aos 5 principais genes (se existirem)
text(plotData[upReg_visual[1:min(5, length(upReg_visual))],], labels=upReg_visual[1:min(5, length(upReg_visual))], col="black", pos=sample(c(1:3), size=min(5, length(upReg_visual)), replace=T), cex=0.8)
text(plotData[downReg_visual[1:5],], labels=downReg_visual[1:5],col="black", pos=sample(c(1:2), size=min(5, length(downReg_visual)), replace=T), cex=0.8)


# -------------------------------------------------------------------------
# MA Plot (customized)
# -------------------------------------------------------------------------

# Mapeamos o plotData para o MA Plot (LogCPM vs LogFC)
plotData <- cbind(allGenes$logCPM, allGenes$logFC); rownames(plotData) <- rownames(allGenes)
plot(plotData, pch=20,col="gray",xlab="Mean Expression Levels (log CPMs)", ylab="Biological Variation (log2 Fold-Change)")
abline(h=c(-1,1),lty=2)

# Usamos os conjuntos para colorir
points(plotData[upReg_visual,], col="tomato", pch=20)
points(plotData[downReg_visual,], col="steelblue", pch=20)

# Adicionar texto aos 5 principais genes (se existirem)
text(plotData[upReg_visual[1:min(5, length(upReg_visual))],], labels=upReg_visual[1:min(5, length(upReg_visual))], col="black", pos=sample(c(1:3), size=min(5, length(upReg_visual)), replace=T), cex=0.8)
text(plotData[downReg_visual[1:5],], labels=downReg_visual[1:5],col="black", pos=sample(c(1:2), size=min(5, length(downReg_visual)), replace=T), cex=0.8)

# -------------------------------------------------------------------------
# genes para o Heatmap (Top 20)
# -------------------------------------------------------------------------

# quantos genes UP e DOWN. 
n_genes_up <- min(10, length(upReg_visual))
n_genes_down <- min(10, length(downReg_visual))

# Combina os top X genes up-reg e os top X genes down-reg
genes_para_heatmap <- c(
  upReg_visual[1:n_genes_up], 
  downReg_visual[1:n_genes_down]
)

# -------------------------------------------------------------------------
# Preparar as Cores
# -------------------------------------------------------------------------

library(gplots) 
library(RColorBrewer)

# Paleta de Cores para a Coluna Lateral (ColSideColors)
# se plotCol nao definido: plotCol <- colorMap[sampleTypes])

# Paleta de Cores para o CORPO do Heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100) 
#gradiente de Azul (baixo), Branco (média) e Vermelho (alta)

# -------------------------------------------------------------------------
# Heatmap
# -------------------------------------------------------------------------
dev.off()

heatmap.2(
  logTPMs[genes_para_heatmap, ], 
  scale = "row",                      # Standardiza as linhas (mostra desvio da média)
  trace = "none",                     # Remove as linhas tracejadas
  density.info = "none",              # Remove o histograma de densidade
  ColSideColors = plotCol,            # Usa as cores Control/SLE para a barra de cima
  col = heatmap_colors,               # Aplica a paleta Azul/Branco/Vermelho
  dendrogram = "column",              # Agrupa apenas as amostras (colunas)
  margins = c(8, 6),
  labCol = short_labels               # Ajusta margens para nomes legíveis
)


