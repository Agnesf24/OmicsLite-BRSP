#Case: Gene Expression Analysis of Lung Cancer
#Dataset: GSE10072 (Lung Adenocarcinoma vs Normal)
#Platform: Microarray (Affymetrix GPL96 = Affymetrix Human Genome U133A)
#Goal: Identify Differentially Expressed Genes (DEG)

#PART A. Work environment preparation
#1. Install BiocManager (manager Bioconductor package)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# 2. Install Bioconductor package (GEOquery & limma)
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)
#Install annotation package based on the platform
BiocManager::install("hgu133a.db", ask = FALSE, update = FALSE)
#3. Install CRAN package for data manipulation and visualization
install.packages(c("pheatmap", "ggplot2", "dplyr"))
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}
#4. Load library packages
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)

#PART B. Download data from GEO
#Download data in ExpressionSet format with its gene annotation (Gene Symbol) 
gset <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#PART C. Pre-processing expression data
# Retrieve gene expression matrix; Row = probe/gen, Column = sample
ex <- exprs(gset)
#Transform quantile result from named to numerical vector
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
#Logical operator:
#>  : more than
#|| : OR
#&& : AND
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
#If LogTransform = TRUE, apply log2
#Value <= 0, log2 is not applicable, transform to NA
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

#PART D. Define sample group
#Retrieve sample metadata that contain sample biological condition information
group_info <- pData(gset)[["source_name_ch1"]]
#Transform text into format that valid for R
groups <- make.names(group_info)
#Transform categorical data to factor
gset$group <- factor(groups)
#Identify unique category inside factor
group_name <- levels(gset$group)
print(group_name)

#PART E. Matrix design (statistical framework)
#Design matrix for linear model without intercept (best practice limma)
design <- model.matrix(~0 + gset$group)
#Name the column
colnames(design) <- levels(gset$group)
#Decide biological comparison
group_cancer <- "Adenocarcinoma.of.the.Lung"
group_normal <- "Normal.Lung.Tissue"
contrast_formula <- paste(group_cancer, "-", group_normal)
print(paste("Analyzed contrast:", contrast_formula))

#PART F. Differential expression analysis (LIMMA)
#Build linear model for each gene
fit <- lmFit(ex, design)
#Define comparison between groups
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
#Apply contrast to model
fit2 <- contrasts.fit(fit, contrast_matrix)
#Stabilize variants estimation
fit2 <- eBayes(fit2)
#Generate DEG end result
#adjust = "fdr" -> multiple testing correction
#p.value = 0.01 -> very significant gene
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)
head(topTableResults)

#PART G. Gene name annotation
#Take probe ID from DEG result
probe_ids <- rownames(topTableResults)
#Probe mapping -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)
#Combine with limma result
topTableResults$PROBEID <- rownames(topTableResults)
topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)
#Check annotation result
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#PART H.1 Expression distribution value (Boxplot)
#Set color based on group
group_colors <- as.numeric(gset$group)
boxplot(
  ex,
  col = c("darkgreen","hotpink")[group_colors],
  las = 2,
  outline = FALSE,
  main = "Gene expression distribution value per sample",
  ylab = "Expression Value (log2)"
)
legend(
  "topright",
  legend = levels(gset$group),
  fill = c("darkgreen","hotpink")[group_colors],
  cex = 0.8
)

#PART H.2 Expression distribution value (density plot)
#Combine expression & group to data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Gene expression density",
    x = "Expression Value (log2)",
    y = "Density"
)

#PART H.3 UMAP (Low dimension visualization)
#Transpose expression matrix
umap_input <- t(ex)
#Run UMAP
umap_result <- umap(umap_input)
#Save result to data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)
#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
)

#PART I.1 Volcano plot visualization
#Volcano plot combine logfc (biological effect) and statistical significance
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Gene status classification
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val <
                      0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val <
                      0.01] <- "DOWN"
#Visualization
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color =
                           status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Adenocarcinoma vs Normal Lung")

#PART I.2 Heatmap visualization
#Choose 50 significant genes based on adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]
top50 <- head(topTableResults, 50)
#Retrieve expression matrix for the chosen genes
mat_heatmap <- ex[top50$PROBEID, ]
#Use Gene Symbol (fallback to Probe ID)
# if SYMBOL empty → probe ID
# if filled → gene symbol
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID, 
  top50$SYMBOL 
)
rownames(mat_heatmap) <- gene_label
#Data cleaning (required so error hclust not happen)
#Delete row with NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]
#Delete genes with 0 variants
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]
#Column annotation (sample group)
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)
#Heatmap visualization
pheatmap(
  mat_heatmap,
  scale = "row", #Z-score per gene
  annotation_col = annotation_col,
  show_colnames = FALSE, #mute sample name
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)
#To be submitted to GO term and KEGG pathway software
print(top50$SYMBOL)

#PART J. Save result
#Save result to CSV file
write.csv(topTableResults, "Result_GSE10072_DEG.csv")
message("Analysis is complete. The result file has been saved.")
