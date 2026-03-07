#Case: Gene Expression Analysis of Arabidopsis thaliana
#Dataset: GSE7642 (Salt-stressed vs Normal)
#Platform: Microarray (Affymetrix GPL198)
#[ATH1-121501] Affymetrix Arabidopsis ATH1 Genome Array
#Goal: Identify Differentially Expressed Genes (DEG)

#PART A. Work environment preparation
#1. Install BiocManager (manager Bioconductor package)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# 2. Install Bioconductor package (GEOquery & limma)
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)
#Install annotation package based on the platform
BiocManager::install("ath1121501.db", ask = FALSE, update = FALSE)
#3. Install  CRAN package for data manipulation and visualization
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
library(ath1121501.db)
library(AnnotationDbi)
library(umap)

#PART B. Download data from GEO
#Download data in ExpressionSet format with its gene annotation (Gene Symbol) 
gset <- getGEO("GSE7642", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

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
#Check column names available in the data
colnames(pData(gset))
#Check some part of the top row data
head(pData(gset))
#Retrieve sample metadata that contain sample biological condition information
group_info <- pData(gset)[["Treatment length:ch1"]]
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
grup_treatment <- c("X30.minutes","X1.hour","X4.hours","X16.hours","X32.hours")
grup_normal <- "none"
contrast_formula <- paste(grup_treatment, "-", grup_normal)
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
  p.value = 0.05
)
head(topTableResults)

#PART G. Gene name annotation
#Take probe ID from DEG result
probe_ids <- rownames(topTableResults)
#Probe mapping -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  ath1121501.db,
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
#Cek anotation result
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#PART H.1 Expression distribution value (Boxplot)
#Set color based on group
group_colors <- as.numeric(gset$group)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Gene expression distribution value per sample",
  ylab = "Expression Value (log2)"
)
legend(
  "bottomright",
  legend = levels(gset$group),
  fill = unique(group_colors),
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

#PART I.1 Volcano plot visualization
#Volcano plot combine logfc (biological effect) and statistical significance
volcano_data1 <- data.frame(
   logFC = topTableResults$X32.hours...none,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Gene status classification
volcano_data1$status <- "NO"
volcano_data1$status[volcano_data1$logFC > 1 & volcano_data1$adj.P.Val <
                       0.05] <- "UP"
volcano_data1$status[volcano_data1$logFC < -1 & volcano_data1$adj.P.Val <
                       0.05] <- "DOWN"

#Visualization
ggplot(volcano_data1, aes(x = logFC, y = -log10(adj.P.Val), color =
                            status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Salt-stressed 32 hours vs Normal")

volcano_data2 <- data.frame(
  logFC = topTableResults$X16.hours...none,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Gene status classification
volcano_data2$status <- "NO"
volcano_data2$status[volcano_data2$logFC > 1 & volcano_data2$adj.P.Val <
                       0.05] <- "UP"
volcano_data2$status[volcano_data2$logFC < -1 & volcano_data2$adj.P.Val <
                       0.05] <- "DOWN"

#Visualization
ggplot(volcano_data2, aes(x = logFC, y = -log10(adj.P.Val), color =
                            status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Salt-stressed 16 hours vs Normal")

volcano_data3 <- data.frame(
  logFC = topTableResults$X4.hours...none,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Gene status classification
volcano_data3$status <- "NO"
volcano_data3$status[volcano_data3$logFC > 1 & volcano_data3$adj.P.Val <
                       0.05] <- "UP"
volcano_data3$status[volcano_data3$logFC < -1 & volcano_data3$adj.P.Val <
                       0.05] <- "DOWN"

#Visualization
ggplot(volcano_data3, aes(x = logFC, y = -log10(adj.P.Val), color =
                            status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Salt-stressed 4 hours vs Normal")

volcano_data4 <- data.frame(
  logFC = topTableResults$X1.hour...none,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Gene status classification
volcano_data4$status <- "NO"
volcano_data4$status[volcano_data4$logFC > 1 & volcano_data4$adj.P.Val <
                       0.05] <- "UP"
volcano_data4$status[volcano_data4$logFC < -1 & volcano_data4$adj.P.Val <
                       0.05] <- "DOWN"

#Visualization
ggplot(volcano_data4, aes(x = logFC, y = -log10(adj.P.Val), color =
                            status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Salt-stressed 1 hour vs Normal")

volcano_data5 <- data.frame(
  logFC = topTableResults$X30.minutes...none,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Gene status classification
volcano_data5$status <- "NO"
volcano_data5$status[volcano_data5$logFC > 1 & volcano_data5$adj.P.Val <
                       0.05] <- "UP"
volcano_data5$status[volcano_data5$logFC < -1 & volcano_data5$adj.P.Val <
                       0.05] <- "DOWN"

#Visualization
ggplot(volcano_data5, aes(x = logFC, y = -log10(adj.P.Val), color =
                            status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Salt-stressed 30 minutes vs Normal")

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

#Top 50 DEGs for GO term and KEGG Pathway analysis
print(top50$SYMBOL)

#PART J. Save result
#Save result to CSV file
write.csv(topTableResults, "Result_GSE10072_DEG.csv")
message("Analysis is complete. The result file has been saved.")
