### R Script for DESeq analysis
setwd("H:/ngs_projects")

library(DESeq2)
library(pheatmap)
library(tidyverse)

#Loading count file and sample info file
data <- read.table("rawcount_input_without_outlier.txt", header = TRUE, row.names = 1)
data <- as.matrix(data)
colnames(data)<- c('HC1', 'HC2',	'HC3',	'HC4',	'HC5',	'HC6',	'HC7',	'LTB1',	'LTB2',	'LTB3',	'LTB4',	'LTB5',	'LTB6',	'LTB7',	'LTB8',	'DS-TB1',	'DS-TB2',	'DS-TB3',	'DS-TB4',	'DS-TB5',	'DS-TB6',	'DR-TB1',	'DR-TB2',	'DR-TB3',	'DR-TB4',	'DR-TB5',	'DR-TB6','DR-TB7')
meta <- read.table("rawcount_info_without_outlier.txt", header = TRUE, row.names = 1)

# Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))
#rownames(meta) <- meta$labels

#Creating DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
dim(dds)

# pre-filtering: removing rows with low gene counts and keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dim(dds)

#Generate the normalized counts and saving for later use
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# set the factor level
dds$condition <- relevel(dds$condition, ref = "HC")


###normalize 
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)


#QC for DE analysis using DESeq2
#Transform normalized counts using the rlog function (vst normalization)
vsd <- vst(dds, blind=FALSE)
write.table(assay(vsd), file="vst_normalized_count.txt", sep="\t", quote=F,col.names=NA)

#Principal components analysis (PCA)
library("ggrepel")
png("pca_plot.png",height=725,width=961)
pca<-plotPCA(vsd, intgroup = "condition") +
  geom_text_repel(aes(label = name), show.legend = FALSE)+
  labs(color = "Group")

print(pca)
dev.off()


#Hierarchical Clustering heatmap
### Extract the rlog matrix from the object
vsd_mat <- assay(vsd)
##Pairwise correlation values for samples
vsd_cor <- cor(vsd_mat)
##Correlation values as a heatmap
png("heatmap.png",height=725,width=1000)
annotation_col <- data.frame(Group = c(rep("HC", 7), rep("LTB", 8), rep("DS-TB", 6), rep("DR-TB", 7)))
row.names(annotation_col) <- colnames(vsd_mat)
pheatmap(vsd_cor,cluster_rows=FALSE,cluster_cols=FALSE,annotation_col=annotation_col,scale = "row")
dev.off()


####heatmap for top 2000  variance genes

var_genes <- apply(vsd_mat, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:2000]
highly_variable_lcpm <- vsd_mat[select_var,]
dim(highly_variable_lcpm)
annotation_col <- data.frame(Group = c(rep("Healthy Controls", 7), rep("Latent Tuberculosis", 8), rep("Drug Sensitive Tuberculosis", 6), rep("Drug Resistant Tuberculosis", 7)))
row.names(annotation_col) <- colnames(vsd_mat)
#png("heatmap_top2000_variance_genes.png",height=725,width=1000)
legend_order <- c("Healthy Controls", "Latent Tuberculosis", "Drug Sensitive Tuberculosis", "Drug Resistant Tuberculosis")
annotation_col$Group <- factor(annotation_col$Group, levels = legend_order)
pheatmap(highly_variable_lcpm,cluster_rows=TRUE,cluster_cols=TRUE,annotation_col=annotation_col,scale = "row",show_rownames=F, color=colorRampPalette(c("green", "black", "red"))(100),annotation_legend = T,annotation_legend_breaks = legend_order)
#dev.off()

#### Dendogram 
png("dendogram.png",height=725,width=961)
hc <- hclust(as.dist(1-cor(vsd_mat, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
sampleTree = as.dendrogram(hc, method="average")
plot(sampleTree)
dev.off()



# Run DESeq 
dds <- DESeq(dds)
#dispersion estimates for our dataset
plotDispEsts(dds)
# Explore Results
res <- results(dds)
summary(res)
# For customizing adjusted p-value
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)



# contrasts
resultsNames(dds)

##volcano plot
library(EnhancedVolcano)
ltb_vs_hc<-results(dds, contrast = c("condition", "LTB","HC"))
png("volcanoplot.png",height=1000,width=1000)
EnhancedVolcano(ltb_vs_hc,
                lab = rownames(ltb_vs_hc),
                x = 'log2FoldChange',
                y = 'pvalue',FCcutoff = 1.5)
dev.off()

res.df <- as.data.frame(ltb_vs_hc)
# writing result foldchange comparison file
write.csv(res.df, file="fold_change_LTB_vs_HC.csv")
#MA Plot
plotMA(res, ylim=c(-2,2))

#heatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=meta, scale = 'row')

