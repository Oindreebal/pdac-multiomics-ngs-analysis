install.packages("BiocManager")
BiocManager::install(c("DESeq2", "biomaRt", "EnhancedVolcano", "pheatmap"))

library(DESeq2)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(EnhancedVolcano)

#########################################
# LOAD DATA
#########################################
counts <- read.csv("D:/BIOINFO/DGE_ad/2_factor/countsdata.csv", row.names = 1)
metadata <- read.csv("D:/BIOINFO/DGE_ad/2_factor/metadata.csv", row.names = 1)

# Ensure counts columns match metadata rownames
all(colnames(counts) %in% rownames(metadata))  # Should return TRUE

# Convert condition and sex to factors
metadata$condition <- factor(metadata$condition)
metadata$sex <- factor(metadata$sex)

dds_condition <- DESeqDataSetFromMatrix(countData = counts,
                                        colData = metadata,
                                        design = ~ condition)

dds_condition <- DESeq(dds_condition)
res_condition <- results(dds_condition, contrast = c("condition", "AD", "control"))

#########################################
# ANNOTATE GENES
#########################################
res_condition$ensembl_gene_id <- gsub("\\..*", "", rownames(res_condition))
res_df_cond <- as.data.frame(res_condition)

# Annotation with biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations_cond <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                          filters = "ensembl_gene_id",
                          values = res_df_cond$ensembl_gene_id,
                          mart = mart)
res_cond_annot <- merge(res_df_cond, annotations_cond, by = "ensembl_gene_id", all.x = TRUE)

#########################################
# RESULT ANALYSIS
#########################################
degs_significant <- res_cond_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

top_up <- degs_significant %>%
  filter(log2FoldChange > 1) %>%
  arrange(padj) %>%
  slice(1:20)

top_down <- degs_significant %>%
  filter(log2FoldChange < -1) %>%
  arrange(padj) %>%
  slice(1:20)

# Save annotated results
write.csv(res_cond_annot, "DEG_AD_vs_Control_annotated.csv", row.names = FALSE)
write.csv(degs_significant, "DEG_AD_vs_Control_significant.csv", row.names = FALSE)
write.csv(top_up, "Top_Upregulated_AD_vs_Control.csv", row.names = FALSE)
write.csv(top_down, "Top_Downregulated_AD_vs_Control.csv", row.names = FALSE)

#########################################
# VISUALISATION
#########################################

# Volcano Plot
res_cond_annot$significance <- "Not Significant"
res_cond_annot$significance[res_cond_annot$log2FoldChange > 1 & res_cond_annot$padj < 0.05] <- "Upregulated"
res_cond_annot$significance[res_cond_annot$log2FoldChange < -1 & res_cond_annot$padj < 0.05] <- "Downregulated"

ggplot(res_cond_annot, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot: AD vs Control", x = "log2 Fold Change", y = "-log10 Adjusted P-Value")

dds_sex <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ sex)

dds_sex <- DESeq(dds_sex)
res_sex <- results(dds_sex, contrast = c("sex", "male", "female"))

res_sex$ensembl_gene_id <- gsub("\\..*", "", rownames(res_sex))
res_df_sex <- as.data.frame(res_sex)

# Annotation with biomaRt
annotations_sex <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                         filters = "ensembl_gene_id",
                         values = res_df_sex$ensembl_gene_id,
                         mart = mart)
res_sex_annot <- merge(res_df_sex, annotations_sex, by = "ensembl_gene_id", all.x = TRUE)

#Result analysis
degs_significant_sex <- res_sex_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

top_up_sex <- degs_significant_sex %>%
  filter(log2FoldChange > 1) %>%
  arrange(padj) %>%
  slice(1:20)

top_down_sex <- degs_significant_sex %>%
  filter(log2FoldChange < -1) %>%
  arrange(padj) %>%
  slice(1:20)

#Save Results
write.csv(res_sex_annot, "DEG_Male_vs_Female_annotated.csv", row.names = FALSE)
write.csv(degs_significant_sex, "DEG_Male_vs_Female_significant.csv", row.names = FALSE)
write.csv(top_up_sex, "Top_Upregulated_Male_vs_Female.csv", row.names = FALSE)
write.csv(top_down_sex, "Top_Downregulated_Male_vs_Female.csv", row.names = FALSE)

# Volcano Plot
res_sex_annot$significance <- "Not Significant"
res_sex_annot$significance[res_sex_annot$log2FoldChange > 1 & res_sex_annot$padj < 0.05] <- "Upregulated"
res_sex_annot$significance[res_sex_annot$log2FoldChange < -1 & res_sex_annot$padj < 0.05] <- "Downregulated"

ggplot(res_sex_annot, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot: Male vs Female", x = "log2 Fold Change", y = "-log10 Adjusted P-Value")

# Make sure to strip Ensembl version suffixes early
rownames(counts) <- gsub("\\..*", "", rownames(counts))

# Rebuild DESeqDataSet with updated rownames
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition + sex)

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Transform counts (for Plotting)
vsd <- vst(dds, blind = FALSE)
# Update vst rownames too
rownames(vsd) <- gsub("\\..*", "", rownames(vsd))

# Ensure sample order matches
stopifnot(all(colnames(vsd) == rownames(metadata)))

#PCA
# PCA Plot: AD vs Control
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA Plot: AD vs Control")
print(pca_plot)  # Show in RStudio

# PCA Plot: Male vs Female
pca_plot2 <- plotPCA(vsd, intgroup = "sex") +
  ggtitle("PCA Plot: Male vs Female")
print(pca_plot2)  # Show in RStudio

#heatmap of top expressed genes
top50 <- head(order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE), 50)

pheatmap(assay(vsd)[top50, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = metadata,
         show_rownames = FALSE,
         main = "Heatmap of Most Expressed Genes")

# heatmap AD vs Control
top_degs_cond <- res_cond_annot %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice(1:50)

deg_ids_cond <- intersect(top_degs_cond$ensembl_gene_id, rownames(vsd))

# Plot
pheatmap(assay(vsd)[deg_ids_cond, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = metadata,
         show_rownames = FALSE,
         main = "Top DEGs: AD vs Control")

#heatmap male vs female
top_degs_cond <- res_sex_annot %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice(1:50)

deg_ids_cond <- intersect(top_degs_cond$ensembl_gene_id, rownames(vsd))

pheatmap(assay(vsd)[deg_ids_cond, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = metadata,
         show_rownames = FALSE,
         main = "Top DEGs: Male vs Female")

#enhanced volcano plots
EnhancedVolcano(res_cond_annot,
                lab = res_cond_annot$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Enhanced Volcano: AD vs Control',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3.5)

EnhancedVolcano(res_sex_annot,
                lab = res_sex_annot$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Enhanced Volcano: Male vs Female',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3.5)

#gene expression bar plot control vs AD
top20_ad <- res_cond_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice(1:20)

expr_vals_ad <- rowMeans(assay(vsd)[top20_ad$ensembl_gene_id, ])
bar_df_ad <- data.frame(
  gene = ifelse(is.na(top20_ad$external_gene_name) | top20_ad$external_gene_name == "",
                top20_ad$ensembl_gene_id,
                top20_ad$external_gene_name),
  mean_expr = expr_vals_ad
)

bar_plot_ad <- ggplot(bar_df_ad, aes(x = reorder(gene, -mean_expr), y = mean_expr, fill = gene)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Mean Expression of Top DEGs: AD vs Control",
       x = "Gene", y = "Mean VST Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
print(bar_plot_ad)

#gene expression bar plot male vs female
top20_cond <- res_sex_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice(1:20)

# Get mean VST expression
expr_vals <- rowMeans(assay(vsd)[top20_cond$ensembl_gene_id, ])
bar_df <- data.frame(
  gene = ifelse(is.na(top20_cond$external_gene_name) | top20_cond$external_gene_name == "",
                top20_cond$ensembl_gene_id,
                top20_cond$external_gene_name),
  mean_expr = expr_vals
)

# Plot
ggplot(bar_df, aes(x = reorder(gene, -mean_expr), y = mean_expr, fill = gene)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Mean Expression of Top 20 DEGs: Male vs Female",  
       x = "Gene", y = "Mean VST Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
