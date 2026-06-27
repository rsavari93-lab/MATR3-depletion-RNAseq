## =============================================================================
## Figure 3 — Quality assessment and reproducibility of bulk RNA-seq datasets
## Panels:
##   a) PCA of VST-transformed expression data (PC1 vs PC2)
##   b) Pairwise Pearson correlation heatmap across all 12 samples
## Output: single combined figure (Fig3_PCA_Correlation_combined.tiff / .pdf)
library(DESeq2)
library(ggplot2)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(ggplotify)
library(grid)

setwd("D:/")

## ---- 0b. Journal-standard font settings (same constants used in Fig 4) -----
FIG_FONT          <- "sans"
TITLE_SIZE        <- 10
TITLE_FACE        <- "bold"
AXIS_TITLE_SIZE   <- 9
AXIS_TEXT_SIZE    <- 8
LEGEND_TITLE_SIZE <- 8
LEGEND_TEXT_SIZE  <- 7.5
TAG_SIZE          <- 11

## ---- 1. Panel a — PCA of VST-transformed data -------------------------------
count_data <- read_excel("GENE.COUNT.xlsm", sheet = "gene_count")
count_matrix <- as.matrix(count_data[, 2:13])
rownames(count_matrix) <- count_data$gene_id

sample_names <- colnames(count_matrix)
col_data <- data.frame(
  row.names = sample_names,
  cellLine  = ifelse(grepl("U87", sample_names), "U87", "SH-SY5Y"),
  condition = ifelse(grepl("siMATR3", sample_names), "KD", "Control")
)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = col_data,
  design    = ~ cellLine + condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
vsd <- vst(dds, blind = TRUE)
vst_data <- assay(vsd)

gene_var <- apply(vst_data, 1, var)
vst_data_filtered <- vst_data[gene_var > 0, ]

pca <- prcomp(t(vst_data_filtered), center = TRUE, scale. = TRUE)
var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)

pca_data <- data.frame(
  Sample    = rownames(pca$x),
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  CellLine  = col_data$cellLine,
  Condition = col_data$condition
)

pA <- ggplot(pca_data, aes(x = PC1, y = PC2, color = CellLine, shape = Condition)) +
  geom_point(size = 2.8, stroke = 0.75) +
  scale_color_manual(
    values = c("SH-SY5Y" = "#1a6faf", "U87" = "#f5c800"),
    name = "Cell line"
  ) +
  scale_shape_manual(
    values = c("Control" = 16, "KD" = 17),
    name = "Condition"
  ) +
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  theme_classic(base_family = FIG_FONT) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    axis.title    = element_text(size = AXIS_TITLE_SIZE, family = FIG_FONT),
    axis.text     = element_text(size = AXIS_TEXT_SIZE, family = FIG_FONT),
    legend.title  = element_text(size = LEGEND_TITLE_SIZE, family = FIG_FONT),
    legend.text   = element_text(size = LEGEND_TEXT_SIZE, family = FIG_FONT),
    legend.position = "right",
    legend.background = element_rect(color = "gray80", linewidth = 0.4),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_blank()
  )

## ---- 2. Panel b — Pairwise Pearson correlation heatmap ----------------------
corr_data <- read_excel("CORRELATION.xlsx", col_names = TRUE)
corr_data <- as.data.frame(corr_data)
rownames(corr_data) <- corr_data[, 1]
corr_data <- corr_data[, -1]
corr_data <- as.matrix(corr_data)

annotation <- data.frame(
  Cell  = ifelse(grepl("U87", colnames(corr_data)), "U87", "SH-SY5Y"),
  Group = ifelse(grepl("siMATR3", colnames(corr_data)), "KD", "Control")
)
rownames(annotation) <- colnames(corr_data)

ann_colors <- list(
  Cell  = c("SH-SY5Y" = "#1a6faf", "U87" = "#f5c800"),
  Group = c("Control" = "#888888", "KD" = "#E84040")
)

pB_grob <- pheatmap(
  corr_data,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(0.6, 1, length.out = 101),
  annotation_col = annotation,
  annotation_row = annotation,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = AXIS_TEXT_SIZE - 2,   # slightly smaller than axis text so numbers fit in cells
  fontsize         = AXIS_TEXT_SIZE,
  fontsize_row     = AXIS_TEXT_SIZE,
  fontsize_col     = AXIS_TEXT_SIZE,
  fontfamily       = FIG_FONT,
  angle_col = 45,
  border_color = "grey70",
  main = "",
  legend = TRUE,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  legend_breaks = c(0.6, 0.7, 0.8, 0.9, 1.0),
  legend_labels = c("0.6", "0.7", "0.8", "0.9", "1.0"),
  silent = TRUE
)

pB <- as.ggplot(pB_grob)

## ---- 3. Combine into one multi-panel figure (patchwork) --------------------
## Adjust the width ratio below if a/b don't line up after rendering:
##   - if panel b looks too wide, decrease the second number (e.g. 1.3 -> 1.15)
##   - if panel b looks too narrow, increase it (e.g. 1.3 -> 1.45)
final_fig <- (pA + pB) +
  plot_layout(widths = c(0.8, 1.7)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = TAG_SIZE, face = "bold", family = FIG_FONT))

## ---- 4. Save -----------------------------------------------------------------
ggsave("Fig3_PCA_Correlation_combined.tiff", plot = final_fig,
       width = 27, height = 12, units = "cm", dpi = 600, compression = "lzw")
ggsave("Fig3_PCA_Correlation_combined.pdf", plot = final_fig,
       width = 27, height = 12, units = "cm", dpi = 600)

message("Done. Saved Fig3_PCA_Correlation_combined.tiff and .pdf")
