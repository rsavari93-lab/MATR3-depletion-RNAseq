## Figure 5 — Differential expression analysis of MATR3-depleted neural models
## Panels:
##   A) Volcano plot — siMATR3 vs siLUC, SH-SY5Y
##   B) Volcano plot — siMATR3 vs siLUC, U87
##   C) Venn diagram — shared DEGs between SH-SY5Y and U87 (84 genes)
##   D) Heatmap — 84 shared DEGs, all 12 samples, Z-scored, row-clustered
## Output: single combined figure (Fig4_DEG_combined.tiff / .pdf)


library(ggplot2)
library(readxl)
library(dplyr)
library(tibble)
library(tidyr)
library(VennDiagram)
library(grid)
library(pheatmap)
library(patchwork)
library(ggplotify)   # converts pheatmap/grid objects into ggplot-compatible objects

## ---- 1. Paths (edit these for your own machine) ----------------------------
setwd("D:/")

file_SH   <- "M3_SHvsLUC_SH_deg_all.xls"        # DEG table, SH-SY5Y
file_U87  <- "M3_U87vsLUC_U87_deg_all.xls"      # DEG table, U87
file_expr <- "MATR3_KD_mRNA_FPKM_matrix.xlsx"   # normalized expression matrix, all 12 samples

## ---- 1b. Robust file reader --------------------------------------------------
## Many bioinformatics pipelines (DESeq2, rMATS, etc.) export plain tab- or
## comma-delimited text but save it with an .xls/.xlsx extension. Excel parsers
## (readxl, openxlsx) then fail because the file isn't actually a binary/zip
## Excel file. This reader checks the real file content first and dispatches
## to the correct parser instead of trusting the extension.
read_table_robust <- function(path) {
  ## Peek at the first bytes to detect the real format
  raw_head <- readBin(path, "raw", n = 8)
  
  is_zip_xlsx <- length(raw_head) >= 4 &&
    raw_head[1] == 0x50 && raw_head[2] == 0x4B   # "PK" -> zip/xlsx
  is_ole_xls  <- length(raw_head) >= 4 &&
    raw_head[1] == 0xD0 && raw_head[2] == 0xCF   # legacy binary .xls
  
  if (is_zip_xlsx || is_ole_xls) {
    ## Genuinely an Excel file -> use readxl, with format-forcing fallback
    out <- tryCatch(read_excel(path), error = function(e) NULL)
    if (!is.null(out)) return(out)
    for (ext in c("xlsx", "xls", "xlsm")) {
      tmp <- tempfile(fileext = paste0(".", ext))
      file.copy(path, tmp, overwrite = TRUE)
      out <- tryCatch(read_excel(tmp), error = function(e) NULL)
      if (!is.null(out)) return(out)
    }
    if (!"openxlsx" %in% installed.packages()[, "Package"]) install.packages("openxlsx")
    out <- tryCatch(as_tibble(openxlsx::read.xlsx(path)), error = function(e) NULL)
    if (!is.null(out)) return(out)
    stop("File looks like Excel but could not be parsed: ", path)
  }
  
  ## Not a real Excel file -> it's plain text. Detect delimiter and read it.
  first_line <- readLines(path, n = 1, warn = FALSE)
  delim <- if (grepl("\t", first_line)) "\t" else if (grepl(";", first_line)) ";" else ","
  
  message("'", path, "' is plain text (not real Excel) -- reading as delimited text (delim = '",
          ifelse(delim == "\t", "\\t", delim), "').")
  
  out <- tryCatch(
    read.delim(path, sep = delim, header = TRUE, check.names = FALSE,
               stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(out)) stop("Could not parse plain-text file: ", path)
  as_tibble(out)
}

## ---- 2. Load DEG tables -----------------------------------------------------
deg_SH  <- read_table_robust(file_SH)
deg_U87 <- read_table_robust(file_U87)

## Standardize the gene-ID column name -> "gene"
## (assumes the first column is the gene identifier; edit if your column is named differently)
names(deg_SH)[1]  <- "gene"
names(deg_U87)[1] <- "gene"

deg_SH  <- deg_SH  %>% mutate(gene = trimws(as.character(gene)))
deg_U87 <- deg_U87 %>% mutate(gene = trimws(as.character(gene)))

## ---- 3. Classify regulation + significant gene sets ------------------------
classify_reg <- function(df) {
  df %>%
    mutate(
      regulation = case_when(
        log2FoldChange >  1 & padj < 0.05 ~ "UP",
        log2FoldChange < -1 & padj < 0.05 ~ "DOWN",
        TRUE ~ "NO"
      ),
      regulation = factor(regulation, levels = c("UP", "DOWN", "NO"))
    )
}

deg_SH  <- classify_reg(deg_SH)
deg_U87 <- classify_reg(deg_U87)

sig_SH  <- deg_SH  %>% filter(regulation %in% c("UP", "DOWN")) %>% pull(gene)
sig_U87 <- deg_U87 %>% filter(regulation %in% c("UP", "DOWN")) %>% pull(gene)

shared_genes <- intersect(sig_SH, sig_U87)
message("Number of shared DEGs: ", length(shared_genes))   # should print 84

## ---- 4. Shared color scheme (consistent across panels) ---------------------
col_up   <- "#E84040"
col_down <- "#3AB54A"
col_no   <- "#2171b5"

panel_theme <- theme_classic(base_size = 9) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    plot.title    = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.title    = element_text(size = 9),
    axis.text     = element_text(size = 8),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 7.5),
    legend.key.size = unit(0.35, "cm"),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.margin = margin(3, 3, 3, 3)
  )

## ---- 5. Panel A — Volcano, SH-SY5Y ------------------------------------------
make_volcano <- function(df, title_text) {
  n_up   <- sum(df$regulation == "UP")
  n_down <- sum(df$regulation == "DOWN")
  n_no   <- sum(df$regulation == "NO")
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = regulation, size = regulation), alpha = 0.9) +
    scale_size_manual(values = c("UP" = 1.2, "DOWN" = 1.2, "NO" = 0.5), guide = "none") +
    scale_color_manual(
      values = c("UP" = col_up, "DOWN" = col_down, "NO" = col_no),
      labels = c(paste("UP", n_up), paste("DOWN", n_down), paste("NO", n_no))
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 1.301, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    labs(title = title_text, x = "log2FoldChange", y = "-log10(padj)",
         color = "padj < 0.05\n|log2FC| > 1") +
    panel_theme
}

pA <- make_volcano(deg_SH,  "siMATR3 vs siLUC — SH-SY5Y")
pB <- make_volcano(deg_U87, "siMATR3 vs siLUC — U87")

## ---- 6. Panel C — Venn diagram (U87 vs SH-SY5Y) -----------------------------
venn_list <- list(U87 = sig_U87, "SH-SY5Y" = sig_SH)

venn_grob <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("#5B9BD5", "#F2C811"),
  alpha = 0.65,
  col = "white",
  cex = 1.1,                # number font size
  cat.cex = 1.1,            # set-label font size
  cat.fontface = "bold",
  cat.col = c("#5B9BD5", "#B8860B"),
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.08,
  lwd = 1
)

pC <- as.ggplot(grid::grobTree(venn_grob)) +
  theme(plot.margin = margin(5, 5, 5, 5))

## ---- 7. Panel D — Heatmap of 84 shared DEGs ---------------------------------
expr <- read_table_robust(file_expr)
names(expr)[1] <- "gene"
expr <- expr %>% mutate(gene = trimws(as.character(gene)))

expr_mat_raw <- expr %>%
  filter(gene %in% shared_genes) %>%
  distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames("gene")

## Keep only actual sample columns (siLUC_*/siMATR3_* ... _mRNA), drop
## annotation columns such as gene_name, gene_chr, gene_biotype, tf_family, etc.
sample_cols <- grep("^si(LUC|MATR3)_", colnames(expr_mat_raw), value = TRUE)
if (length(sample_cols) == 0) {
  stop("No sample columns matched pattern '^si(LUC|MATR3)_' in '", file_expr,
       "'. Found columns: ", paste(colnames(expr_mat_raw), collapse = ", "))
}
message("Using ", length(sample_cols), " sample columns out of ",
        ncol(expr_mat_raw), " total columns (dropped annotation columns: ",
        paste(setdiff(colnames(expr_mat_raw), sample_cols), collapse = ", "), ")")

expr_mat <- as.matrix(expr_mat_raw[, sample_cols])

## ---- 7b. Guard against non-numeric columns/values --------------------------
## Common causes: a stray annotation column left in, "NA"/"" strings, or
## decimal commas instead of dots. Coerce safely and warn if anything is lost.
if (!is.numeric(expr_mat)) {
  message("Expression matrix is not purely numeric -- attempting safe coercion.")
  
  storage_before <- expr_mat
  expr_mat_num <- suppressWarnings(apply(storage_before, 2, function(col) {
    as.numeric(gsub(",", ".", trimws(col)))
  }))
  rownames(expr_mat_num) <- rownames(storage_before)
  colnames(expr_mat_num) <- colnames(storage_before)
  
  n_new_na <- sum(is.na(expr_mat_num) & !is.na(storage_before) & storage_before != "")
  if (n_new_na > 0) {
    bad_cols <- colnames(storage_before)[colSums(
      is.na(expr_mat_num) & !is.na(storage_before) & storage_before != ""
    ) > 0]
    warning(n_new_na, " value(s) could not be converted to numeric and became NA.",
            " Affected column(s): ", paste(bad_cols, collapse = ", "),
            ". Check these columns in '", file_expr, "' -- they may contain text,",
            " merged cells, or a non-sample column (e.g. gene symbol/annotation)",
            " that should be removed or set as rownames instead.")
  }
  expr_mat <- expr_mat_num
}

if (any(is.na(expr_mat))) {
  message(sum(is.na(expr_mat)), " NA value(s) remain in the expression matrix; ",
          "these genes/samples may show blank cells in the heatmap.")
}

## log2-transform (add pseudocount) then row Z-score, matching Fig 4D legend
expr_log <- log2(expr_mat + 1)
expr_z   <- t(scale(t(expr_log)))   # row-wise Z-score

## Sample annotation — EDIT sample name patterns if yours differ
sample_names <- colnames(expr_z)
cell_line <- ifelse(grepl("SH", sample_names, ignore.case = TRUE), "SH-SY5Y", "U87")
condition <- ifelse(grepl("siMATR3|MATR3", sample_names, ignore.case = TRUE), "MATR3 KD", "Control")

annot_col <- data.frame(CellLine = cell_line, Condition = condition,
                        row.names = sample_names)

annot_colors <- list(
  CellLine  = c("SH-SY5Y" = "#1F4E79", "U87" = "#F2C811"),
  Condition = c("Control" = "#5B9BD5", "MATR3 KD" = "#C0392B")
)

pD_grob <- pheatmap(
  expr_z,
  annotation_col   = annot_col,
  annotation_colors = annot_colors,
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  fontsize_row     = 4.5,
  fontsize_col     = 7,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  border_color = NA,
  treeheight_row = 25,
  treeheight_col = 15,
  main = "Shared differentially expressed genes",
  silent = TRUE
)

pD <- as.ggplot(pD_grob)

## ---- 8. Combine into one multi-panel figure (patchwork) --------------------
top_row    <- pA + pB + plot_layout(widths = c(1, 1))
bottom_row <- pC + pD + plot_layout(widths = c(0.85, 1.6))

final_fig <- (top_row / bottom_row) +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

## ---- 9. Save -----------------------------------------------------------------
ggsave("Fig4_DEG_combined.tiff", plot = final_fig, width = 22, height = 18, units = "cm", dpi = 300, compression = "lzw")
ggsave("Fig4_DEG_combined.pdf",  plot = final_fig, width = 22, height = 18, units = "cm", dpi = 300)

message("Done. Saved Fig4_DEG_combined.tiff and .pdf")

getwd()
