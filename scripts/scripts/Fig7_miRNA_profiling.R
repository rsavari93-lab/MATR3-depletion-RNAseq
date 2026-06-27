## Figure 7 — MicroRNA expression profiling following MATR3 knockdown
## Panels:
##   a) Volcano plot — miRNA, siMATR3 vs siLUC, SH-SY5Y (let-7c-3p highlighted)
##   b) Volcano plot — miRNA, siMATR3 vs siLUC, U87 (let-7c-3p highlighted)
##   c) Bar chart — let-7c-3p normalized read counts, both cell lines
## Output: single combined figure (Fig7_miRNA_combined.tiff / .pdf)
## ============================================================================
library(ggplot2)
library(readxl)
library(patchwork)

setwd("D:/")

## ---- 0b. Journal-standard font settings (same constants used in Fig 3, 4, 5) ----
FIG_FONT          <- "sans"
TITLE_SIZE        <- 10
TITLE_FACE        <- "bold"
AXIS_TITLE_SIZE   <- 9
AXIS_TEXT_SIZE    <- 8
LEGEND_TITLE_SIZE <- 8
LEGEND_TEXT_SIZE  <- 7.5
TAG_SIZE          <- 11

## ---- 1. Load + classify miRNA DEG tables ------------------------------------
mirna_sh  <- read_excel("M3_SH_SvsLUC_SH_S.Differential_analysis_results.xlsx")
mirna_u87 <- read_excel("M3_U87_SvsLUC_U87_S.Differential_analysis_results.xlsx")

mirna_sh$log2FoldChange  <- as.numeric(mirna_sh$log2FoldChange)
mirna_sh$pval            <- as.numeric(mirna_sh$pval)
mirna_u87$log2FoldChange <- as.numeric(mirna_u87$log2FoldChange)
mirna_u87$pval           <- as.numeric(mirna_u87$pval)

mirna_sh  <- mirna_sh[!is.na(mirna_sh$pval) & !is.na(mirna_sh$log2FoldChange), ]
mirna_u87 <- mirna_u87[!is.na(mirna_u87$pval) & !is.na(mirna_u87$log2FoldChange), ]

mirna_sh$regulation <- "NO"
mirna_sh$regulation[mirna_sh$log2FoldChange >  1 & mirna_sh$pval < 0.05] <- "UP"
mirna_sh$regulation[mirna_sh$log2FoldChange < -1 & mirna_sh$pval < 0.05] <- "DOWN"
mirna_sh$regulation <- factor(mirna_sh$regulation, levels = c("UP", "DOWN", "NO"))

mirna_u87$regulation <- "NO"
mirna_u87$regulation[mirna_u87$log2FoldChange >  1 & mirna_u87$pval < 0.05] <- "UP"
mirna_u87$regulation[mirna_u87$log2FoldChange < -1 & mirna_u87$pval < 0.05] <- "DOWN"
mirna_u87$regulation <- factor(mirna_u87$regulation, levels = c("UP", "DOWN", "NO"))

up_sh   <- sum(mirna_sh$regulation == "UP")
down_sh <- sum(mirna_sh$regulation == "DOWN")
no_sh   <- sum(mirna_sh$regulation == "NO")
up_u87   <- sum(mirna_u87$regulation == "UP")
down_u87 <- sum(mirna_u87$regulation == "DOWN")
no_u87   <- sum(mirna_u87$regulation == "NO")

## ---- 2. Shared volcano theme (matches Fig 3/4/5 standard) -------------------
volcano_theme <- theme_classic(base_size = AXIS_TEXT_SIZE, base_family = FIG_FONT) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    plot.title    = element_text(size = TITLE_SIZE, hjust = 0.5, face = TITLE_FACE, family = FIG_FONT),
    axis.title    = element_text(size = AXIS_TITLE_SIZE, family = FIG_FONT),
    axis.text     = element_text(size = AXIS_TEXT_SIZE, family = FIG_FONT),
    legend.title  = element_text(size = LEGEND_TITLE_SIZE, family = FIG_FONT),
    legend.text   = element_text(size = LEGEND_TEXT_SIZE, family = FIG_FONT),
    legend.background = element_rect(fill = "white", color = "grey80")
  )

## ---- 3. Panel a — Volcano, SH-SY5Y -------------------------------------------
pA <- ggplot(mirna_sh, aes(x = log2FoldChange, y = -log10(pval))) +
  geom_point(aes(color = regulation, size = regulation), alpha = 0.9) +
  scale_size_manual(values = c("UP" = 2.5, "DOWN" = 2.5, "NO" = 0.8), guide = "none") +
  scale_color_manual(
    values = c("UP" = "#E84040", "DOWN" = "#3AB54A", "NO" = "#2171b5"),
    labels = c(paste("UP", up_sh), paste("DOWN", down_sh), paste("NO", no_sh))
  ) +
  geom_point(data = mirna_sh[mirna_sh$sRNA == "hsa-let-7c-3p", ],
             aes(x = log2FoldChange, y = -log10(pval)),
             color = "#E84040", size = 3, shape = 22, stroke = 1.1, fill = NA) +
  geom_text(data = mirna_sh[mirna_sh$sRNA == "hsa-let-7c-3p", ],
            aes(x = log2FoldChange, y = -log10(pval), label = "let-7c-3p"),
            hjust = -0.3, vjust = 0.5, size = AXIS_TEXT_SIZE / .pt,
            fontface = "bold.italic", family = FIG_FONT) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(title = "siMATR3 vs siLUC — SH-SY5Y",
       x = "log2FoldChange", y = "-log10(p-value)",
       color = "p-value < 0.05\n|log2FC| > 1") +
  volcano_theme

## ---- 4. Panel b — Volcano, U87 -----------------------------------------------
pB <- ggplot(mirna_u87, aes(x = log2FoldChange, y = -log10(pval))) +
  geom_point(aes(color = regulation, size = regulation), alpha = 0.9) +
  scale_size_manual(values = c("UP" = 2.5, "DOWN" = 2.5, "NO" = 0.8), guide = "none") +
  scale_color_manual(
    values = c("UP" = "#E84040", "DOWN" = "#3AB54A", "NO" = "#2171b5"),
    labels = c(paste("UP", up_u87), paste("DOWN", down_u87), paste("NO", no_u87))
  ) +
  geom_point(data = mirna_u87[mirna_u87$sRNA == "hsa-let-7c-3p", ],
             aes(x = log2FoldChange, y = -log10(pval)),
             color = "#E84040", size = 3, shape = 22, stroke = 1.1, fill = NA) +
  geom_text(data = mirna_u87[mirna_u87$sRNA == "hsa-let-7c-3p", ],
            aes(x = log2FoldChange, y = -log10(pval), label = "let-7c-3p"),
            hjust = -0.3, vjust = 0.5, size = AXIS_TEXT_SIZE / .pt,
            fontface = "bold.italic", family = FIG_FONT) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(title = "siMATR3 vs siLUC — U87",
       x = "log2FoldChange", y = "-log10(p-value)",
       color = "p-value < 0.05\n|log2FC| > 1") +
  volcano_theme

## ---- 5. Panel c — let-7c-3p expression bar chart -----------------------------
bardata <- data.frame(
  Group = c("SH-SY5Y", "SH-SY5Y", "U87", "U87"),
  Condition = c("Control", "MATR3 KD", "Control", "MATR3 KD"),
  TPM = c(79.5, 297.6, 45.2, 160.4)
)
bardata$Group <- factor(bardata$Group, levels = c("SH-SY5Y", "U87"))
bardata$Condition <- factor(bardata$Condition, levels = c("Control", "MATR3 KD"))

pC <- ggplot(bardata, aes(x = Group, y = TPM, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Control" = "#1a6faf", "MATR3 KD" = "#f5c800")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 350)) +
  annotate("text", x = 1.175, y = 315, label = "*", size = 6) +
  annotate("segment", x = 0.825, xend = 1.525, y = 308, yend = 308) +
  annotate("text", x = 2.175, y = 185, label = "*", size = 6) +
  annotate("segment", x = 1.825, xend = 2.525, y = 178, yend = 178) +
  labs(title = "let-7c-3p expression",
       x = "Cell Line", y = "Normalized read count", fill = NULL) +
  theme_classic(base_size = AXIS_TEXT_SIZE, base_family = FIG_FONT) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    plot.title    = element_text(size = TITLE_SIZE, hjust = 0.5, face = "bold.italic", family = FIG_FONT),
    axis.title.x  = element_text(size = AXIS_TITLE_SIZE, family = FIG_FONT, margin = margin(t = 8)),
    axis.title.y  = element_text(size = AXIS_TITLE_SIZE, family = FIG_FONT),
    axis.text.x   = element_text(size = AXIS_TEXT_SIZE, family = FIG_FONT, margin = margin(t = 5)),
    axis.text.y   = element_text(size = AXIS_TEXT_SIZE, family = FIG_FONT),
    legend.title  = element_text(size = LEGEND_TITLE_SIZE, family = FIG_FONT),
    legend.text   = element_text(size = LEGEND_TEXT_SIZE, family = FIG_FONT)
  )

## ---- 6. Combine into one multi-panel figure (patchwork) --------------------
final_fig <- (pA | pB | pC) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = TAG_SIZE, face = "bold", family = FIG_FONT),
        plot.background = element_rect(fill = "white", color = NA))

## ---- 7. Save -----------------------------------------------------------------
ggsave("Fig7_miRNA_combined.tiff", plot = final_fig,
       width = 24, height = 9, units = "cm", dpi = 300, compression = "lzw")
ggsave("Fig7_miRNA_combined.pdf", plot = final_fig,
       width = 24, height = 9, units = "cm", dpi = 300)

message("Done. Saved Fig7_miRNA_combined.tiff and .pdf")
