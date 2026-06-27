## =============================================================================
## Figure 6 — Alternative splicing analysis following MATR3 knockdown
## Panels:
##   a) Bar chart of significant splicing events per class (SE, A3SS, A5SS, MXE, RI)
##   b) Venn diagram of genes with shared splicing alterations (SH-SY5Y vs U87)
## Output: single combined figure (Fig5_Splicing_combined.tiff / .pdf)

library(ggplot2)
library(ggforce)
library(patchwork)

setwd("D:/")

## ---- 0b. Journal-standard font settings (same constants used in Fig 3 & 4) --
FIG_FONT          <- "sans"
TITLE_SIZE        <- 10
TITLE_FACE        <- "bold"
AXIS_TITLE_SIZE   <- 9
AXIS_TEXT_SIZE    <- 8
LEGEND_TITLE_SIZE <- 8
LEGEND_TEXT_SIZE  <- 7.5
TAG_SIZE          <- 11

## ---- 1. Panel a — Splicing event bar chart ----------------------------------
data <- data.frame(
  Event = rep(c("SE", "A3SS", "A5SS", "MXE", "RI"), 2),
  Count = c(208, 13, 13, 17, 5, 679, 79, 65, 79, 59),
  CellLine = c(rep("SH-SY5Y", 5), rep("U87", 5))
)
data$Event <- factor(data$Event, levels = c("SE", "A3SS", "A5SS", "MXE", "RI"))
data$CellLine <- factor(data$CellLine, levels = c("SH-SY5Y", "U87"))

pA <- ggplot(data, aes(x = Event, y = Count, fill = CellLine)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.65),
           width = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = c("SH-SY5Y" = "#1a6faf", "U87" = "#f5c800"),
    labels = c("SH-SY5Y (neuroblastoma)", "U87 (glioblastoma)")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 750), breaks = seq(0, 700, 100)) +
  labs(x = "Alternative splicing event type",
       y = "Significant events (FDR < 0.05)", fill = NULL) +
  theme_classic(base_size = AXIS_TEXT_SIZE, base_family = FIG_FONT) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line     = element_line(color = "black", linewidth = 0.5),
    axis.title.x  = element_text(size = AXIS_TITLE_SIZE, color = "black",
                                  family = FIG_FONT, margin = margin(t = 8)),
    axis.title.y  = element_text(size = AXIS_TITLE_SIZE, color = "black",
                                  family = FIG_FONT, margin = margin(r = 8)),
    axis.text     = element_text(size = AXIS_TEXT_SIZE, color = "black", family = FIG_FONT),
    legend.position = "top",
    legend.key.size = unit(0.4, "cm"),
    legend.title  = element_text(size = LEGEND_TITLE_SIZE, family = FIG_FONT),
    legend.text   = element_text(size = LEGEND_TEXT_SIZE, family = FIG_FONT),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.margin = margin(15, 20, 10, 15)
  )

## ---- 2. Panel b — Venn diagram of shared splicing-associated genes ----------
venn_data <- data.frame(
  x = c(-0.3, 0.3),
  y = c(0, 0),
  r = c(0.5, 0.5),
  label = c("U87", "SH-SY5Y")
)

pB <- ggplot() +
  geom_circle(data = venn_data,
              aes(x0 = x, y0 = y, r = r, fill = label),
              alpha = 0.6, color = "grey40", linewidth = 0.8) +
  scale_fill_manual(values = c("U87" = "#1a6faf", "SH-SY5Y" = "#f5c800")) +
  annotate("text", x = -0.42, y = 0, label = "784", size = AXIS_TEXT_SIZE / .pt,
           fontface = TITLE_FACE, color = "white", family = FIG_FONT) +
  annotate("text", x = 0, y = 0, label = "79", size = AXIS_TEXT_SIZE / .pt,
           fontface = TITLE_FACE, color = "white", family = FIG_FONT) +
  annotate("text", x = 0.42, y = 0, label = "162", size = AXIS_TEXT_SIZE / .pt,
           fontface = TITLE_FACE, color = "black", family = FIG_FONT) +
  annotate("text", x = -0.45, y = 0.62, label = "U87", size = AXIS_TITLE_SIZE / .pt,
           fontface = TITLE_FACE, color = "#1a6faf", family = FIG_FONT) +
  annotate("text", x = 0.45, y = 0.62, label = "SH-SY5Y", size = AXIS_TITLE_SIZE / .pt,
           fontface = TITLE_FACE, color = "#c9a800", family = FIG_FONT) +
  labs(title = "Shared splicing-associated genes") +
  coord_fixed() +
  theme_void(base_family = FIG_FONT) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = TITLE_SIZE, face = TITLE_FACE,
                               hjust = 0.5, family = FIG_FONT,
                               margin = margin(b = 8))
  )

## ---- 3. Combine into one multi-panel figure (patchwork) --------------------
## Adjust the width ratio below if a/b don't look balanced after rendering:
##   - if panel b (Venn) looks too small, increase its number (e.g. 0.8 -> 1)
##   - if panel b looks too large, decrease it (e.g. 0.8 -> 0.65)
final_fig <- (pA + pB) +
  plot_layout(widths = c(1, 0.8)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = TAG_SIZE, face = "bold", family = FIG_FONT))

## ---- 4. Save -----------------------------------------------------------------
ggsave("Fig5_Splicing_combined.tiff", plot = final_fig,
       width = 20, height = 10, units = "cm", dpi = 300, compression = "lzw")
ggsave("Fig5_Splicing_combined.pdf", plot = final_fig,
       width = 20, height = 10, units = "cm", dpi = 300)

message("Done. Saved Fig5_Splicing_combined.tiff and .pdf")
