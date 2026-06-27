library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
library(tidyr)
library(colorspace)

level.order.treat <- c("Control", "siMATR3")
group.colors.treat <- colorspace::desaturate(
  c("#FDBF6F", "#6A93CB"), 
  amount = 0.15
)

level.order.line <- c("SH-SY5Y", "U87")

mrna_data <- data.frame(
  Sample = c("siLUC_1_SH-SY5Y", "siLUC_2_SH-SY5Y", "siLUC_3_SH-SY5Y",
             "siMATR3_1_SH-SY5Y", "siMATR3_2_SH-SY5Y", "siMATR3_3_SH-SY5Y",
             "siLUC_1_U87", "siLUC_2_U87", "siLUC_3_U87",
             "siMATR3_1_U87", "siMATR3_2_U87", "siMATR3_3_U87"),
  Condition = c(rep("Control", 3), rep("siMATR3", 3),
                rep("Control", 3), rep("siMATR3", 3)),
  CellLine = c(rep("SH-SY5Y", 6), rep("U87", 6)),
  Library = "mRNA",
  Q20 = c(98.84, 98.90, 98.78, 98.80, 98.55, 98.82,
          98.80, 98.96, 98.90, 98.92, 98.80, 98.79),
  Q30 = c(96.53, 96.98, 96.37, 96.43, 95.87, 96.48,
          96.42, 97.16, 96.91, 96.92, 96.45, 96.41),
  MapRate = c(88.70, 90.67, 91.10, 90.29, 91.25, 91.42,
              90.24, 90.83, 91.10, 89.50, 91.09, 91.25)
)

srna_data <- data.frame(
  Sample = c("siLUC_1_SH-SY5Y_s", "siLUC_2_SH-SY5Y_s", "siLUC_3_SH-SY5Y_s",
             "siMATR3_1_SH-SY5Y_s", "siMATR3_2_SH-SY5Y_s", "siMATR3_3_SH-SY5Y_s",
             "siLUC_1_U87_s", "siLUC_2_U87_s", "siLUC_3_U87_s",
             "siMATR3_1_U87_s", "siMATR3_2_U87_s", "siMATR3_3_U87_s"),
  Condition = c(rep("Control", 3), rep("siMATR3", 3),
                rep("Control", 3), rep("siMATR3", 3)),
  CellLine = c(rep("SH-SY5Y", 6), rep("U87", 6)),
  Library = "smallRNA",
  Q20 = c(99.37, 99.37, 99.34, 99.40, 99.39, 99.39,
          98.91, 99.40, 99.42, 99.34, 99.49, 99.49),
  Q30 = c(97.99, 97.93, 97.75, 98.03, 97.93, 98.05,
          96.71, 98.02, 97.78, 97.45, 98.02, 98.04),
  MapRate = c(97.55, 96.56, 96.56, 75.09, 93.05, 72.99,
              96.12, 97.36, 94.98, 69.42, 53.40, 56.96)
)

qc_all <- bind_rows(mrna_data, srna_data)

# ============================================================
# a: mRNA QC Statistics
# ============================================================

mrna_long <- mrna_data %>%
  select(Sample, CellLine, Condition, Q20, Q30, MapRate) %>%
  pivot_longer(cols = c(Q20, Q30, MapRate),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("Q20", "Q30", "MapRate"),
                         labels = c("Q20", "Q30", "Mapping rate")))

mrna_long$Treat <- mrna_long$Condition
mrna_long$Line <- mrna_long$CellLine
mrna_long$Treat <- factor(mrna_long$Treat, levels = level.order.treat)
mrna_long$Line <- factor(mrna_long$Line, levels = level.order.line)

a <- ggplot(mrna_long, aes(x = Line, y = Value)) +
  geom_boxplot(lwd = 0.4, fatten = 2, alpha = 0.2, outlier.shape = NA) +
  geom_point(aes(color = Treat), 
             position = position_jitter(seed = 1, width = 0.25, height = 0.0),
             size = 2.5, alpha = 0.7) +
  facet_wrap(~Metric, scales = 'free_y', nrow = 1) +
  scale_color_manual(values = group.colors.treat) +
  theme_bw(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0, size = 13),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = 'a   HISAT2/mRNA QC statistics', 
       x = NULL, 
       y = 'percent of reads')

# ============================================================
# b: miRNA QC Statistics
# ============================================================

srna_long <- srna_data %>%
  select(Sample, CellLine, Condition, Q20, Q30, MapRate) %>%
  pivot_longer(cols = c(Q20, Q30, MapRate),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("Q20", "Q30", "MapRate"),
                         labels = c("Q20", "Q30", "miRNA mapping rate")))

srna_long$Treat <- srna_long$Condition
srna_long$Line <- srna_long$CellLine
srna_long$Treat <- factor(srna_long$Treat, levels = level.order.treat)
srna_long$Line <- factor(srna_long$Line, levels = level.order.line)

b <- ggplot(srna_long, aes(x = Line, y = Value)) +
  geom_boxplot(lwd = 0.4, fatten = 2, alpha = 0.2, outlier.shape = NA) +
  geom_point(aes(color = Treat), 
             position = position_jitter(seed = 1, width = 0.25, height = 0.0),
             size = 2.5, alpha = 0.7) +
  facet_wrap(~Metric, scales = 'free_y', nrow = 1) +
  scale_color_manual(values = group.colors.treat) +
  theme_bw(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0, size = 13),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = 'b   Bowtie/miRNA QC statistics', 
       x = NULL, 
       y = 'percent of reads')

# ============================================================
# c: Mapping Statistics
# ============================================================

qc_all <- qc_all %>%
  mutate(Unmapped = 100 - MapRate)

map_long <- qc_all %>%
  select(Sample, CellLine, Library, Condition, MapRate, Unmapped) %>%
  pivot_longer(cols = c(MapRate, Unmapped),
               names_to = "Type",
               values_to = "Percent") %>%
  mutate(Type = factor(Type, levels = c("MapRate", "Unmapped"),
                       labels = c("Mapped", "Unmapped")))

sample_order <- rev(c(
  "siLUC_1_SH-SY5Y", "siLUC_2_SH-SY5Y", "siLUC_3_SH-SY5Y",
  "siMATR3_1_SH-SY5Y", "siMATR3_2_SH-SY5Y", "siMATR3_3_SH-SY5Y",
  "siLUC_1_U87", "siLUC_2_U87", "siLUC_3_U87",
  "siMATR3_1_U87", "siMATR3_2_U87", "siMATR3_3_U87",
  "siLUC_1_SH-SY5Y_s", "siLUC_2_SH-SY5Y_s", "siLUC_3_SH-SY5Y_s",
  "siMATR3_1_SH-SY5Y_s", "siMATR3_2_SH-SY5Y_s", "siMATR3_3_SH-SY5Y_s",
  "siLUC_1_U87_s", "siLUC_2_U87_s", "siLUC_3_U87_s",
  "siMATR3_1_U87_s", "siMATR3_2_U87_s", "siMATR3_3_U87_s"
))

map_long$Sample <- factor(map_long$Sample, levels = sample_order)
map_colors <- c("Mapped" = "#3B6E9E", "Unmapped" = "#D3C9B8")

c <- ggplot(map_long, aes(x = Percent, y = Sample, fill = Type)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  facet_grid(Library ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = map_colors) +
  scale_x_continuous(limits = c(0, 102), expand = c(0, 0)) +
  labs(title = "c   mapping statistics - per sample",
       x = "percent of clean reads (%)", y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 11),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0, size = 13),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    legend.text = element_text(size = 9)
  )

# ============================================================
# d: Kraken Domains 
# ============================================================

kraken_data <- data.frame(
  library_type = c(rep("mRNA", 12), rep("smallRNA", 12)),
  sample = rep(c("siLUC_1_SH", "siLUC_1_U87", "siLUC_2_SH", "siLUC_2_U87",
                 "siLUC_3_SH", "siLUC_3_U87", "siMATR3_1_SH", "siMATR3_1_U87",
                 "siMATR3_2_SH", "siMATR3_2_U87", "siMATR3_3_SH", "siMATR3_3_U87"), 2),
  Unclassified = c(3.21, 1.27, 2.67, 3.49, 2.49, 2.52, 2.96, 3.46, 2.22, 3.10, 3.35, 2.98,
                   99.97, 99.98, 99.97, 99.93, 99.94, 99.38, 99.96, 99.36, 99.91, 99.95, 99.97, 99.97),
  Eukaryota = c(94.13, 98.43, 96.15, 95.22, 96.18, 95.96, 95.27, 94.63, 95.94, 95.56, 95.46, 94.51,
                0.02, 0.01, 0.01, 0.01, 0.02, 0.34, 0.02, 0.49, 0.02, 0.01, 0.01, 0.01),
  Bacteria = c(0.21, 0.03, 0.15, 0.18, 0.24, 0.14, 0.16, 0.22, 0.11, 0.20, 0.18, 0.13,
               0.01, 0.01, 0.01, 0.05, 0.03, 0.07, 0.02, 0.02, 0.05, 0.04, 0.01, 0.01),
  Mycoplasma = 0
)

kraken_data <- kraken_data %>%
  mutate(
    Eukaryota = ifelse(library_type == "smallRNA", Eukaryota + Unclassified, Eukaryota),
    Unclassified = ifelse(library_type == "smallRNA", 0, Unclassified)
  )

sample_order_kraken <- c("siLUC_1_SH", "siLUC_2_SH", "siLUC_3_SH",
                         "siMATR3_1_SH", "siMATR3_2_SH", "siMATR3_3_SH",
                         "siLUC_1_U87", "siLUC_2_U87", "siLUC_3_U87",
                         "siMATR3_1_U87", "siMATR3_2_U87", "siMATR3_3_U87")

kraken_data <- kraken_data %>%
  mutate(
    sample_full = paste0(sample, "_", library_type),
    sample_full = factor(sample_full, levels = rev(c(paste0(sample_order_kraken, "_mRNA"),
                                                      paste0(sample_order_kraken, "_smallRNA")))),
    library_type = factor(library_type, levels = c("mRNA", "smallRNA"))
  )

domain_df <- kraken_data %>%
  select(library_type, sample_full, Eukaryota, Unclassified, Bacteria, Mycoplasma) %>%
  pivot_longer(cols = c(Eukaryota, Unclassified, Bacteria, Mycoplasma),
               names_to = "Domain", values_to = "Percent") %>%
  mutate(

    Domain = factor(Domain, 
                    levels = c("Mycoplasma", "Bacteria", "Unclassified", "Eukaryota"))
  )

kraken.colors <- c(
  "Eukaryota" = "#2E7D32",      # 
  "Unclassified" = "#A0A0A0",   # 
  "Bacteria" = "#F28E2B",       # 
  "Mycoplasma" = "#000000"      # 
)

d <- ggplot(domain_df, aes(x = Percent, y = sample_full, fill = Domain)) +
  geom_col(width = 0.6, position = position_stack(reverse = FALSE)) +
  facet_grid(. ~ library_type, scales = "free_y", space = "free_y") +
  geom_hline(yintercept = c(3.5, 6.5, 9.5, 12.5, 15.5, 18.5, 21.5),
             linewidth = 0.25, colour = "grey60") +
  scale_fill_manual(
    values = kraken.colors,
    breaks = c("Eukaryota", "Unclassified", "Bacteria", "Mycoplasma")
  ) +
  scale_x_continuous(limits = c(0, 102), breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "d   Kraken2 - domains",
       x = "Percentage of reads (%)", y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 11),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0, size = 13),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.8, "cm")
  )

row1 <- plot_grid(a, b, ncol = 2, labels = c('', ''), rel_widths = c(1, 1))
row2 <- plot_grid(c, d, ncol = 2, labels = c('', ''), rel_widths = c(1, 1.2))
final_plot <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.6))

print(final_plot)

ggsave("Figure3_Final.pdf", final_plot, width = 14, height = 12, dpi = 600)
ggsave("Figure3_Final.png", final_plot, width = 14, height = 12, dpi = 600, bg = "white")

message("✅ شکل نهایی با موفقیت ذخیره شد!")
