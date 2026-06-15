library(ggplot2)
library(readxl)

setwd("D:/")
data <- read_excel("M3_SHvsLUC_SH_deg_all.xlsx")

data$regulation <- "NO"
data$regulation[data$log2FoldChange > 1 & data$padj < 0.05] <- "UP"
data$regulation[data$log2FoldChange < -1 & data$padj < 0.05] <- "DOWN"
data$regulation <- factor(data$regulation, levels=c("UP","DOWN","NO"))

up_count <- sum(data$regulation == "UP")
down_count <- sum(data$regulation == "DOWN")
no_count <- sum(data$regulation == "NO")

p <- ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=regulation, size=regulation), alpha=0.9) +
  scale_size_manual(values=c("UP"=2,"DOWN"=2,"NO"=0.8), guide="none") +
  scale_color_manual(
    values=c("UP"="#E84040","DOWN"="#3AB54A","NO"="#2171b5"),
    labels=c(paste("UP",up_count), paste("DOWN",down_count), paste("NO",no_count))
  ) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept=1.301, linetype="dashed", color="grey50") +
  labs(title="siMATR3 vs siLUC — SH-SY5Y",
       x="log2FoldChange", y="-log10(padj)",
       color="padj < 0.05\n|log2FoldChange| > 1") +
  theme_classic() +
  theme(
    plot.background=element_rect(fill="white", color=NA),
    panel.background=element_rect(fill="white"),
    plot.title=element_text(size=12, hjust=0.5, face="bold"),
    axis.title=element_text(size=11),
    axis.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.text=element_text(size=9),
    legend.background=element_rect(fill="white", color="grey80"),
    legend.margin=margin(5,5,5,5)
  )

ggsave("Volcano_SH.tiff", plot=p, width=10, height=9, units="cm", dpi=300)
