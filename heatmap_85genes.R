library(pheatmap)
library(readxl)
library(RColorBrewer)

setwd("D:/")
data <- read_excel("GENE.COUNT.xlsm")

shared_genes <- c("RMI2","RPS27L","SCAMP2","YRDC","SNRPB2","RHOG","KCNJ9",
                  "FAM89A","TTPAL","CMTM3","SMCO4","ATP6V0E1","STMN4","UBE2E3",
                  "ARL1","PSTPIP2","MANF","RPIA","CX3CL1","QDPR","GINS2","UBE2W",
                  "RHPN2","FTL","ZFAND3","TSPAN7","ACBD5","UMAD1","PRPS2","ZNF268",
                  "SFXN1","CCNG1","TSEN34","INTS5","PMP22","ZNF271P","EEF1AKMT3",
                  "CYB5B","MYCBP","CDK4","PRLR","RAC1","LAMP1","LSM4","TNS3",
                  "LMNB1","TCF7","GABRA3","SLITRK6","AC016394.1","SYPL1","JPT2",
                  "AL159169.2","PIK3C2B","GNG3","ASF1B","CMTM6","CMAHP","CFL2",
                  "PGM2L1","PINK1-AS","SDCBP","ZNF35","RABL3","AP003486.1",
                  "AC026740.1","PLIN3","MFSD4A","POLE3","VSTM2L","DTL","PAK6",
                  "AC093218.1","AC026368.1","RAB3A","SPTB","GLRX5","IGSF3",
                  "SUCLG2","PARD6B","LGI1","PTGES","POM121","EPHB3")

data_filtered <- data[data$gene_name %in% shared_genes,]

mat <- as.matrix(data_filtered[, c("siLUC_1_SH","siLUC_2_SH","siLUC_3_SH",
                                    "siMATR3_1_SH","siMATR3_2_SH","siMATR3_3_SH",
                                    "siLUC_1_U87","siLUC_2_U87","siLUC_3_U87",
                                    "siMATR3_1_U87","siMATR3_2_U87","siMATR3_3_U87")])
rownames(mat) <- data_filtered$gene_name

colnames(mat) <- c("siLUC SH 1","siLUC SH 2","siLUC SH 3",
                   "siMATR3 SH 1","siMATR3 SH 2","siMATR3 SH 3",
                   "siLUC U87 1","siLUC U87 2","siLUC U87 3",
                   "siMATR3 U87 1","siMATR3 U87 2","siMATR3 U87 3")

mat_log <- log2(mat + 1)
mat_scaled <- t(scale(t(mat_log)))

col_annotation <- data.frame(
  Condition = c(rep("Control",3), rep("MATR3 KD",3),
                rep("Control",3), rep("MATR3 KD",3)),
  CellLine  = c(rep("SH-SY5Y",6), rep("U87",6))
)
rownames(col_annotation) <- colnames(mat_scaled)

ann_colors <- list(
  Condition = c("Control"="#888888", "MATR3 KD"="#E84040"),
  CellLine  = c("SH-SY5Y"="#1a6faf", "U87"="#f5c800")
)

tiff("Heatmap_85genes.tiff", width=8, height=12, units="in", res=300)
pheatmap(mat_scaled,
         annotation_col    = col_annotation,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         show_colnames     = TRUE,
         fontsize_row      = 7,
         fontsize_col      = 8,
         angle_col         = 45,
         border_color      = NA,
         main              = "Shared differentially expressed genes")
dev.off()
