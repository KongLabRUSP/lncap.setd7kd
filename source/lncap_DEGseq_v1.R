# |----------------------------------------------------------------------------------|
# | Project: Study of LNCaP WT/KD cells                                              |
# | Script: Analysis of microarray data from a LNCaP WT/KD experiment                |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 02/28/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/lncap_DEGseq_v1.R")
# Source: 
# https://bioconductor.org/packages/release/bioc/html/DEGseq.html

# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

require(data.table)
require(ggplot2)
require(DEGseq)
# require(VennDiagram)
# require(gridExtra)

# Part I: Data----
load("data/lncap_setd7_normilized_annotated.RData")
dt1

# Make all gene names unique----
setkey(dt1, SYMBOL)
dt1[, N := 1:.N,
    by = SYMBOL]
dt1$SYMBOL[dt1$N > 1] <- paste(dt1$SYMBOL[dt1$N > 1],
                               dt1$N[dt1$N > 1],
                               sep = "_")
dt1[, N := NULL]
dt1

# Remove TNF samples
dt1 <- droplevels(subset(dt1,
                         select = -c(6, 9)))
dt1
summary(dt1)

# Remove all unnamed genes (LOC...)
gene.loc <- unique(dt1$SYMBOL[substr(dt1$SYMBOL, 1, 3) == "LOC"])
dt1 <- subset(dt1,
              !(SYMBOL %in% gene.loc))

# # Remove genes with low counts----
# summary(dt1[, -1])
# tmp <- rowSums(dt1[, -1])
# # Remove if total across 5 samples is no more than 20
# dt1 <- droplevels(subset(dt1,
#                          tmp > 20))
# dt1
# # 13,954 genes left, down from 24,421 genes

# DEGseq----
# WT_PEITC vs.WT====
DEGexp(geneExpMatrix1 = dt1[, c(3, 5:8)], 
       geneCol1 = 1, 
       expCol1 = 2, 
       
       geneExpMatrix2 = dt1[, c(3, 5:8)], 
       geneCol2 = 1, 
       expCol2 = 3,
       
       foldChange = 1,
       qValue = 0.1,
       thresholdKind = 5, 
       rawCount = FALSE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp")
  
# Read in the output file
dt.wtpeitc.wt <- fread("tmp/output_score.txt")
dt.wtpeitc.wt

# Write as CSV----
write.csv(dt.wtpeitc.wt,
          file = paste("tmp/lncap_DEGseq_wt_peitc-wt_.csv",
                       sep = ""),
          row.names = FALSE)

# KD_PEITC vs.KD----
DEGexp(geneExpMatrix1 = dt1[, c(3, 5:8)], 
       geneCol1 = 1, 
       expCol1 = 4, 
       
       geneExpMatrix2 = dt1[, c(3, 5:8)], 
       geneCol2 = 1, 
       expCol2 = 5,
       
       foldChange = 1,
       qValue = 0.1,
       thresholdKind = 5, 
       rawCount = FALSE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp")

# Read in the output file
dt.kdpeitc.kd <- fread("tmp/output_score.txt")
dt.kdpeitc.kd

# Write as CSV----
write.csv(dt.kdpeitc.kd,
          file = paste("tmp/lncap_DEGseq_kd_peitc-kd_.csv",
                       sep = ""),
          row.names = FALSE)

# # Significant genes----
# tmp <- dt2[dt2$`Signature(q-value(Storey et al. 2003) < 0.1)`, ]
# tmp
# write.csv(tmp,
#           file = "tmp/mes13_rnaseq_DEGseq2_signif_out.csv",
#           row.names = FALSE)

# sink()