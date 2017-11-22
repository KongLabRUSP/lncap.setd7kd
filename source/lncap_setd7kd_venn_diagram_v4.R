# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCap/SETD7 KD experiment             |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/21/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_lncap_setd7kd_venn_diagram_v2.txt")

require(data.table)
require(ggplot2)
require(VennDiagram)
require(gridExtra)

# Part I: Data----
# Preprocessed data was created by 'lncap_setd7kd_data_preprocessing_v1.R'
load("tmp/lncap_setd7_normilized_annotated.RData")
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

# Rename columns----
colnames(dt1)[5:8] <- c("WT_PEITC",
                        "WT",
                        "KD_PEITC",
                        "KD")

# Differences----
dt1$`WT PEITC vs. WT` <- dt1$WT_PEITC - dt1$WT
dt1$`KD PEITC vs. KD` <- dt1$KD_PEITC - dt1$KD
dt1$`WT vs. KD` <- dt1$WT - dt1$KD
dt1

g1.1 <- dt1$SYMBOL[dt1$`WT PEITC vs. WT`>= 0.5]
g1.2 <- dt1$SYMBOL[dt1$`KD PEITC vs. KD` <= -0.5]

g1.3 <- dt1$SYMBOL[dt1$`WT PEITC vs. WT`<= -0.5]
g1.4 <- dt1$SYMBOL[dt1$`KD PEITC vs. KD` >= 0.5]

g2.1 <- dt1$SYMBOL[dt1$`KD PEITC vs. KD`>= 0.5]
g2.2 <- dt1$SYMBOL[dt1$`WT vs. KD` <= -0.5]

g2.3 <- dt1$SYMBOL[dt1$`KD PEITC vs. KD`<= -0.5]
g2.4 <- dt1$SYMBOL[dt1$`WT vs. KD` >= 0.5]

# Save the lists
dtt <- qpcR:::cbind.na(`WT PEITC vs. WT Up` = g1.1,
                       `WT PEITC vs. WT Down` = g1.2,
                       `WT PEITC vs. WT Down` = g1.3,
                       `WT PEITC vs. WT Up` = g1.4,
                       `KD PEITC vs. KD Up` = g2.1,
                       `WT vs. KD Down` = g2.2,
                       `KD PEITC vs. KD Down` = g2.3,
                       `WT vs. KD Up` = g2.4)
write.csv(dtt,
          file = "tmp/gene_list.csv",
          row.names = FALSE)

# Draw Venn diagrams: WT vs. KD----
png(filename = "tmp/venn_lncap_setd7_pos.png",
    height = 5,
    width = 5,
    units = 'in',
    res = 300)

venn.diagram(x = list(`WT PEITC vs. WT` = g1.1,
                      `KD PEITC vs. KD`= g1.2),
             col = c("green",
                     "red"),
             alpha = rep(0.5, 2),
             compression = "lzw+p",
             cat.just = list(c(0, 0),
                             c(1, 0)),
             main = "",
             filename = "tmp/venn_wt.up_vs_kd.dn.png",
             units = 'in',
             height = 5,
             width = 5)

venn.diagram(x = list(`WT PEITC vs. WT` = g1.3,
                      `KD PEITC vs. KD`= g1.4),
             col = c("red",
                      "green"),
             alpha = rep(0.5, 2),
             compression = "lzw+p",
             cat.just = list(c(1, -1),
                             c(0, -1)),
             main = "",
             filename = "tmp/venn_wt.dn_vs_kd.up.png",
             units = 'in',
             height = 5,
             width = 5)

# Draw Venn diagrams: All vs. KD----
venn.diagram(x = list(`KD PEITC vs. KD` = g2.1,
                      `WT vs. KD`= g2.2),
             col = c("green",
                     "red"),
             alpha = rep(0.5, 2),
             compression = "lzw+p",
             cat.just = list(c(0, -5),
                             c(2, -5)),
             main = "",
             filename = "tmp/venn_ctr.up_vs_trt.dn.png",
             units = 'in',
             height = 5,
             width = 5)

venn.diagram(x = list(`KD PEITC vs. KD` = g2.3,
                      `WT vs. KD`= g2.4),
             col = c("red",
                      "green"),
             alpha = rep(0.5, 2),
             compression = "lzw+p",
             cat.just = list(c(1, 0),
                             c(-1, 0)),
             main = "",
             filename = "tmp/venn_ctr.dn_vs_trt.up.png",
             units = 'in',
             height = 5,
             width = 5)

# sink()