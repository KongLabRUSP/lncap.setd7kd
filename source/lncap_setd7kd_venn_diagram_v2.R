# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCap/SETD7 KD experiment             |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/20/2017                                                              |
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

dt1$`LNCaP - Setd7` <- dt1$LNCaP - dt1$Setd7
dt1$`LNCaP PEITC - Setd7 PEITC` <- dt1$LNCaP_PEITC - dt1$Setd7_PEITC
dt1$`LNCaP TNF - Setd7 TNF` <- dt1$LNCaP_TNF - dt1$Setd7_TNF
dt1

# Genes with more than 2-fold change----
gene1.pls <- unique(dt1$SYMBOL[dt1$`LNCaP - Setd7` >= 1])
gene2.pls <- unique(dt1$SYMBOL[dt1$`LNCaP PEITC - Setd7 PEITC` >= 1])
gene3.pls <- unique(dt1$SYMBOL[dt1$`LNCaP TNF - Setd7 TNF` >= 1])

gene1.mns <- unique(dt1$SYMBOL[dt1$`LNCaP - Setd7` <= -1])
gene2.mns <- unique(dt1$SYMBOL[dt1$`LNCaP PEITC - Setd7 PEITC` <= -1])
gene3.mns <- unique(dt1$SYMBOL[dt1$`LNCaP TNF - Setd7 TNF` <= -1])

# CHECKPOINT----
subset(dt1,
       SYMBOL %in% gene1)

# Draw Venn diagrams----
p1 <- venn.diagram(x = list(`LNCaP - Setd7` = gene1.pls,
                            `LNCaP PEITC - Setd7 PEITC`= gene2.pls,
                            `LNCaP TNF - Setd7 TNF` = gene3.pls),
                   filename = NULL,
                   fill = c("light blue",
                            "grey", 
                            "light green"),
                   alpha = rep(0.5, 3),
                   compression = "lzw+p",
                   main = "Genes With At Least 2-Fold Positive Change")
grid.arrange(gTree(children = p1))

png(filename = "tmp/venn_lncap_setd7_pos.png",
    height = 5,
    width = 5,
    units = 'in',
    res = 300)
grid.arrange(gTree(children = p1))
graphics.off()

p2 <- venn.diagram(x = list(`LNCaP - Setd7` = gene1.mns,
                            `LNCaP PEITC - Setd7 PEITC`= gene2.mns,
                            `LNCaP TNF - Setd7 TNF` = gene3.mns),
                   filename = NULL,
                   fill = c("light blue",
                            "grey", 
                            "light green"),
                   alpha = rep(0.5, 3),
                   compression = "lzw+p",
                   main = "Genes With At Least 2-Fold Negative Change")
grid.arrange(gTree(children = p2))

png(filename = "tmp/venn_lncap_setd7_neg.png",
    height = 5,
    width = 5,
    units = 'in',
    res = 300)
grid.arrange(gTree(children = p2))
graphics.off()

# Save the table of genes----
dtt <- data.table(gene = unique(c(gene1.pls,
                                  gene2.pls,
                                  gene3.pls,
                                  gene1.mns,
                                  gene2.mns,
                                  gene3.mns)))
dtt <- merge(dtt, 
             data.table(gene = gene1.pls,
                        a1 = 1),
             by = "gene",
             all = TRUE)
dtt <- merge(dtt, 
             data.table(gene = gene2.pls,
                        a2 = 1),
             by = "gene",
             all = TRUE)
dtt <- merge(dtt, 
             data.table(gene = gene3.pls,
                        a3 = 1),
             by = "gene",
             all = TRUE)
dtt <- merge(dtt, 
             data.table(gene = gene1.mns,
                        b1 = 1),
             by = "gene",
             all = TRUE)
dtt <- merge(dtt, 
             data.table(gene = gene2.mns,
                        b2 = 1),
             by = "gene",
             all = TRUE)
dtt <- merge(dtt, 
             data.table(gene = gene3.mns,
                        b3 = 1),
             by = "gene",
             all = TRUE)
colnames(dtt) <- c("All log2>1 Genes",
                   "LNCaP > Setd7",
                   "LNCaP PEITC > Setd7 PEITC",
                   "LNCaP TNF > Setd7 TNF",
                   "LNCaP < Setd7",
                   "LNCaP PEITC < Setd7 PEITC",
                   "LNCaP TNF < Setd7 TNF")
dtt[is.na(dtt)] <- 0
colSums(dtt[, -1])
dtt
write.csv(dtt,
          file = "tmp/lncap_setd7_cahnged_genes.csv",
          row.names = FALSE)

# sink()