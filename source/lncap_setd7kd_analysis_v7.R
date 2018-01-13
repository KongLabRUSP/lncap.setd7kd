# |----------------------------------------------------------------------------------|
# | Project: Study of LNCaP WT/KD cells                                              |
# | Script: Analysis of microarray data from a LNCaP WT/KD experiment                |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 12/02/2017                                                              |
# | Modified: 12/29/2017, Chao's request to flip WT and KD in the comparison         |
# |----------------------------------------------------------------------------------|
# NOTE: remove TNF samples
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_lncap_setd7kd_analysis_v6.txt")

require(data.table)
require(ggplot2)
require(VennDiagram)
require(gridExtra)

# Part I: Data----
# Preprocessed data was created by 'WT_KDKD_data_preprocessing_v1.R'
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

# Rename columns----
colnames(dt1)[5:8] <- c("WT_PEITC",
                        "WT",
                        "KD_PEITC",
                        "KD")
summary(dt1)

# Remove all unnamed genes (LOC...)
gene.loc <- unique(dt1$SYMBOL[substr(dt1$SYMBOL, 1, 3) == "LOC"])
dt1 <- subset(dt1,
              !(SYMBOL %in% gene.loc))

# Part II: hitmaps
# Differences----
dt1$`WT PEITC vs. WT` <- log2(dt1$WT_PEITC/dt1$WT)
hist(dt1$`WT PEITC vs. WT`, 100)

dt1$`KD PEITC vs. KD` <- log2(dt1$KD_PEITC/dt1$KD)
hist(dt1$`KD PEITC vs. KD`, 100)

dt1$`KD vs. WT` <- log2(dt1$KD/dt1$WT)
hist(dt1$`KD vs. WT`, 100)

# Means----
dt1$`Mean(WT PEITC, WT)` <- (dt1$WT_PEITC + dt1$WT)/2
dt1$`Mean(KD PEITC, KD)` <- (dt1$KD_PEITC + dt1$KD)/2
dt1$`Mean(KD, WT)` <- (dt1$KD + dt1$WT)/2

dt1

# Part II: SD Estimation----
# a. WT PEITC vs. WT----
plot(dt1$`WT PEITC vs. WT` ~ dt1$`Mean(WT PEITC, WT)`)
dt1$abs.wtpeitc.wt <- abs(dt1$`WT PEITC vs. WT`)


m1 <- smooth.spline(x = dt1$`Mean(WT PEITC, WT)`,
                    y = dt1$abs.wtpeitc.wt)
plot(dt1$abs.wtpeitc.wt ~ dt1$`Mean(WT PEITC, WT)`,
     main = "SD Estimation with Smoothing Spline",
     xlab = "Means",
     ylab = "Absolue Differences")
lines(m1, 
      col = "red",
      lw = 3)

res1 <- data.table(`Mean(WT PEITC, WT)` = m1$x,
                                `SD(WT PEITC, WT)` = m1$y)
res1

dt1 <- merge(dt1,
             res1,
             by = "Mean(WT PEITC, WT)")
dt1$`t(WT PEITC, WT)` <- dt1$`WT PEITC vs. WT`/dt1$`SD(WT PEITC, WT)`
hist(dt1$`t(WT PEITC, WT)`, 100)
qqnorm(dt1$`t(WT PEITC, WT)`)

hist(abs(dt1$`t(WT PEITC, WT)`), 100)
# 5% cut-off
cutoff <- quantile(abs(dt1$`t(WT PEITC, WT)`),
                   probs = 0.975)
cutoff 
dt1$`p(WT PEITC, WT)` <- abs(dt1$`t(WT PEITC, WT)`) >= cutoff
table(dt1$`p(WT PEITC, WT)`)

# b. KD PEITC vs. KD----
plot(dt1$`KD PEITC vs. KD` ~ dt1$`Mean(KD PEITC, KD)`)
dt1$abs.kdpeitc.kd <- abs(dt1$`KD PEITC vs. KD`)

m2 <- smooth.spline(x = dt1$`Mean(KD PEITC, KD)`,
                    y = dt1$abs.kdpeitc.kd)
plot(dt1$abs.kdpeitc.kd ~ dt1$`Mean(KD PEITC, KD)`,
     main = "SD Estimation with Smoothing Spline",
     xlab = "Means",
     ylab = "Absolue Differences")
lines(m2, 
      col = "red",
      lw = 3)

res2 <- data.table(`Mean(KD PEITC, KD)` = m2$x,
                   `SD(KD PEITC, KD)` = m2$y)
res2

dt1 <- merge(dt1,
             res2,
             by = "Mean(KD PEITC, KD)")
dt1$`t(KD PEITC, KD)` <- dt1$`KD PEITC vs. KD`/dt1$`SD(KD PEITC, KD)`
hist(dt1$`t(KD PEITC, KD)`, 100)
qqnorm(dt1$`t(KD PEITC, KD)`)

hist(abs(dt1$`t(KD PEITC, KD)`), 100)
# 5% cut-off
cutoff <- quantile(abs(dt1$`t(KD PEITC, KD)`),
                   probs = 0.975)
cutoff 
dt1$`p(KD PEITC, KD)` <- abs(dt1$`t(KD PEITC, KD)`) >= cutoff
table(dt1$`p(KD PEITC, KD)`)

# c. KD vs. WT----
plot(dt1$`KD vs. WT` ~ dt1$`Mean(KD, WT)`)
dt1$abs.kd.wt <- abs(dt1$`KD vs. WT`)

m3 <- smooth.spline(x = dt1$`Mean(KD, WT)`,
                    y = dt1$abs.kd.wt)
plot(dt1$abs.kd.wt ~ dt1$`Mean(KD, WT)`,
     main = "SD Estimation with Smoothing Spline",
     xlab = "Means",
     ylab = "Absolue Differences")
lines(m3, 
      col = "red",
      lw = 3)

res3 <- data.table(`Mean(KD, WT)` = m3$x,
                   `SD(KD, WT)` = m3$y)
res3

dt1 <- merge(dt1,
             res3,
             by = "Mean(KD, WT)")
dt1$`t(KD, WT)` <- dt1$`KD vs. WT`/dt1$`SD(KD, WT)`
hist(dt1$`t(KD, WT)`, 100)
qqnorm(dt1$`t(KD, WT)`)

hist(abs(dt1$`t(KD, WT)`), 100)
# 5% cut-off
cutoff <- quantile(abs(dt1$`t(KD, WT)`),
                   probs = 0.975)
cutoff 
dt1$`p(KD, WT)` <- abs(dt1$`t(KD, WT)`) >= cutoff
table(dt1$`p(KD, WT)`)

# Part III: differences----
# a. Diffs in KD vs. diffs in WT----
dt.diff <- melt.data.table(dt1,
                           id.vars = "SYMBOL",
                           measure.vars = 9:10,
                           variable.name = "Group",
                           value.name = "Diff of Logs")
dt.diff$SYMBOL <- factor(dt.diff$SYMBOL)
dt.diff$Group <- factor(dt.diff$Group,
                        levels = c("WT PEITC vs. WT",
                                   "KD PEITC vs. KD"))
dt.diff
length(unique(dt.diff$SYMBOL))

# Sort by differences of differences----
dt1$wt.vs.kd <- dt1$`WT PEITC vs. WT` - dt1$`KD PEITC vs. KD`
hist(dt1$wt.vs.kd, 100)

peitc.nopeitc <- subset(dt1,
                        (`WT PEITC vs. WT`>= 0.3 &
                           `KD PEITC vs. KD` <= -0.3)  |
                          (`WT PEITC vs. WT` <= -0.3 &
                             `KD PEITC vs. KD` >= 0.3))
summary(peitc.nopeitc)

peitc.nopeitc <- peitc.nopeitc[, c("SYMBOL", 
                                   "wt.vs.kd")]
peitc.nopeitc
peitc.nopeitc <- peitc.nopeitc[order(peitc.nopeitc$wt.vs.kd), ]
ordered.genes <- as.character(peitc.nopeitc$SYMBOL)

# gene.keep <- ordered.genes[c(1:20,
#                              (length(ordered.genes) - 19):length(ordered.genes))]

# Get all genes----
gene.keep <- ordered.genes

peitc.nopeitc <- subset(dt.diff,
                        SYMBOL %in% gene.keep)

peitc.nopeitc$SYMBOL <- factor(peitc.nopeitc$SYMBOL,
                               levels = gene.keep)
peitc.nopeitc
summary(peitc.nopeitc)

# Plot----
p01 <- ggplot(data = peitc.nopeitc) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = `Diff of Logs`)) +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       # limit = c(-2, 2),
                       name = "Diff of Logs") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Sorted by Most Differences in LNCaP\nChanges in WT vs. Changes in KD") +
  theme(plot.title = element_text(hjust = 0.5))
p01

tiff(filename = "tmp/all_anno_genes_diffs_wt.vs.kd.tiff",
     height = 9,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p01)
graphics.off() 

# b. Diffs from KD----
dt.diff <- melt.data.table(dt1,
                           id.vars = "SYMBOL",
                           measure.vars = 10:11,
                           variable.name = "Group",
                           value.name = "Diff of Logs")
dt.diff$SYMBOL <- factor(dt.diff$SYMBOL)
dt.diff$Group <- factor(dt.diff$Group,
                        levels = c("KD vs. WT",
                                   "KD PEITC vs. KD"))
dt.diff
length(unique(dt.diff$SYMBOL))

dt1$from.kd <- dt1$`KD vs. WT` - dt1$`KD PEITC vs. KD`
hist(dt1$from.kd, 100)

kdpeits.kd <- subset(dt1,
                     (`KD PEITC vs. KD`>= 0.4 &
                        `KD vs. WT` <= -0.4)  |
                       (`KD PEITC vs. KD` <= -0.4 &
                          `KD vs. WT` >= 0.4))

kdpeits.kd <- kdpeits.kd[, c("SYMBOL", 
                             "from.kd")]
kdpeits.kd
kdpeits.kd <- kdpeits.kd[order(kdpeits.kd$from.kd), ]
ordered.genes <- as.character(kdpeits.kd$SYMBOL)

# gene.keep <- ordered.genes[c(1:20,
#                              (length(ordered.genes) - 19):length(ordered.genes))]

# All genes
gene.keep <- ordered.genes

kdpeits.kd <- subset(dt.diff,
                     SYMBOL %in% gene.keep)

kdpeits.kd$SYMBOL <- factor(kdpeits.kd$SYMBOL,
                            levels = gene.keep)
kdpeits.kd

# Plot----
summary(kdpeits.kd$`Diff of Logs`)
# Reset bounds to +/-2
# kdpeits.kd$`Diff of Logs`[kdpeits.kd$`Diff of Logs` > 2] <- 2
# kdpeits.kd$`Diff of Logs`[kdpeits.kd$`Diff of Logs` < -2] <- -2

p02 <- ggplot(data = kdpeits.kd) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = `Diff of Logs`)) +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       # limit = c(-2, 2),
                       name = "Diff of Logs") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Sorted by Most Differences in LNCaP\nChanges From KD") +
  theme(plot.title = element_text(hjust = 0.5))
p02

tiff(filename = "tmp/all_anno_genes_diffs_from.kd.tiff",
     height = 9,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p02)
graphics.off() 

# Save values of KD, WT and KD PEICT (Chao, 01/03/2018)----
dt.vals <- droplevels(subset(dt1,
                             SYMBOL %in% kdpeits.kd$SYMBOL))

write.csv(dt.vals,
          file = "tmp/KD vs WT and KD PEICT Genes and Values.csv")

# Part IV: genes selected for qPCR----
dt.val <- dt1[dt1$SYMBOL %in% c("LRRTM1",
                                "FCRL6", 
                                "THSD7A",
                                "APOM", 
                                "ARMC12",
                                "TKTL1", 
                                "ZNF235",
                                "TPM3",
                                "OR9A2", 
                                "DCAF4L1",
                                "CBLN2",
                                "TAS2R19"), ]
tt1 <- data.table(gene = dt.val$SYMBOL,
                  `WT PEITC - WT` = dt.val$WT_PEITC - dt.val$WT,
                  `WT PEITC vs. WT` = 2^dt.val$`WT PEITC vs. WT`,
                  `p(WT PEITC, WT)` = dt.val$`p(WT PEITC, WT)`,
                  `KD PEITC - KD` = dt.val$KD_PEITC - dt.val$KD,
                  `KD PEITC vs. KD` = 2^dt.val$`KD PEITC vs. KD`,
                  `p(KD PEITC, KD)` = dt.val$`p(KD PEITC, KD)`,
                  `KD - WT` = dt.val$KD - dt.val$WT,
                  `KD vs. WT` = 2^dt.val$`KD vs. WT`,
                  `p(KD, WT)` = dt.val$`p(KD, WT)`)
tt1
write.csv(tt1,
          file = "tmp/qPCR Genes.csv")
# sink()