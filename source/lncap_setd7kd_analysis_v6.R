# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCaP WT/KD experiment                |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 12/02/2017                                                              |
# |----------------------------------------------------------------------------------|
# NOTE: remove TNF samples
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_WT_KDKD_analysis_v5.txt")

require(data.table)
require(ggplot2)
require(VennDiagram)
require(gridExtra)

# Part I: Data----
# Preprocessed data was created by 'WT_KDKD_data_preprocessing_v1.R'
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

# Part II: p-Values for single-replica samples----
# a. WT PEITC vs. WT----
dt1$`Diff(WT PEITC, WT)` <- dt1$WT_PEITC - dt1$WT
dt1$`Mean(WT PEITC, WT)` <- (dt1$WT_PEITC + dt1$WT)/2

tiff(filename = "tmp/lncap_marray_wtpeitc_vs_wt_diff_mean.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt1$`Diff(WT PEITC, WT)` ~ dt1$`Mean(WT PEITC, WT)`,
     xlab = "Means",
     ylab = "Differences",
     main = "LNCaP Microarray: WT PEITC vs. WT")
graphics.off()

# Regularize by the above within margin of epsilon
epln <- 0.1
dt1$wtpeitc.wt.sd <- NA

for (i in 1:nrow(dt1)) {
  tmp <- subset(dt1,
                (`Mean(WT PEITC, WT)` <= dt1$`Mean(KD PEITC, KD)`[i] + epln) &
                  (`Mean(WT PEITC, WT)` >= dt1$`Mean(KD PEITC, KD)`[i] - epln))
  dt1$wtpeitc.wt.sd[i] <- sd(tmp$`Diff(WT PEITC, WT)`)
}

hist(dt1$`Diff(WT PEITC, WT)`, 100)

m1 <- loess(dt1$wtpeitc.wt.sd ~ dt1$`Mean(WT PEITC, WT)`)
dt1$wtpeitc.wt.sd.fit <- predict(m1,
                                 newdata = data.frame(`Mean(WT PEITC, WT)` = dt1$`Mean(WT PEITC, WT)`))
dt1
dt1 <- dt1[order(dt1$`Mean(WT PEITC, WT)`), ]

plot(dt1$sd ~ dt1$`Mean(KD PEITC, KD)`,
     xlab = "Means",
     ylab = "SD",
     main = "WT PEITC vs. WT")
lines(dt1$wtpeitc.wt.sd.fit ~ dt1$`Mean(WT PEITC, WT)`,
      col = "red",
      lw = 2)

# Assuming t follows normal distribution (it does not!)
# dt1$t <- dt1$`WT PEITC vs. WT`/dt1$sd
# qqnorm(dt1$t)
# abline(0, 1)

dt1$t <- dt1$`WT PEITC vs. WT`/dt1$`Fitted SD`
qqnorm(dt1$t)
abline(0, 1)

dt1$p <- 2*pnorm(-abs(dt1$t))
hist(dt1$p)

tmp <- subset(dt1,
              p <= 0.05,
              select = c(1:6, 9, 12, 17:20))
tmp <- tmp[order(tmp$SYMBOL), ]
tmp
write.csv(tmp, file = "tmp/pvals.csv")





# # Part II: long data----
# dtl <- melt.data.table(dt1,
#                        id.vars = "SYMBOL",
#                        measure.vars = 5:8,
#                        variable.name = "Group",
#                        value.name = "Expression")
# dtl$SYMBOL <- factor(dtl$SYMBOL)
# dtl$Group <- factor(dtl$Group,
#                     levels = c("WT",
#                                "KD",
#                                "WT_PEITC",
#                                "KD_PEITC"))
# dtl
# summary(dtl)


# Hitmap of differences----
# Differences----
dt1$`WT PEITC vs. WT` <- dt1$WT_PEITC - dt1$WT
dt1$`KD PEITC vs. KD` <- dt1$KD_PEITC - dt1$KD
dt1$`WT vs. KD` <- dt1$WT - dt1$KD

# Means----
dt1$`Mean(WT PEITC, WT)` <- (dt1$WT_PEITC + dt1$WT)/2
dt1$`Mean(KD PEITC, KD)` <- (dt1$KD_PEITC + dt1$KD)/2
dt1$`Mean(WT, KD)` <- (dt1$WT + dt1$KD)/2

dt1

# a. Diffs in WT vs. diffs in KD----
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
                        (`WT PEITC vs. WT`>= 0.5 &
                           `KD PEITC vs. KD` <= -0.5)  |
                          (`WT PEITC vs. WT` <= -0.5 &
                             `KD PEITC vs. KD` >= 0.5))
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
                       limit = c(-2, 2),
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
                        levels = c("WT vs. KD",
                                   "KD PEITC vs. KD"))
dt.diff
length(unique(dt.diff$SYMBOL))

dt1$from.kd <- dt1$`WT vs. KD` - dt1$`KD PEITC vs. KD`
hist(dt1$from.kd, 100)

kdpeits.kd <- subset(dt1,
                     (`KD PEITC vs. KD`>= 0.5 &
                        `WT vs. KD` <= -0.5)  |
                       (`KD PEITC vs. KD` <= -0.5 &
                          `WT vs. KD` >= 0.5))

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
# Reset upper bound to 2
kdpeits.kd$`Diff of Logs`[kdpeits.kd$`Diff of Logs` > 2] <- 2

p02 <- ggplot(data = kdpeits.kd) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = `Diff of Logs`)) +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       limit = c(-2, 2),
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











# p-Values for single-replica samples----
# Plot differences vs. means in WT PEITC vs. WT----
plot(dt1$`WT PEITC vs. WT` ~ dt1$`Mean(WT PEITC, WT)`,
     xlab = "Means",
     ylab = "Differences",
     main = "WT PEITC vs. WT")
# Regularize by the above within margin of epsilon
epln <- 0.1
dt1$sd <- NA

for (i in 1:nrow(dt1)) {
  tmp <- subset(dt1,
                (`Mean(WT PEITC, WT)` <= dt1$`Mean(KD PEITC, KD)`[i] + epln) &
                  (`Mean(WT PEITC, WT)` >= dt1$`Mean(KD PEITC, KD)`[i] - epln))
  dt1$sd[i] <- sd(tmp$`WT PEITC vs. WT`)
}

hist(dt1$`WT PEITC vs. WT`, 100)

m1 <- loess(dt1$sd ~ dt1$`Mean(WT PEITC, WT)`)
dt1$`Fitted SD` <- predict(m1,
                           newdata = data.frame(`Mean(WT PEITC, WT)` = dt1$`Mean(WT PEITC, WT)`))
dt1
dt1 <- dt1[order(dt1$`Mean(WT PEITC, WT)`), ]

plot(dt1$sd ~ dt1$`Mean(KD PEITC, KD)`,
     xlab = "Means",
     ylab = "SD",
     main = "WT PEITC vs. WT")
lines(dt1$`Fitted SD` ~ dt1$`Mean(WT PEITC, WT)`,
      col = "red",
      lw = 2)

# Assuming t follows normal distribution (it does not!)
# dt1$t <- dt1$`WT PEITC vs. WT`/dt1$sd
# qqnorm(dt1$t)
# abline(0, 1)

dt1$t <- dt1$`WT PEITC vs. WT`/dt1$`Fitted SD`
qqnorm(dt1$t)
abline(0, 1)

dt1$p <- 2*pnorm(-abs(dt1$t))
hist(dt1$p)

tmp <- subset(dt1,
              p <= 0.05,
              select = c(1:6, 9, 12, 17:20))
tmp <- tmp[order(tmp$SYMBOL), ]
tmp
write.csv(tmp, file = "tmp/pvals.csv")

# sink()