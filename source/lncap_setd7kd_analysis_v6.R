# |----------------------------------------------------------------------------------|
# | Project: Study of LNCaP WT/KD cells                                              |
# | Script: Analysis of microarray data from a LNCaP WT/KD experiment                |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 12/02/2017                                                              |
# | Modified: 12/29/2017, Chao's request to flip WT and KD in the comparison         |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
sink(file = "tmp/log_lncap_setd7kd_analysis_v6.txt")

require(data.table)
require(ggplot2)
require(VennDiagram)
require(gridExtra)

# Part I: Data----
# Preprocessed data was created by 'WT_KDKD_data_preprocessing_v2.R'
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

# Part II: hitmaps
# Differences----
dt1$`WT PEITC vs. WT` <- dt1$WT_PEITC - dt1$WT
dt1$`KD PEITC vs. KD` <- dt1$KD_PEITC - dt1$KD
dt1$`KD vs. WT` <- dt1$KD - dt1$WT

# Means----
dt1$`Mean(WT PEITC, WT)` <- (dt1$WT_PEITC + dt1$WT)/2
dt1$`Mean(KD PEITC, KD)` <- (dt1$KD_PEITC + dt1$KD)/2
dt1$`Mean(KD, WT)` <- (dt1$KD + dt1$WT)/2

dt1

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
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))
p01

tiff(filename = "tmp/all_anno_genes_diffs_wt.vs.kd.tiff",
     height = 9,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p01)
graphics.off() 

# Heatmap: all 4 groups----
tmp <- droplevels(subset(dt1,
                         SYMBOL %in% unique(peitc.nopeitc$SYMBOL)))
dtl <- melt.data.table(tmp,
                       id.vars = "SYMBOL",
                       measure.vars = 5:8,
                       variable.name = "Group",
                       value.name = "Readout")
dtl$SYMBOL <- factor(dtl$SYMBOL,
                     levels = gene.keep)
dtl

# Plot all annotated genes found in all 4 samples----
p01.1 <- ggplot(data = dtl) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = Readout),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = range(dtl$Readout),
                       name = "Readout") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p01.1

tiff(filename = "tmp/all_anno_genes_wt.vs.kd.tiff",
     height = 9,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p01.1)
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
                     (`KD PEITC vs. KD`>= 0.75 &
                        `KD vs. WT` <= -0.75)  |
                       (`KD PEITC vs. KD` <= -0.75 &
                          `KD vs. WT` >= 0.75))

kdpeits.kd <- kdpeits.kd[, c("SYMBOL", 
                             "from.kd")]
kdpeits.kd
kdpeits.kd <- kdpeits.kd[order(kdpeits.kd$from.kd), ]
ordered.genes <- as.character(kdpeits.kd$SYMBOL)

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
kdpeits.kd$`Diff of Logs`[kdpeits.kd$`Diff of Logs` > 2] <- 2
kdpeits.kd$`Diff of Logs`[kdpeits.kd$`Diff of Logs` < -2] <- -2

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
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))
p02

tiff(filename = "tmp/all_anno_genes_diffs_from.kd.tiff",
     height = 9,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p02)
graphics.off() 

# Heatmap: all 4 groups----
tmp <- droplevels(subset(dt1,
                         SYMBOL %in% unique(kdpeits.kd$SYMBOL)))
dtl <- melt.data.table(tmp,
                       id.vars = "SYMBOL",
                       measure.vars = 5:8,
                       variable.name = "Group",
                       value.name = "Readout")
dtl$SYMBOL <- factor(dtl$SYMBOL,
                     levels = gene.keep)
dtl

# Plot all annotated genes found in all 4 samples----
p02.1 <- ggplot(data = dtl) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = Readout),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = range(dtl$Readout),
                       name = "Readout") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p02.1

tiff(filename = "tmp/all_anno_genes_from.kd.tiff",
     height = 9,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p02.1)
graphics.off()

# Save values of KD, WT and KD PEICT (Chao, 01/03/2018)----
dt.vals <- droplevels(subset(dt1,
                             SYMBOL %in% kdpeits.kd$SYMBOL))

write.csv(dt.vals,
          file = "tmp/KD vs WT and KD PEICT Genes and Values.csv")

# # Reverse direction KD vs. WT----
# kdpeits.kd$`Diff of Logs`[kdpeits.kd$Group == "KD vs. WT"] <- (-1)*kdpeits.kd$`Diff of Logs`[kdpeits.kd$Group == "KD vs. WT"]
# levels(kdpeits.kd$Group)[1] <- "WT vs. KD"
# 
# p03 <- ggplot(data = kdpeits.kd) +
#   geom_tile(aes(x =  Group,
#                 y = SYMBOL,
#                 fill = `Diff of Logs`)) +
#   scale_fill_gradient2(high = "green",
#                        mid = "black",
#                        low = "red",
#                        # limit = c(-2, 2),
#                        name = "Diff of Logs") +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete("Gene",
#                    expand = c(0, 0)) +
#   ggtitle("Sorted by Most Differences in LNCaP\nChanges From KD") +
#   theme(plot.title = element_text(hjust = 0.5))
# p03
# 
# tiff(filename = "tmp/all_anno_genes_diffs_from.kd_rev.tiff",
#      height = 9,
#      width = 7,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p03)
# graphics.off() 

# Plot all gene expressions
dt1[, 5:8]
sqrt(nrow(dt1))
179*178 - nrow(dt1)

# Add 58 rows of NAs to make a 179*178 matrix
tmp <- data.table(matrix(NA,
                         nrow = 179*178 - nrow(dt1),
                         ncol = 4))
colnames(tmp) <- colnames(dt1)[5:8]
tmp
tmp <- rbindlist(list(dt1[, 5:8],
                      tmp))
tmp$rownum <- rev(rep(1:179, 178))
tmp$colnum <- rep(1:178, each = 179)
tmp <- melt.data.table(tmp,
                       id.vars = c("rownum",
                                   "colnum"),
                       measure.vars = 1:4)
tmp

p04<- ggplot(data = tmp) +
  facet_wrap(~ variable) +
  geom_tile(aes(x = colnum,
                y = rownum,
                fill = value)) +
  scale_fill_gradient2(high = "red",
                       limit = range(tmp$value,
                                     na.rm = TRUE),
                       name = "log2(expr)") +
  scale_x_discrete("",
                   expand = c(0, 0)) +
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("Normalized Expressions of 31,804 Genes") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        plot.title = element_text(hjust = 0.5))
p04

tiff(filename = "tmp/all_anno_genes_norm_expr.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p04)
graphics.off()

# Plot 150 genes with teh most differential expressions (i.e. variability)
# a. highest range
dt1[, 5:8]
dt1$rng.diff <- apply(X = dt1[, 5:8],
                      MARGIN = 1,
                      FUN = function(a) {
                        return(diff(range(a, 
                                          na.rm = TRUE)))
                      })
hist(log(dt1$rng.diff), 
     100,
     xlab = "Log of Range Difference",
     main = "Distribution of Range Differences")

dt1$std <- apply(X = dt1[, 5:8],
                 MARGIN = 1,
                 FUN = function(a) {
                   return(sd(a, 
                             na.rm = TRUE))
                 })
hist(log(dt1$std), 
     100,
     xlab = "Log of SD",
     main = "Distribution of SD")

plot(log(dt1$rng.diff) ~ log(dt1$std))

# Top 150 genes (by SD)
ndx <- which(rank(dt1$std) > (nrow(dt1) - 150))
ndx

tmp <- dt1[ndx, ]
tmp
summary(tmp)
hist(dt1$std, 
     100,
     xlab = "SD",
     main = "Distribution of SD")
abline(v = min(tmp$std),
       lty = 2)
arrows(x0 = min(tmp$std),
       y0 = 1500,
       x1 = max(tmp$std),
       y1 = 1500)

# Hitmap----
dtl <- melt.data.table(tmp,
                       id.vars = "SYMBOL",
                       measure.vars = 5:8,
                       variable.name = "Group",
                       value.name = "Readout")
dtl$SYMBOL <- factor(dtl$SYMBOL)
dtl

dtl$Readout <- scale(dtl$Readout)

# Plot all annotated genes found in all 4 samples----
p05 <- ggplot(data = dtl) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = Readout)) +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high = "blue",
                       limit = range(dtl$Readout),
                       name = "Normalized Expr") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p05

tiff(filename = "tmp/top_150_most_diff_expr.tiff",
     height = 9,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p05)
graphics.off()

sink()