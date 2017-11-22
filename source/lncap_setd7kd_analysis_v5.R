# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCaP WT/KD experiment                |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/21/2017                                                              |
# |----------------------------------------------------------------------------------|
# NOTE: remove TNF samples
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_WT_KDKD_analysis_v5.txt")

require(data.table)
require(ggplot2)

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

# Part II: long data----
dtl <- melt.data.table(dt1,
                       id.vars = "SYMBOL",
                       measure.vars = 5:8,
                       variable.name = "Group",
                       value.name = "Expression")
dtl$SYMBOL <- factor(dtl$SYMBOL)
dtl$Group <- factor(dtl$Group,
                    levels = c("WT",
                               "KD",
                               "WT_PEITC",
                               "KD_PEITC"))
dtl
summary(dtl)

# Part III: compare cellines----
# Hitmap of differences----
dt1$`WT PEITC vs. WT` <- dt1$WT_PEITC - dt1$WT
dt1$`KD PEITC vs. KD` <- dt1$KD_PEITC - dt1$KD
dt1$`WT vs. KD` <- dt1$WT - dt1$KD
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

tmp <- subset(dt1,
              (`WT PEITC vs. WT`>= 0.5 &
                 `KD PEITC vs. KD` <= -0.5)  |
                (`WT PEITC vs. WT` <= -0.5 &
                   `KD PEITC vs. KD` >= 0.5))
summary(tmp)

tmp <- tmp[, c("SYMBOL", 
               "wt.vs.kd")]
tmp
tmp <- tmp[order(tmp$wt.vs.kd), ]
ordered.genes <- as.character(tmp$SYMBOL)

gene.keep <- ordered.genes[c(1:20,
                             (length(ordered.genes) - 19):length(ordered.genes))]
# gene.keep <- ordered.genes

tmp <- subset(dt.diff,
              SYMBOL %in% gene.keep)

tmp$SYMBOL <- factor(tmp$SYMBOL,
                         levels = gene.keep)
tmp
summary(tmp)

# Plot----
p01 <- ggplot(data = tmp) +
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
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p01)
graphics.off() 

# b. Diffs form KD----
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

tmp <- subset(dt1,
              (`KD PEITC vs. KD`>= 0.5 &
                 `WT vs. KD` <= -0.5)  |
                (`KD PEITC vs. KD` <= -0.5 &
                   `WT vs. KD` >= 0.5))

tmp <- tmp[, c("SYMBOL", 
               "from.kd")]
tmp
tmp <- tmp[order(tmp$from.kd), ]
ordered.genes <- as.character(tmp$SYMBOL)

gene.keep <- ordered.genes[c(1:20,
                             (length(ordered.genes) - 19):length(ordered.genes))]
# gene.keep <- ordered.genes

tmp <- subset(dt.diff,
              SYMBOL %in% gene.keep)

tmp$SYMBOL <- factor(tmp$SYMBOL,
                     levels = gene.keep)
tmp

# Plot----
summary(tmp$`Diff of Logs`)
# Reset upper bound to 2
tmp$`Diff of Logs`[tmp$`Diff of Logs` > 2] <- 2

p02 <- ggplot(data = tmp) +
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
     height = 10,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p02)
graphics.off() 

# sink()