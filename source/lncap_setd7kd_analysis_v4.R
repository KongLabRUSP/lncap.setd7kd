# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCap/SETD7 KD experiment             |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/20/2017                                                              |
# |----------------------------------------------------------------------------------|
# NOTE: remove TNF samples
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_lncap_setd7kd_analysis_v4.txt")

require(data.table)
require(ggplot2)

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

# Part II: long data----
dtl <- melt.data.table(dt1,
                       id.vars = "SYMBOL",
                       measure.vars = 5:8,
                       variable.name = "Group",
                       value.name = "Expression")
dtl$SYMBOL <- factor(dtl$SYMBOL)
dtl$Group <- factor(dtl$Group,
                    levels = c("LNCaP",
                               "Setd7",
                               "LNCaP_PEITC",
                               "Setd7_PEITC"))
dtl
summary(dtl)

# Part III: compare cellines----
# Hitmap of differences----
dt1$`LNCaP - Setd7` <- dt1$LNCaP - dt1$Setd7
dt1$`LNCaP PEITC - Setd7 PEITC` <- dt1$LNCaP_PEITC - dt1$Setd7_PEITC
dt1$`Setd7 PEITC - Setd7` <- dt1$Setd7_PEITC - dt1$Setd7
dt1

dt.diff <- melt.data.table(dt1,
                           id.vars = "SYMBOL",
                           measure.vars = 9:11,
                           variable.name = "Group",
                           value.name = "Diff of Logs")
dt.diff$SYMBOL <- factor(dt.diff$SYMBOL)
dt.diff$Group <- factor(dt.diff$Group,
                        levels = c("LNCaP - Setd7",
                                   "LNCaP PEITC - Setd7 PEITC",
                                   "Setd7 PEITC - Setd7"))
dt.diff
length(unique(dt.diff$SYMBOL))

# Sort by differences in controls----
tmp <- subset(dt.diff,
              Group == "LNCaP - Setd7")
tmp <- tmp[order(tmp$`Diff of Logs`), ]
ordered.genes <- as.character(tmp$SYMBOL)

gene.keep <- ordered.genes[c(1:20,
                             (length(ordered.genes) - 19):length(ordered.genes))]

tmp <- subset(dt.diff,
              SYMBOL %in% gene.keep)

tmp$SYMBOL <- factor(tmp$SYMBOL,
                         levels = gene.keep)
tmp

# Plot----
p01 <- ggplot(data = tmp) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = `Diff of Logs`),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       limit = range(dt.diff$`Diff of Logs`),
                       name = "Diff of Logs") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Sorted by Most Changed in LNCaP vs. Setd7") +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p01

tiff(filename = "tmp/all_anno_genes_diffs_1.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p01)
graphics.off() 

# Plot all 4 groups for these genes----
tmp <- subset(dtl,
              SYMBOL %in% gene.keep)
tmp

p11 <- ggplot(data = tmp) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = Expression),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = range(dtl$Expression),
                       name = "Expression") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Genes Most Changed in LNCaP vs. Setd7") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p11

tiff(filename = "tmp/all_anno_genes_by_dctrl.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p11)
graphics.off()

# Sort by differences in treated----
tmp <- subset(dt.diff,
              Group == "LNCaP PEITC - Setd7 PEITC")
tmp <- tmp[order(tmp$`Diff of Logs`), ]
ordered.genes <- as.character(tmp$SYMBOL)

gene.keep <- ordered.genes[c(1:20,
                             (length(ordered.genes) - 19):length(ordered.genes))]

tmp <- subset(dt.diff,
              SYMBOL %in% gene.keep)

tmp$SYMBOL <- factor(tmp$SYMBOL,
                     levels = gene.keep)
tmp

# Plot----
p02 <- ggplot(data = tmp) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = `Diff of Logs`),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       limit = range(dt.diff$`Diff of Logs`),
                       name = "Diff of Logs") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Sorted by Most Changed in LNCaP PEITC vs. Setd7 PEITC") +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p02

tiff(filename = "tmp/all_anno_genes_diffs_2.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p02)
graphics.off()

# Plot all 4 groups for these genes----
tmp <- subset(dtl,
              SYMBOL %in% gene.keep)
tmp

p12 <- ggplot(data = tmp) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = Expression),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = range(dtl$Expression),
                       name = "Expression") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Sorted by Most Changed in LNCaP PEITC vs. Setd7 PEITC") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p12

tiff(filename = "tmp/all_anno_genes_by_dtrt.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p12)
graphics.off()

# Sort by differences in Setd7----
tmp <- subset(dt.diff,
              Group == "Setd7 PEITC - Setd7")
tmp <- tmp[order(tmp$`Diff of Logs`), ]
ordered.genes <- as.character(tmp$SYMBOL)

gene.keep <- ordered.genes[c(1:20,
                             (length(ordered.genes) - 19):length(ordered.genes))]

tmp <- subset(dt.diff,
              SYMBOL %in% gene.keep)

tmp$SYMBOL <- factor(tmp$SYMBOL,
                     levels = gene.keep)
tmp

# Plot----
p03 <- ggplot(data = tmp) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = `Diff of Logs`),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       limit = range(dt.diff$`Diff of Logs`),
                       name = "Diff of Logs") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Genes Most Changed in Setd7 PEITC vs. Setd7") +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p03

tiff(filename = "tmp/all_anno_genes_diffs_3.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p03)
graphics.off()

# Plot all 4 groups for these genes----
tmp <- subset(dtl,
              SYMBOL %in% gene.keep)
tmp

p13 <- ggplot(data = tmp) +
  geom_tile(aes(x =  Group,
                y = SYMBOL,
                fill = Expression),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = range(dtl$Expression),
                       name = "Expression") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene Name",
                   expand = c(0, 0)) +
  ggtitle("Genes Most Changed in Setd7 PEITC vs. Setd7") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p13

tiff(filename = "tmp/all_anno_genes_by_dsetd7.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p13)
graphics.off()

# sink()