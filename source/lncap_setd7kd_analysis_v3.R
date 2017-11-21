# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCap/SETD7 KD experiment             |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/20/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_lncap_setd7kd_analysis_v3.txt")

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

# Heatmap: all 4 groups----
dtl <- melt.data.table(dt1,
                       id.vars = "SYMBOL",
                       measure.vars = 5:10,
                       variable.name = "Group",
                       value.name = "Expression")
dtl$SYMBOL <- factor(dtl$SYMBOL)
dtl$Group <- factor(dtl$Group,
                    levels = c("LNCaP",
                               "Setd7",
                               "LNCaP_PEITC",
                               "Setd7_PEITC",
                               "LNCaP_TNF",
                               "Setd7_TNF"))
dtl
summary(dtl)

# Plot all annotated genes found in all 4 samples----
set.seed(2017)
tmp <- subset(dtl,
              SYMBOL %in% sample(unique(SYMBOL), 40))
tmp

p0 <- ggplot(data = tmp) +
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
  ggtitle("Random 40 Annotated Genes Found Across All 4 Samples") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p0

tiff(filename = "tmp/all_anno_genes.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p0)
graphics.off()

# Hitmap of differences----
dt1$`LNCaP - Setd7` <- dt1$LNCaP - dt1$Setd7
dt1$`LNCaP PEITC - Setd7 PEITC` <- dt1$LNCaP_PEITC - dt1$Setd7_PEITC
dt1$`LNCaP TNF - Setd7 TNF` <- dt1$LNCaP_TNF - dt1$Setd7_TNF
dt1

dt.diff <- melt.data.table(dt1,
                           id.vars = "SYMBOL",
                           measure.vars = 11:13,
                           variable.name = "Group",
                           value.name = "Diff of Logs")
dt.diff$SYMBOL <- factor(dt.diff$SYMBOL)
dt.diff$Group <- factor(dt.diff$Group,
                        levels = c("LNCaP - Setd7",
                                   "LNCaP PEITC - Setd7 PEITC",
                                   'LNCaP TNF - Setd7 TNF'))
dt.diff
length(unique(dt.diff$SYMBOL))


# Sort by differences in controls----
tmp <- subset(dt.diff,
              Group == "LNCaP - Setd7")
tmp <- tmp[order(tmp$`Diff of Logs`), ]
ordered.genes <- as.character(tmp$SYMBOL)

tmp <- subset(dt.diff,
              SYMBOL %in% ordered.genes[c(1:20,
                                          (length(ordered.genes) - 19):length(ordered.genes))])

tmp$SYMBOL <- factor(tmp$SYMBOL,
                         levels = ordered.genes[c(1:20,
                                                  (length(ordered.genes) - 19):length(ordered.genes))])
tmp

# # Sort by differences in treated----
# tmp <- subset(dt.diff,
#               Group == "LNCaP PEITC - Setd7 PEITC")
# tmp <-tmp[order(tmp$`Diff of Logs`), ]
# length(unique(tmp$SYMBOL))
# ordered.genes <- as.character(tmp$SYMBOL)
# 
# dt.diff$SYMBOL <- factor(dt.diff$SYMBOL,
#                          levels = ordered.genes)

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
  scale_y_discrete("Gene ID",
                   expand = c(0, 0)) +
  ggtitle("Celline Differences") +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1),
        # legend.position = "top",
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

# Sort by differences in treated----
tmp <- subset(dt.diff,
              Group == "LNCaP PEITC - Setd7 PEITC")
tmp <- tmp[order(tmp$`Diff of Logs`), ]
ordered.genes <- as.character(tmp$SYMBOL)

tmp <- subset(dt.diff,
              SYMBOL %in% ordered.genes[c(1:20,
                                          (length(ordered.genes) - 19):length(ordered.genes))])

tmp$SYMBOL <- factor(tmp$SYMBOL,
                     levels = ordered.genes[c(1:20,
                                              (length(ordered.genes) - 19):length(ordered.genes))])
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
  scale_y_discrete("Gene ID",
                   expand = c(0, 0)) +
  ggtitle("Celline Differences") +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1),
        # legend.position = "top",
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

# sink()