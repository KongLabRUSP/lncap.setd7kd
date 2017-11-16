# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCap/SETD7 KD experiment             |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/14/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_lncap_setd7kd_analysis_v1.txt")

require(data.table)
require(ggplot2)

# Part I: Data----
# Move up one directory
wd <- getwd()
setwd("..")
DATA_HOME <- paste(getwd(),
                   "data/lncap.setd7kd",
                   sep = "/")
# Reset working directory
setwd(wd)
getwd()

lData <- dir(DATA_HOME)
lData

dt1 <- fread(paste(DATA_HOME,
                   "LNCaP vs LNCaP PEITC.csv",
                   sep = "/"),
             skip = 4,
             header = TRUE,
             colClasses = "character")[, c(1, 3, 6, 9)]
names(dt1)[1:2] <- c("LNCaP",
                     "LNCaP PEITC")
summary(dt1)
dt1 <- subset(dt1,
              `Other ID` != "")
dt1

dt2 <- fread(paste(DATA_HOME,
                   "Setd7 vs Setd7 PEITC.csv",
                   sep = "/"),
             skip = 5,
             header = TRUE,
             colClasses = "character")[, c(1, 3, 6, 9)]
names(dt2)[1:2] <- c("Setd7",
                     "Setd7 PEITC")

dt2 <- subset(dt2,
              `Other ID` != "")
dt2

# dtt <- merge(dt1,
#              dt2,
#              by = "Gene ID",
#              all = TRUE)
dtt <- merge(dt1,
             dt2,
             by = "Other ID")
dtt$LNCaP <- as.numeric(dtt$LNCaP)
dtt$`LNCaP PEITC` <- as.numeric(dtt$`LNCaP PEITC`)
dtt$Setd7 <- as.numeric(dtt$Setd7)
dtt$`Setd7 PEITC` <- as.numeric(dtt$`Setd7 PEITC`)
summary(dtt)
dtt

# Get mapping file: affimetrix ID to gene name----
dt1.map <- fread(paste(DATA_HOME,
                       "LNCaP vs LNCaP PEITC ID-gene name.csv",
                       sep = "/"),
                 header = TRUE,
                 colClasses = "character")
dt1.map <- unique(dt1.map)
dt1.map <- subset(dt1.map,
                  Symbol != "")

dt2.map <- fread(paste(DATA_HOME,
                       "Setd7 vs Setd7 PEITC ID-gene name.csv",
                       sep = "/"),
                 header = TRUE,
                 colClasses = "character")
dt2.map <- unique(dt2.map)
dt2.map <- subset(dt2.map,
                  Symbol != "")

dt3.map <- fread(paste(DATA_HOME,
                       "LNCaP vs Setd7 ID-gene name.csv",
                       sep = "/"),
                 header = TRUE,
                 colClasses = "character")
dt3.map <- unique(dt3.map)
dt3.map <- subset(dt3.map,
                  Symbol != "")

dt.map <- merge(dt1.map,
                dt2.map,
                by = names(dt1.map),
                all = TRUE)
dt.map <- merge(dt.map,
                dt3.map,
                by = names(dt.map),
                all = TRUE)
names(dt.map)[1] <- "Other ID"
dt.map

# Merge data and mappiing
dtt <- merge(dt.map, 
             dtt,
             by = "Other ID")
summary(dtt)

# NOTE: 2 genes are not mapped uniquely; rename one of the copies
dtt$Symbol[dtt$`Other ID` == "17114581"] <- "RAB31_2"
dtt$Symbol[dtt$`Other ID` == "16726880"] <- "NEAT1_2"
length(unique(dtt$Symbol))

# Heatmap: all 4 groups----
dtl <- melt.data.table(dtt,
                       id.vars = "Symbol",
                       measure.vars = c(6, 7, 9, 10),
                       variable.name = "Group",
                       value.name = "Readout")
dtl$Symbol <- factor(dtl$Symbol)
dtl$Group <- factor(dtl$Group,
                    levels = c("LNCaP",
                               "Setd7",
                               "LNCaP PEITC",
                               "Setd7 PEITC"))
dtl

# Plot all annotated genes found in all 4 samples----
p0 <- ggplot(data = dtl) +
  geom_tile(aes(x =  Group,
                y = Symbol,
                fill = Readout),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = range(dtl$Readout),
                       name = "Readout") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene Name",
                   expand = c(0, 0)) +
  ggtitle("All Annotated Genes Found Across All 4 Samples") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p0

tiff(filename = "tmp/all_anno_genes.tiff",
     height = 10,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p0)
graphics.off()

# Hitmap of differences----
dtt$`LNCaP - Setd7` <- dtt$LNCaP - dtt$Setd7
dtt$`LNCaP PEITC - Setd7 PEITC` <- dtt$`LNCaP PEITC` - dtt$`Setd7 PEITC`

dt.diff <- melt.data.table(dtt,
                           id.vars = "Symbol",
                           measure.vars = 12:13,
                           variable.name = "Group",
                           value.name = "Diff of Logs")
dt.diff$Symbol <- factor(dt.diff$Symbol)
dt.diff$Group <- factor(dt.diff$Group,
                        levels = c("LNCaP - Setd7",
                                   "LNCaP PEITC - Setd7 PEITC"))
dt.diff
length(unique(dt.diff$Symbol))
dt.diff[duplicated(dt.diff$Symbol), ]

# Sort by differences in controls
tmp <- subset(dt.diff,
              Group == "LNCaP - Setd7")
tmp <-tmp[order(tmp$`Diff of Logs`), ]
length(unique(tmp$Symbol))
ordered.genes <- as.character(tmp$Symbol)

dt.diff$Symbol <- factor(dt.diff$Symbol,
                         levels = ordered.genes)

# Plot all annotated genes found in all 4 samples----
p00 <- ggplot(data = dt.diff) +
  geom_tile(aes(x =  Group,
                y = Symbol,
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
p00

tiff(filename = "tmp/all_anno_genes_diffs.tiff",
     height = 10,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p00)
graphics.off()

# Sort by differences in treated----
tmp <- subset(dt.diff,
              Group == "LNCaP PEITC - Setd7 PEITC")
tmp <-tmp[order(tmp$`Diff of Logs`), ]
length(unique(tmp$Symbol))
ordered.genes <- as.character(tmp$Symbol)

dt.diff$Symbol <- factor(dt.diff$Symbol,
                         levels = ordered.genes)

# Plot all annotated genes found in all 4 samples----
p01 <- ggplot(data = dt.diff) +
  geom_tile(aes(x =  Group,
                y = Symbol,
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
     height = 10,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p01)
graphics.off()

# Keep only the genes that changed----
dtt1 <- subset(dtt,
               ((`LNCaP - Setd7` > 0) & 
                  (`LNCaP PEITC - Setd7 PEITC` < 0)) |
                 ((`LNCaP - Setd7` < 0) & 
                    (`LNCaP PEITC - Setd7 PEITC` > 0)))
# Differences----
dt.diff <- melt.data.table(dtt1,
                           id.vars = "Symbol",
                           measure.vars = 12:13,
                           variable.name = "Group",
                           value.name = "Diff of Logs")
dt.diff$Symbol <- factor(dt.diff$Symbol)
dt.diff$Group <- factor(dt.diff$Group,
                        levels = c("LNCaP - Setd7",
                                   "LNCaP PEITC - Setd7 PEITC"))
dt.diff
length(unique(dt.diff$Symbol))
dt.diff[duplicated(dt.diff$Symbol), ]

# Sort
tmp <- subset(dt.diff,
              Group == "LNCaP - Setd7")
tmp <- tmp[order(tmp$`Diff of Logs`), ]
length(unique(tmp$Symbol))
ordered.genes <- as.character(tmp$Symbol)

dt.diff$Symbol <- factor(dt.diff$Symbol,
                         levels = ordered.genes)

# Plot all annotated genes found in all 4 samples----
p02 <- ggplot(data = dt.diff) +
  geom_tile(aes(x =  Group,
                y = Symbol,
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

# # Genes orered by ratios---
# dtt$Ratio.x <- as.numeric(dtt$Ratio.x)
# range(dtt$Ratio.x, na.rm = TRUE)
# 
# dtt$Ratio.y <- as.numeric(dtt$Ratio.y)
# range(dtt$Ratio.y, na.rm = TRUE)
# 
# # # Select 20 genes at random to plot
# # keepID <- sample(levels(dtl$GeneID), 20)
# # tmp <- droplevels(dtl[GeneID %in% keepID, ])
# 
# # Top LNCaP genes (by ratio)----
# tmp <- subset(dtt,
#               !is.na(dtt$Ratio.x))
# top20lc <- tmp$`Gene ID`[order(tmp$Ratio.x,
#                                decreasing = TRUE)][1:20]
# tmp <- droplevels(dtl[`Gene ID` %in% top20lc, ])
# # tmp <- droplevels(dtl[GeneID %in% bot20lc, ])
# 
# p1 <- ggplot(data = tmp) +
#   geom_tile(aes(x =  Group,
#                 y = `Gene ID`,
#                 fill = Readout),
#             color = "black") +
#   scale_fill_gradient2(high = "red",
#                        limit = range(dtl$Readout),
#                        name = "Microarray Readout") +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete("Gene ID",
#                    expand = c(0, 0)) +
#   ggtitle("Top 20 Ratio LNCaP vs. LNCaP-PEITC") +
#   theme(axis.text.x = element_text(angle = 30,
#                                    hjust = 1),
#         # legend.position = "top",
#         plot.title = element_text(hjust = 0.5))
# p1
# 
# tiff(filename = "tmp/LNCaP_LNCaP-PEITC.tiff",
#      height = 5,
#      width = 5,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p1)
# graphics.off()
# 
# # Top Setd7 genes (by ratio)----
# tmp <- subset(dtt,
#               !is.na(dtt$Ratio.y))
# top20s7 <- tmp$`Gene ID`[order(tmp$Ratio.y,
#                                decreasing = TRUE)][1:20]
# tmp <- droplevels(dtl[GeneID %in% top20s7, ])
# # tmp <- droplevels(dtl[GeneID %in% bot20lc, ])
# 
# p2 <- ggplot(data = tmp) +
#   geom_tile(aes(x =  Group,
#                 y = GeneID,
#                 fill = Readout),
#             color = "black") +
#   scale_fill_gradient2(high = "red",
#                        limit = range(dtl$Readout),
#                        name = "Microarray Readout") +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete("Gene ID",
#                    expand = c(0, 0)) +
#   ggtitle("Top 20 Ratio Setd7 vs. Setd7-PEITC") +
#   theme(axis.text.x = element_text(angle = 30,
#                                    hjust = 1),
#         # legend.position = "top",
#         plot.title = element_text(hjust = 0.5))
# p2
# 
# tiff(filename = "tmp/Setd7_Setd7-PEITC.tiff",
#      height = 5,
#      width = 5,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p2)
# graphics.off()
# 
# # Entrez ID to gene name
# entraz_id <- list()
# for (i in 1:length(gene_names)) {
#   out <- try(unlist(mget(x = gene_names[i],
#                          envir = org.Mm.egALIAS2EG)),
#              silent = TRUE)
#   if (class(out)[1] != "try-error") {
#     entraz_id[i] <- out
#   } else {
#     entraz_id[i] <- NA
#   }
# }
# 
# entraz_id <- list()
# 
# gene_names <- c("BC137383")
# for (i in 1:length(gene_names)) {
#   out <- try(unlist(mget(x = gene_names[i],
#                          envir = org.Mm.egACCNUM2EG)),
#              silent = TRUE)
#   if (class(out)[1] != "try-error") {
#     entraz_id[i] <- out
#   } else {
#     entraz_id[i] <- NA
#   }
# }
# 
# entraz_id
# 
# org.Mm.eg.db
# org.Mm.eg_dbInfo()
# org.Mm.egACCNUM
# org.Mm.egACCNUM2EG