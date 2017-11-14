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

# PArt I: Data----
lData <- dir("data")
lData

dt1 <- fread(paste("data",
                   lData[1],
                   sep = "/"),
             skip = 4,
             header = TRUE)[, c(1, 3, 5, 6, 9)]
dt1$`Other ID` <- as.numeric(dt1$`Other ID`)
dt1 <- subset(dt1,
              !is.na(`Other ID`))
names(dt1)[1:2] <- c("LNCaP",
                     "LNCaP PEITC")
dt1

dt4 <- fread(paste("data",
                   lData[4],
                   sep = "/"),
             skip = 5,
             header = TRUE)[, c(1, 3, 5, 6, 9)]
dt4$`Other ID` <- as.numeric(dt4$`Other ID`)
dt4 <- subset(dt4,
              !is.na(`Other ID`))
names(dt4)[1:2] <- c("Setd7",
                     "Setd7 PEITC")
dt4

dtt <- merge(dt1,
             dt4,
             by = "Other ID")
names(dtt) [1] <- "GeneID"
dtt

# Heatmap: all 4 groups----
dtl <- melt.data.table(dtt,
                       id.vars = "GeneID",
                       measure.vars = c(2, 3, 6, 7),
                       variable.name = "Group",
                       value.name = "Readout")
dtl$GeneID <- factor(dtl$GeneID)
dtl$Group <- factor(dtl$Group,
                    levels = unique(dtl$Group))
dtl$Readout <- as.numeric(dtl$Readout)
range(dtl$Readout)
dtl

# Genes orered by ratios---
dtt$Ratio.x <- as.numeric(dtt$Ratio.x)
range(dtt$Ratio.x)

dtt$Ratio.y <- as.numeric(dtt$Ratio.y)
range(dtt$Ratio.y)

top20lc <- dtt$GeneID[order(dtt$Ratio.x,
                            decreasing = TRUE)][1:20]
# bot20lc <- dtt$GeneID[order(dtt$Ratio.x,
#                             decreasing = FALSE)][1:20]
top20s7 <- dtt$GeneID[order(dtt$Ratio.y,
                            decreasing = TRUE)][1:20]
# bot20s7<- dtt$GeneID[order(dtt$Ratio.y,
#                            decreasing = FALSE)][1:20]

# # Select 20 genes at random to plot
# keepID <- sample(levels(dtl$GeneID), 20)
# tmp <- droplevels(dtl[GeneID %in% keepID, ])

# Top LNCaP genes (by ratio)----
tmp <- droplevels(dtl[GeneID %in% top20lc, ])
# tmp <- droplevels(dtl[GeneID %in% bot20lc, ])

p1 <- ggplot(data = tmp) +
  geom_tile(aes(x =  Group,
                y = GeneID,
                fill = Readout),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = range(dtl$Readout),
                       name = "Microarray Readout") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene ID",
                   expand = c(0, 0)) +
  ggtitle("Top 20 Ratio LNCaP vs. LNCaP-PEITC") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "tmp/LNCaP_LNCaP-PEITC.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Top Setd7 genes (by ratio)----
tmp <- droplevels(dtl[GeneID %in% top20s7, ])
# tmp <- droplevels(dtl[GeneID %in% bot20lc, ])

p2 <- ggplot(data = tmp) +
  geom_tile(aes(x =  Group,
                y = GeneID,
                fill = Readout),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = range(dtl$Readout),
                       name = "Microarray Readout") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene ID",
                   expand = c(0, 0)) +
  ggtitle("Top 20 Ratio Setd7 vs. Setd7-PEITC") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p2

tiff(filename = "tmp/Setd7_Setd7-PEITC.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()