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
# source("http://bioconductor.org/biocLite.R")
# biocLite("oligo")
# biocLite("pd.hugene.2.0.st")
# biocLite("hugene20sttranscriptcluster.db")
# biocLite("edgeR")
# biocLite("affycoretools")
# biocLite("VariantAnnotation")

require(data.table)
require(ggplot2)
require(oligo)
require(hugene20sttranscriptcluster.db)
require(affycoretools)

# Part I: Data----
lData <- dir("data/Human Cel files")
lData

l1 <- read.celfiles(filenames = paste("data/Human Cel files",
                                 lData,
                                 sep = "/"))
l1
class(l1)

# Annotation
l1@annotation
?pd.hugene.2.0.st

# Sample names
sampleNames(l1) <- c("LNCaP_PEITC",
                     "LNCaP_TNF",
                     "LNCaP",
                     "Setd7_PEITC",
                     "Setd7_TNF",
                     "Setd7")
l1

# Normilazed/bg corected genes----
# genePS <- rma(l1,
#               target = "probeset")
# NOTE: no mapping can be done for target probeset
genePS <- rma(l1,
              target = "core")

# genePS <- rma(l1,
#               target = "core",
#               background = FALSE,
#               normalize = FALSE)
genePS
dt1 <- exprs(genePS)
dt1 <- data.table(PROBEID = rownames(dt1),
                  dt1)
dt1

# Annotate----
# dt2 <- annotateEset(genePS,
#                     hugene20sttranscriptcluster.db)
dt2 <- annotateEset(object = genePS,
                    x = l1@annotation)

anno <- data.table(dt2@featureData@data)
anno$PROBEID <- as.character(anno$PROBEID)
anno
summary(anno$SYMBOL)
length(unique(anno$SYMBOL))
# hugene20sttranscriptcluster.db: 28,019 genes annotated
# pd.hugene.2.0.st: 34,662

# Merge expressions and annotation----
dt1 <- merge(anno, 
             dt1,
             by = "PROBEID")
dt1

# Keep only the annotated genes----
dt1$SYMBOL <- as.character(dt1$SYMBOL)
length(unique(dt1$SYMBOL))

dt1 <- subset(dt1,
              !is.na(SYMBOL))
dt1
summary(dt1)


tmp <- subset(dt1,
              SYMBOL %in% c("RARRES3",
                            "GBP2"))
tmp
tmp <- subset(dt1,
              ID %in% c("NM_173505"))
tmp

# Raw expressions----
head(l1@assayData$exprs)
head(l1@featureData@data)
# Same as:
dt0 <- oligo::exprs(l1)
l1@phenoData@data
head(dt0)
summary(dt0)

head(dt1)
summary(dt1)

# Keep only the annotated portion----
dt0 <- data.table(exprs(l1))

# Raw data distribution----
set.seed(2017)
tmp <- melt.data.table(dt0[sample(x = 1:nrow(dt1),
                                  size = 10000)])
# Remove TNF----
tmp <- droplevels(subset(tmp,
                         !(variable %in% c("LNCaP_TNF",
                                           "Setd7_TNF"))))

p1 <- ggplot(data = tmp) +
  geom_boxplot(aes(x = variable,
                   y = log2(value),
                   fill = variable)) +
  scale_x_discrete("Group") + 
  scale_y_continuous("log2(Expression)") + 
  ggtitle("Random Sample of 10,000 Raw Expressions") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
print(p1)

tiff(filename = "tmp/lncap_setd7_raw_expr.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Normilazed and annotated data distribution----
tmp <- melt.data.table(dt1,
                       measure.vars = 5:10)

# Remove TNF----
tmp <- droplevels(subset(tmp,
                         !(variable %in% c("LNCaP_TNF",
                                           "Setd7_TNF"))))

p2 <- ggplot(data = tmp) +
  geom_boxplot(aes(x = variable,
                   y = log2(value),
                   fill = variable)) +
  scale_x_discrete("Group") + 
  scale_y_continuous("log2(Expression)") + 
  ggtitle("All Normilazed Annotated Genes") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
print(p2)

tiff(filename = "tmp/lncap_setd7_norm_anno_expr.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Pseudo-images of raw microarray chip----
# par(mfrow=c(2,3),
#     mai = c(0.1, 0.1, 0.1, 0.1))
l1@manufacturer
l1@intensityFile
head(l1@assayData$exprs) # expressions matrix
l1@phenoData@varMetadata
l1@phenoData@data
l1@phenoData@dimLabels
l1@featureData@data
l1@experimentData
l1@annotation
l1@protocolData
l1@.__classVersion__

# How is this plot constructed (tech)?
oligo::image(l1,
             which = 1)

for (i in 1:6) {
  tiff(filename = paste("tmp/lncap_setd7_raw_expr_",
                        sampleNames(l1)[i],
                        ".tiff",
                        sep = ""),
       height = 5,
       width = 5,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  
  oligo::image(l1,
               which = i,
               transfo = rank)
  
  graphics.off()
}

# MA Plot----
for (i in c(1, 3:4, 6)) {
  tiff(filename = paste("tmp/lncap_setd7_expr_maplot_",
                        i,
                        ".tiff",
                        sep = ""),
       height = 5,
       width = 5,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  MAplot(l1, 
         which = i)
  graphics.off()
}

# Save as CSV----
write.csv(dt1,
          file = "tmp/lncap_setd7_normilized_annotated.csv",
          row.names = FALSE)

# Save R data----
save(dt1,
     file = "data/lncap_setd7_normilized_annotated.RData")

# sink()