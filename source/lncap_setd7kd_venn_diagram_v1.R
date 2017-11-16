# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCap/SETD7 KD experiment             |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/15/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_lncap_setd7kd_venn_diagram_v1.txt")

require(data.table)
require(ggplot2)
require(VennDiagram)
require(gridExtra)

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

# Get mapping file: Affimetrix ID to gene name----
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

# Map LNCaP----
# Merge data and mappiing
dt1 <- merge(dt.map, 
             dt1,
             by = "Other ID")
gene1 <- unique(dt1$`Other ID`)

dt2 <- merge(dt.map, 
             dt2,
             by = "Other ID")
gene2 <- unique(dt2$`Other ID`)

# CHECKPOINT----
dtt <- merge(dt1, 
             dt2, 
             by = "Other ID")
nrow(dtt) == sum(gene1 %in% gene2)

# Draw Venn diagram----
p1 <- venn.diagram(x = list(LNCaP = gene1,
                            Setd7 = gene2),
                   filename = NULL,
                   fill = c("light blue", "grey"),
                   alpha = c(0.5, 0.5),
                   compression = "lzw+p",
                   main = "Number of Genes Found")
grid.arrange(gTree(children = p1))

png(filename = "tmp/venn_lncap_setd7.png",
    height = 5,
    width = 5,
    units = 'in',
    res = 300)
grid.arrange(gTree(children = p1))
graphics.off()