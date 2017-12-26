# |----------------------------------------------------------------------------------|
# | Project: Analysis of microarray data from a LNCap/SETD7 KD experiment            |
# | Script:                                                                          |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 12/19/2017                                                              |
# |----------------------------------------------------------------------------------|
require(data.table)

# Part I: Data----
# Preprocessed data was created by 'WT_KDKD_data_preprocessing_v1.R'
load("tmp/lncap_setd7_normilized_annotated.RData")
dt1

dt1$wt.kd <- dt1$LNCaP - dt1$Setd7
dt1$kdpeitc.kd <- dt1$Setd7_PEITC - dt1$Setd7
dt1$wt.kd.vs.kdpeitc.kd <- dt1$wt.kd - dt1$kdpeitc.kd

# Chao's list of genes (12/19/2017)----
dt2 <- fread("data/pick up genes from the heatmap.csv")
dt2

wt.kd.kdpeitc.kd.up.dn <- dt1[dt1$SYMBOL %in% dt2$`WT vs KD-KD PEITC vs KD Up vs Down`]
wt.kd.kdpeitc.kd.up.dn
write.csv(wt.kd.kdpeitc.kd.up.dn,
          file = "tmp/wt.kd.kdpeitc.kd.up.dn.csv")

wt.kd.kdpeitc.kd.dn.up <- dt1[dt1$SYMBOL %in% dt2$`WT vs KD-KD PEITC vs KD Down vs Up`]
wt.kd.kdpeitc.kd.dn.up
write.csv(wt.kd.kdpeitc.kd.dn.up,
          file = "tmp/wt.kd.kdpeitc.kd.dn.up.csv")

wtpietc.wt.kdpeitc.kd.up.dn <- dt1[dt1$SYMBOL %in% dt2$`WT PEITC vs WT-KD PEITC vs KD Up vs Down`]
wtpietc.wt.kdpeitc.kd.up.dn
write.csv(wtpietc.wt.kdpeitc.kd.up.dn,
          file = "tmp/wtpietc.wt.kdpeitc.kd.up.dn.csv")

wtpietc.wt.kdpeitc.kd.dn.up <- dt1[dt1$SYMBOL %in% dt2$`WT PEITC vs WT-KD PEITC vs KD Down vs Up`]
wtpietc.wt.kdpeitc.kd.dn.up
write.csv(wtpietc.wt.kdpeitc.kd.dn.up,
          file = "tmp/wtpietc.wt.kdpeitc.kd.dn.up.csv")
