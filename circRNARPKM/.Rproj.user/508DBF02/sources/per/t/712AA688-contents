### ---------------
###
### Create: Shuo Wang
### Date: 2020-04-12
### Email: mywanuo@163.com
### Southern Medical University,Guang Zhou / Major in Urology Surgery 
### Update Log: 2020-09-26 11:17:47  Modified versions
###
### ---------------

rm(list=ls())
options(stringsAsFactors = F)

library(pheatmap)
library(gplots)

load( file = './Rdata/PRAD-circRNA-sample.Rdata')


#Data Processing,normalization by scale
exprSet_scale <- as.matrix((scale(exprSet))) 
exprSet_scale[exprSet_scale > 2] = 2
exprSet_scale[exprSet_scale< -2] = -2

meta <- phe
meta$T_stage2=ifelse(meta$T_stage1=="high",'T2','T1')
meta$BCR_stage = ifelse(meta$BCR=="1",'YES','NO')

# construct annotation_col of pheatmapfunction
annotation_col = data.frame(
  BCR = factor(meta$BCR_stage),
  `T stage` = meta$T_stage2, 
  `Gleason score` = meta$Gleason)
rownames(annotation_col) = colnames(exprSet)
annotation_col <- annotation_col[order(annotation_col$BCR),, drop = FALSE]
exprSet <- exprSet[,rownames(annotation_col)]
colnames(annotation_col) = c("BCR","T stage","Gleason score")

pheatmap(exprSet_scale, annotation_col = annotation_col)
pheatmap(exprSet_scale,annotation_col=annotation_col,labels_col = '',
         cluster_cols = FALSE, cellwidth = 3, cellheight = 15
         ,filename = "./Figure/step01KMcurves/pheatmap28_log_scale.pdf"
         )
