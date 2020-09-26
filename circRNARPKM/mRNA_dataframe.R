rm(list=ls())
options(stringsAsFactors = F)

library(openxlsx)
library(tidyverse)
load(file = './data/ICGC_exp.Rdata')
#排序
exprm <- expdata%>%
  arrange(submitted_sample_id,gene_id)

#行列名
colnam <- exprm[!duplicated(exprm$submitted_sample_id),]
colnam <- as.character(colnam$submitted_sample_id)
colnam
rownam <- exprm[exprm$submitted_sample_id=='CPCG0534-F1',]
rownam <- as.character(rownam$gene_id)
rownam

#构建矩阵
exprSet <- data.frame(geneID=rownam)
#exprm$normalized_read_count
expra <- lapply(colnam,function(x){
  e <- exprm[exprm$submitted_sample_id==x,]
  exprSet$x <- e$normalized_read_count
})
expr <- as.data.frame(expra)
colnames(expr) <- colnam
expr <- round(expr,digits =8)
expr[1:4,1:4]
expr <- cbind(rownam,expr)

expr2 <- expr[!duplicated(expr$rownam),]

exprSet <- expr2
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
exprSet[1:4,1:4]
save(exprSet,file='./data/GSE113120_mRNA_expr.Rdata')
load(file='./data/GSE113120_mRNA_expr.Rdata')

write.csv(exprSet,file = './data/GSE113120_mRNAexpr.csv',row.names = T)

expr <- rbind(genename =colnames(exprSet),exprSet)
expr[1:4,1:4]
expr[1,] <- gsub('-F1','',expr[1,])

exprSet <- expr[-1,]
exprSet[1:4,1:4]
colnames(exprSet) <- expr[1,]
