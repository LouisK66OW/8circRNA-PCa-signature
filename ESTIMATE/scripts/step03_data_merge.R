
### ---------------
###
### Create: Shuo Wang
### Date: 2020-04-25
### Email: mywanuo@163.com
### Southern Medical University,Guang Zhou / Major in Urology Surgery 
### Update Log: "2020-09-27 15:18:26 CST"  Modified versions
###
### ---------------

#Merge the BCR and immune score data
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
load('./Rdata/RPKM_8circRNA.Rdata')

data=read.table('./Rdata/step02_scores.txt',sep="\t",header = T,check.names = F)
colnames(meta)
colnames(data)
metaBCR <- meta[,1:3]
colnames(metaBCR)
newmeta <- metaBCR %>%
  left_join(data,c('ID'))
table(is.na(newmeta[,2:6]))
write.table(newmeta,file="./Rdata/step03_surInput.txt",sep="\t",quote=F,col.names=T)
