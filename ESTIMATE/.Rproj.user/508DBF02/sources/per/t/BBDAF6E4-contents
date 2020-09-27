##脚本前的注释
### ---------------
###
### Create: Shuo Wang
### Date: 2020-04-25
### Email: mywanuo@163.com
### Southern Medical University,Guang Zhou / Major in Urology Surgery 
### Update Log: "2020-09-27 15:18:26 CST"  Modified versions
###
### ---------------

rm(list = ls())
options(stringsAsFactors = F)
if(T){
  library(estimate)
  library(tidyverse)
  library(limma)
  library(estimate)
  }

#Read File
load("./Rdata/step02_expr_mRNA.Rdata")

max(data)
data=log2(data+1)
data=data[rowMeans(data)>0,]
dim(data)
data=normalizeBetweenArrays(data)
data=rbind(ID=colnames(data),data)
#Output the organized matrix file
write.table(data,file="./Rdata/step02_normalize.txt",sep="\t",quote=F,col.names=F)

#Running estimate package
filterCommonGenes(input.f="./Rdata/step02_normalize.txt", 
                  output.f="./Rdata/step02_commonGenes.gct", 
                  id="GeneSymbol")

#Calculating immune score
estimateScore(input.ds = "./Rdata/step02_commonGenes.gct",
              output.ds="./Rdata/step02_estimateScore.gct", 
              platform="illumina")

#Output the score of each samples
scores=read.table("./Rdata/step02_estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores=scores[,1:3]
out=rbind(ID=colnames(scores),scores)
write.table(out,file="./Rdata/step02_scores.txt",sep="\t",quote=F,col.names=F)