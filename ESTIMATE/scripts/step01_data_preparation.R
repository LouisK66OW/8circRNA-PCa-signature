
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

#Read File
library(tidyverse)
library(openxlsx)
load('./Rdata/GSE113120_mRNA_count.Rdata')
exprSet <- rownames_to_column(exprSet,var = "geneID")

load('./Rdata/sysdata.rda')
expran <- exprSet%>%
  left_join(biotype,c('geneID'='geneSymbol'))
expran <- expran[!is.na(expran$group),]
expran <- expran[!duplicated(expran$geneID),]
expran <- expran[,c(1,146:150,2:145)]
expran[1:4,1:5]
as.data.frame(table(expran$biotype))
as.data.frame(table(expran$group))
#选出mRNA 和lncRNA
expr_mRNA <- expran[expran$group=='protein_coding',]
expr_mRNA <- expr_mRNA[,c(1,7:150)]
save(expr_mRNA,file='./Rdata/step01_expr_mRNA_annotation.Rdata')
