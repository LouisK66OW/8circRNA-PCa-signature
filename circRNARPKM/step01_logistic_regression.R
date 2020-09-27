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

library(survival)
library(survminer)

#=============================================================================
#survival analysis
#=============================================================================
load( file = "./Rdata/PRAD-circRNA.Rdata")
mySurv=with(phe,Surv(TBCR, BCR))
log_rank_p <- apply(expr, 1 , function(gene){
  phe$expr_group=ifelse(as.numeric(gene)>as.numeric(median(gene)),'high','low')  
  data.survdiff=survdiff(mySurv~expr_group,data=phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
log_rank_p=sort(log_rank_p)
boxplot(log_rank_p)
table(log_rank_p<0.01)
log_rank_p[log_rank_p<0.01]

table(log_rank_p<0.05)
log_rank_p[log_rank_p<0.05]

#choosing circRNAs with P-value < 0.05
exprSet <- expr[names(log_rank_p[log_rank_p<0.05]),]
save(expr,exprSet,phe,file = "./Rdata/PRAD-circRNA-sample.Rdata")
