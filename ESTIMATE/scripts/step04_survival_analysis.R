
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
  library(survival)
  library(survminer)
  library(tidyverse)
  library(ggplot2)
}

#Read File
rt=read.table("./Rdata/step03_surInput.txt",header=T,sep="\t",check.names=F,row.names=1) 
colnames(rt)

outTab=data.frame()
for(score in colnames(rt[,4:ncol(rt)])){
  a=rt[,score]<= median(rt[,score])
  diff=survdiff(Surv(TBCR, BCR) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  fit=survfit(Surv(TBCR, BCR) ~ a, data = rt)
  
  outTab=rbind(outTab,cbind(score=score,pval=pValue) )
  #Draw survival curve
  if(pValue<0.001){
    pValue="<0.001"
  }else{
    pValue=paste0("=",sprintf("%0.3f",pValue))
  }
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=F,
                     pval=paste0("p",pValue),
                     pval.size=6,
                     risk.table=T,
                     legend.labs=c("high","low"),
                     legend.title=gsub("Score"," Score",score),
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("#CC0000", "#2f5688"),
                     risk.table.height=.25)
  #print the survival curve as pdf
  pdf(file=paste0("./Figures/step04_sur.",score,".pdf"), width = 6, height = 5.5,onefile = FALSE)
  print(surPlot)
  dev.off()
}

