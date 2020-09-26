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

library(pROC) 
library(ggplot2)
load(  file = 
         file.path('Rdata/','TCGA_RPKM_8circRNA_ROC_patch.Rdata'))

meta$T_stage2=ifelse(meta$T_stage1 == 'high',1,0)
#plot ROC curve by roc function from pROC packages
roc.list <- roc(BCR ~ Rscore + PSA +Gleason+T_stage2, data = meta)
g.list <- ggroc(roc.list)
g.list
p <- g.list+theme_classic()+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed")+
  #add AUC values 
  annotate("text",x = .25, y = .32,label = paste("Riskscore = ",round(roc.list[[1]]$auc,3)),color = "red")+
  annotate("text",x = .25, y = .24,label = paste("PSA = ",round(roc.list[[2]]$auc,3)),color = "green")+
  annotate("text",x = .25, y = .16,label = paste("Gleason  = ",round(roc.list[[3]]$auc,3)),color = "blue")+
  annotate("text",x = .25, y = .08,label = paste("Tstage  = ",round(roc.list[[4]]$auc,3)),color = "purple")

pdf(file = "./Figure/step02LASSO/ROCcurve.pdf")
p
dev.off()           
