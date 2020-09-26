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

#=============================================================
#Data processing
#=============================================================
if(T){
  library(ggplot2)
  library(glmnet)
  library(ROCR)
  library(caret)
  library(ggsci)
  library(ggpubr) 
  library(survival)
  library(survminer)
  library(cowplot)
  library(ggplotify)
  library(pheatmap)
  library(Hmisc)
}
load( file ="./Rdata/PRAD-circRNA-sample.Rdata")

## Selecting ircRNAs to construct a coxph model
e=t(exprSet[c('circ_30029','circ_117300','circ_176436','circ_112897', 'circ_17720',
              'circ_178252','circ_115617','circ_14736'),])

dat=cbind(phe,e)
colnames(dat) 
s=Surv(TBCR, BCR) ~ circ_30029+circ_117300+circ_176436+circ_112897+circ_17720+circ_178252+circ_115617+circ_14736
model <- coxph(s, data = dat )
summary(model,data=dat)
options(scipen=1)
ggforest(model, data =dat, 
         main = "Hazard ratio", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)+
  theme_classic2()

new_dat=dat
fp <- predict(model,new_dat,type="risk");boxplot(fp)
fp <- predict(model,new_dat,type="expected") ;boxplot(fp)
plot(fp,phe$TBCR)
fp <- predict(model,new_dat) ;boxplot(fp)
basehaz(model) 

options(scipen=200)
with(new_dat,rcorr.cens(fp,Surv(TBCR, BCR)  ))

fp_dat=data.frame(s=1:length(fp),v=as.numeric(sort(fp )))
sur_dat=data.frame(s=1:length(fp),
                   t=phe[names(sort(fp )),'TBCR'] ,
                   e=phe[names(sort(fp )),'BCR']  ) 
fp <- predict(model,dat) ;boxplot(fp)

fp_dat=data.frame(patientid=1:length(fp),fp=as.numeric(sort(fp)))
fp_dat$riskgroup = ifelse(fp_dat$fp>=median(fp_dat$fp),'high','low')
sur_dat=data.frame(patientid=1:length(fp),
                   time=phe[names(sort(fp)),'TBCR'] ,
                   event=phe[names(sort(fp )),'BCR']  ) 
sur_dat$event=ifelse(sur_dat$event=="1",'Biochemical Recurrence','Biochemical Recurrence Free')
exp_dat=dat[names(sort(fp)),(ncol(dat)-7):ncol(dat)]

###Draw a scatter plot of survival status----
if(T){
  p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
    scale_colour_manual(values = c("#2f5688","#CC0000"))+
    theme_bw()+labs(x="Patient ID(increasing risk score)",y="Risk score")+
    geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)+
    theme(legend.position="top")
  p1
}

###Draw high and low risk group display point diagram----
if(T){
  p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
    scale_colour_manual(values = c("#2f5688","#CC0000"))+
    labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
    geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)+
    theme(legend.position = "top")
  p2
}
#Plot heat map of 8 circRNAs
if(T){
  mycolors <- colorRampPalette(c("#2f5688", "white","#CC0000" ), bias = 1.2)(50)
  tmp=t(scale(exp_dat))
  a = 1
  tmp[tmp > a] = a
  tmp[tmp < -a] = -a
  #p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = T)
  colbreaks <- c(#将列分开
    72,
    72
  )
  
  p3=pheatmap(tmp,
              col= mycolors,
              show_colnames = F,cluster_cols = F,cluster_rows =F,
              gaps_col=colbreaks)
}
#Put pictures together
plots = list(p1,p2,as.ggplot(as.grob(p3)))
library(gridExtra)
lay1 = rbind(c(rep(1,7),NA),
             c(rep(2,7)),
             c(rep(3,7))) #layout
#> Warning in rbind(c(rep(1, 7), NA), c(rep(2, 7)), c(rep(3, 7))): number of
#> columns of result is not a multiple of vector length (arg 2)
pdf(file = "./Figure/step02LASSO/Risk_factor_correlation_diagram.pdf")
grid.arrange(grobs = plots, 
             layout_matrix = lay1,
             heigths = c(1, 2,3))
dev.off()
