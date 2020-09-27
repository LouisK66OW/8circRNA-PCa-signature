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
  library(survival)
  library(survminer)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
}

#Load File
load('./Rdata/RPKM_8circRNA.Rdata')
rt=read.table("./Rdata/step03_surInput.txt",header=T,sep="\t",check.names=F,row.names=1) 
meta_fen <- meta[,c(1,18)]
newmeta <- rt%>%
  left_join(meta_fen,by='ID')

#Plot the boxplot of immune score between two Risk score group.
#StromalScore
p1 <- ggplot(newmeta,aes(x=Rscore_group,y=StromalScore,fill=Rscore_group))+
  stat_boxplot(geom ="errorbar")+
  geom_boxplot(outlier.shape = 21,colour = "black",outlier.alpha=0.2,varwidth = TRUE)+
  scale_fill_manual(values = c("#CC0000", "#2f5688"))+
  stat_compare_means( label = "p.format",label.x=1.85,method="wilcox.test")+
  theme_bw()
# ImmuneScore
p2 <- ggplot(newmeta,aes(x=Rscore_group,y=ImmuneScore,fill=Rscore_group))+
  stat_boxplot(geom ="errorbar")+
  geom_boxplot(outlier.shape = 21,colour = "black",outlier.alpha=0.2,varwidth = TRUE)+
  scale_fill_manual(values = c("#CC0000", "#2f5688"))+
  stat_compare_means( label = "p.format",label.x=1.85,method="wilcox.test")+
  theme_bw()
# ESTIMATEScore
p3 <- ggplot(newmeta,aes(x=Rscore_group,y=ESTIMATEScore,fill=Rscore_group))+
  stat_boxplot(geom ="errorbar")+
  geom_boxplot(outlier.shape = 21,colour = "black",outlier.alpha=0.2,varwidth = TRUE)+
  scale_fill_manual(values = c("#CC0000", "#2f5688"))+
  stat_compare_means( label = "p.format",label.x=1.85,method="wilcox.test")+
  theme_bw()

#Print boxplot as pdf
pdf(file="./Figures/step05_boxplot_ESTIMATEScore.pdf")
p1;p2;p3
dev.off()
