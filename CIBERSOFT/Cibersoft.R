##
### ---------------
###
### Create: Shuo Wang
### Date: 2020-04-25
### Email: mywanuo@163.com
### Southern Medical University,Guang Zhou / Major in Urology Surgery 
### Update Log: 2020-09-28 17:15:19 CST  Modified versions
###
### ---------------

#visulization------------------------
rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

rt=read.table("./CIBERSORT.Output_net100.txt",header=T,sep="\t",check.names=F)
#the data computed and downloaded from http://cibersort.stanford.edu/
colnames(rt)

cibersort_raw <-rt%>%
  rename("Patients" = "Input Sample") %>%
  select(-c("P-value","Pearson Correlation","RMSE"))
# The first column name "Input Sample" is changed to "Patiens",delect "P-value","Pearson Correlation","RMSE" column 

#The histogram shows the percentage of cells
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
cibersort_barplot <- cibersort_raw %>%
  gather(key = Cell_type,value = Proportion,2:23)
#Use the RColorBrewer package to configure the color scheme, use the key-value correspondence
# in the gather function to reconstruct the correspondence between the cell name and the ratio

#Load meta clinical information
load(file='TCGA_RPKM_8circRNA.Rdata')
metanew <- meta[,c('ID','Rscore_group')]
colnames(metanew) <- c('Patients','R_score')
cibersort_barplot <- cibersort_barplot %>%
  left_join(metanew, c("Patients"))#add the Risk score group to  the cibersort_barplot

#plot bar plot of cells
ggplot(cibersort_barplot,aes(Patients,Proportion,fill = Cell_type)) + 
  geom_bar(position = "stack",stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))
# ggsave('barplot.pdf', width = 20, height = 20, units = "cm")

#Plot boxplot of cells with Risk score group
p <- ggplot(cibersort_barplot,aes(x=Cell_type,y=Proportion,fill=R_score))+
  stat_boxplot(geom ="errorbar")+
  geom_boxplot(outlier.shape = 21,colour = "black",outlier.alpha=0.2,varwidth = TRUE)+
  scale_fill_manual(values = c("#CC0000", "#2f5688"))+
  scale_y_continuous(limits = c(0,0.6),breaks=seq(0,0.6,0.1))+
  stat_compare_means(aes(group=R_score), label = "p.format",label.y=0.575,method="wilcox.test")+
  theme_bw()+
  theme(
    legend.position = "top",
    legend.background=element_blank(),
    legend.key = element_blank(),
    legend.margin=margin(0,0,0,0,"mm"),
    axis.text.x=element_text(size=rel(1.1),face="bold"),
    axis.line.x = element_line(size = 0.5, colour = "black"),
    axis.line.y = element_line(size = 0.5, colour = "black"),
    legend.text=element_text(size=rel(1.1)),
    legend.title=element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )+
  guides(color=FALSE)+
  coord_flip()
p
# ggsave("boxplot_100.pdf", p, width = 17, height = 17, units = "cm")
  