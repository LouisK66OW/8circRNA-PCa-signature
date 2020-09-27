
### ---------------
###
### Create: Shuo Wang
### Date: 2020-04-25
### Email: mywanuo@163.com
### Southern Medical University,Guang Zhou / Major in Urology Surgery 
### Update Log: "2020-09-27 15:18:26 CST"  Modified versions
###
### ---------------

#Differential analysis---DEseq2

rm(list=ls())
options(stringsAsFactors = F)

if(T){
  library(DESeq2)
  library(tidyverse)  
  library(ggpubr)
  library(tidyverse)
}

load( file = './Rdata/step06_mRNA_count.Rdata')
group_list=meta[,18]
exprSet=na.omit(expr)

### ---------------
###
###run DESeq2 
###
### ---------------

if(T){
  colData <- data.frame(row.names=colnames(exprSet), 
                         group_list=group_list) 
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  tmp_f=c('./Rdata/mRNA-DESeq2-dds.Rdata')
}
  if(!file.exists(tmp_f)){
    dds <- DESeq(dds)
    save(dds,file = tmp_f)
  }
  load(file = tmp_f)
  res <- results(dds, 
                 contrast=c("group_list","high","low"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG =as.data.frame(resOrdered)
  DESeq2_DEG = na.omit(DEG)
  
  nrDEG=DESeq2_DEG[,c(2,6)]
  colnames(nrDEG)=c('log2FoldChange','pvalue')  
  
  #set logFC cutoff-value
  logFC_cutoff = 1
  nrDEG$change = as.factor(ifelse(nrDEG$pvalue < 0.05 & abs(nrDEG$log2FoldChange) > logFC_cutoff,
                                     ifelse(nrDEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))
  
  nrDEG$Group = "not-significant"
  nrDEG$Group[which( (nrDEG$pvalue < 0.05) & (nrDEG$log2FoldChange > 1) )] = "up-regulated"
  nrDEG$Group[which( (nrDEG$pvalue < 0.05) & (nrDEG$log2FoldChange < -1) )] = "down-regulated"
  nrDEG$pvalue = -log10(nrDEG$pvalue)
  
  #Plot Volcano Plot
  p1 <-  ggscatter(nrDEG, 
                  x = "log2FoldChange", 
                  y ="pvalue",
                  ylab="-log10(P.value)",
                  size=1.5,
                  color = "change",
                  palette = c("#2f5688","#BBBBBB","#CC0000") 
  )+scale_y_continuous(breaks = c(0,5,10,15),
                       labels = c(0,5,10,15),
                       limits = c(0, 15)) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.8,alpha = 1 / 2) +
    geom_hline(yintercept = 1.30102,lty=4,col="grey",lwd=0.8,alpha = 1 / 2) +
    theme_classic() 
  p1
  #Print boxplot as pdf
  ggsave(p1,filename ='./Figures/step06_DEseq2_volcano.pdf')

  