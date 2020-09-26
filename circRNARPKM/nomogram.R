### ---------------
###
### Create: Shuo Wang
### Date: 2020-04-12
### Email: mywanuo@163.com
### Southern Medical University,Guang Zhou / Major in Urology Surgery 
### Update Log: 2020-09-26 11:17:47  Modified versions
###
### ---------------

rm(list = ls())
options(stringsAsFactors = F)

library(survival)
library(rms)

load(file ='Rdata/TCGA_RPKM_8circRNA_ROC_patch.Rdata')
# "packing" the data in accordance with the nomogram requirements
ddist <- datadist(meta)
options(datadist="ddist")

#Build COX model, draw Nomo diagram
res.cox <- psm(Surv(TBCR, BCR) ~ Age+PSA+`Clinical_T-category`+Gleason+Rscore, data = meta,dist='lognormal')
res.cox <- psm(Surv(TBCR, BCR) ~ PSA+`Clinical_T-category`+Gleason+Rscore, data = meta,dist='lognormal')
surv <- Survival(res.cox)

nom.cox <- nomogram(res.cox,lp=F,
                    fun=function(x) surv(10, x),
                    funlabel="5-year BCR-free Probability",
                    maxscale=10,
                    fun.at=c(0.01,seq(0.1,0.9,by=0.2),0.95))

plot(nom.cox)
pdf("./Figure/step03nomogram/nomogram_noAgeLinerprediction.pdf")
plot(nom.cox)
dev.off()

# Calculate C-index
#Method 1
rcorrcens(Surv(TBCR,BCR) ~ predict(res.cox), data =  meta)
#Method 2
f <- coxph(Surv(TBCR, BCR==1)~PSA+`Clinical_T-category`+Gleason+Rscore, data = meta) 
sum.surv <- summary(f)
c_index <- sum.surv$concordance
c_index

sessionInfo() 

# R version 3.6.2 (2019-12-12)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18362)
