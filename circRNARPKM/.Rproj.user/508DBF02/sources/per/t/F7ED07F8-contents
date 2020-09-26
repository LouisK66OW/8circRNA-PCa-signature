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
if(T){
  library(ggplot2)
  library(glmnet)
  library(ROCR)
  library(caret)
  library(ggsci)
  library(ggpubr) 
  library(survival)
  library(survminer)
}

load( file ="./Rdata/PRAD-circRNA-sample.Rdata")


#Cross-validation with ten-folds
set.seed(2)
x1 <- t(exprSet)
meta <- phe
BCR <- meta$BCR
TBCR <- meta$TBCR
datafra <- data.frame(x1,BCR,TBCR)
require(caret)
folds <- createFolds(y=datafra[,29],k=10)

# for-loop structure to caculate with lasso
auc_valuemin <-as.numeric()
auc_value1se <- as.numeric()
  for(i in 1:10){
    seedi = i + 11
    set.seed(seedi)
    fold_test1 <- datafra[folds[[i]],] #folds[[i]]as Test Set
    fold_test  <- fold_test1[,1:28]
    fold_train1 <- datafra[-folds[[i]],] # The rest as Validation set 
    fold_train <- fold_train1[,1:28]
    x <- log2(fold_train+1)
    x <- as.matrix(x)
    y <- fold_train1$BCR
    cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
    lambda_min1se <- data.frame(cv_fit$lambda.min,cv_fit$lambda.1se) 
    model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
    model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
    choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
    choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
    fold_test <- as.matrix(log2(fold_test+1))
    lasso.prob <- predict(cv_fit, newx=fold_test, s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
    re=cbind(fold_test1$BCR,lasso.prob)
    re=as.data.frame(re)
    colnames(re)=c('event','prob_min','prob_1se')
    re$event=as.factor(re$event)
    # Training set model predicts test set
    #min
    pred_min <- prediction(re[,2], re[,1])
    auc_min = performance(pred_min,"auc")@y.values[[1]]
    #1se
    pred_1se <- prediction(re[,3], re[,1])
    auc_1se = performance(pred_1se,"auc")@y.values[[1]]
    auc_valuemin <- append(auc_valuemin,auc_min)
    auc_value1se <- append(auc_value1se,auc_1se)
  }
#The auc values by min or 1se
print(auc_valuemin)

print(auc_value1se)

#Obtain the signatures ------

num <- as.numeric(4)
fold_test1 <- datafra[folds[[num]],] 
fold_test  <- fold_test1[,1:28]
fold_train1 <- datafra[-folds[[num]],] 
fold_train <- fold_train1[,1:28]
x <- log2(fold_train+1)
x <- as.matrix(x)
y <- fold_train1$BCR

set.seed(num+11)
model_lasso <- glmnet(x, y, family="binomial", nlambda=50, alpha=1)
print(model_lasso)
head(coef(model_lasso, s=c(model_lasso$lambda[29],0.009)))
plot(model_lasso, xvar="lambda", label=TRUE)
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
plot(cv_fit)
lambda_min1se <- data.frame(cv_fit$lambda.min,cv_fit$lambda.1se) 
#obtain the signature and coefficient by min or 1se
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_min
model_lasso_min$beta[as.numeric(model_lasso_min$beta)!=0]

model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
choose_gene_1se
model_lasso_1se$beta[as.numeric(model_lasso_1se$beta)!=0]

#Riskscore----
#Kaplan-Meier curve between Risk Score and Biochemical Recurrence
load(  file = 
         file.path('Rdata/','TCGA_RPKM_8circRNA.Rdata')
)
library(survival)
library(survminer)
sfit <- survfit(Surv(TBCR,BCR )~Rscore_group, data=meta)
sfit
summary(sfit)

pdf(file = "./Figure/step02LASSO/Risk_score_KMcurve.pdf")
ggsurvplot(sfit,palette = c("#CC0000", "#2f5688"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
dev.off()

#Kaplan-Meier curves of each circRNA in signature. 
gs=c('circ_30029','circ_117300','circ_176436','circ_112897', 'circ_17720',
     'circ_178252','circ_115617','circ_14736') 
splots <- lapply(gs, function(g){
  phe$gene=ifelse(as.numeric(exprSet[g,])>median(as.numeric(exprSet[g,])),'high','low')
  table(phe$gene)
  sfit1=survfit(Surv(TBCR, BCR)~gene, data=phe)
  ggsurvplot(sfit1,pval =TRUE, data = phe, risk.table = TRUE,
             palette=c("#CC0000", "#2f5688"))+
    labs(title = g)
})
pdf(file = "./Figure/step02LASSO/8circs_KMcurve.pdf")
splots
dev.off()