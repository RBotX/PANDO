library(gbm)
# load mtcars data
library(dummies)
require(pROC)
library(dplyr)
library(plotmo)
library(ElemStatLearn)
#data(mtcars)
# fit GBM
#df=mtcars


# target="loss"
df=read.csv("/home/dan/Thesis/convergence/spambase.data",stringsAsFactors=TRUE,header = FALSE) 
colnames(df)[ncol(df)]="y" 
subsampidx = sample(1:nrow(df), size = round(0.30*nrow(df)), replace=FALSE)
df=df[subsampidx,]
rownames(df)=NULL

trainIndex = sample(1:nrow(df), size = round(0.70*nrow(df)), replace=FALSE)
train = df[trainIndex ,]
test = df[-trainIndex ,]
ntrees=3000
gbmFit<-gbm(
             y~.,
             data=train,
             distribution = "bernoulli",
             interaction.depth=1,
             shrinkage = 0.01,
             bag.fraction=1,
             train.fraction = 0.8,
             verbose = TRUE,
             n.trees = ntrees)

gbmpred<-predict(gbmFit,test[,-which(colnames(test)=="y")],n.trees=ntrees,type="response")

stumps=c()
for(i in 1:ntrees){
  
  ### the variables from pretty.gbm.tree are zero indexed
  #stumps = rbind(stumps,c(pretty.gbm.tree(gbmFit,i.tree=i)[1,'SplitVar'], pretty.gbm.tree(gbmFit,i.tree=i)[1,'SplitCodePred']))
  stumps = rbind(stumps,c(pretty.gbm.tree(gbmFit,i.tree=i)[1,'SplitVar'], pretty.gbm.tree(gbmFit,i.tree=i)[1,'SplitCodePred']))
}


### generate matrix of splits for each variable
zs = stumps[!duplicated(stumps),] ## count distinct stumps
## for each type of stump, create the 
k=1
XX = matrix(NA,nrow=nrow(df),ncol= nrow(zs)*2)
colnames(XX) = 1:ncol(XX)
for(i in 1:nrow(zs)){
  splitVar=(zs[i,1])
  splitVal=zs[i,2]

  for(j in 1:nrow(df)){
    XX[j,k] = if(df[j,splitVar]<=splitVal) 1 else 0
    XX[j,k+1] = if(df[j,splitVar]>splitVal) 1 else 0  
  }
  
  colnames(XX)[k] = paste0("cart_svar_",splitVar,"_sval_",splitVal,"_1")
  colnames(XX)[k+1] = paste0("cart_svar_",splitVar,"_sval_",splitVal,"_2")
  k = k+2
}

### generate more decision stumps variables based on deciles:

ntile=9
#splitVars = unique(zs[,1])
splitVars = 1:(ncol(df)-1) ### try to split on all variables
YY=matrix(NA,nrow=nrow(df),ncol=length(splitVars)*2*(ntile+1) )
colnames(YY)=1:ncol(YY)
l=1
for(splitVar in splitVars){
    splitVals = as.numeric(quantile(df[,splitVar], seq(0, 1, 1/ntile))) 
    
    
    for(j in 1:nrow(df)){
      if(j %% 100 == 0){
        cat("finished ",j," row out of ",nrow(df),"\n")
        cat("l is ",l,"\n")
      }
      tile=0
      for(splitVal in splitVals){
        YY[j,l+tile] = if(df[j,splitVar]<=splitVal) 1 else 0
        YY[j,l+tile+1] = if(df[j,splitVar]>splitVal) 1 else 0  
        #cat("inserted ",j," ",l+tile,"\n")
        colnames(YY)[l+tile] = paste0("tree_svar_",splitVar,"sval_",splitVal,"_",tile/2,"_1")
        colnames(YY)[l+tile+1] = paste0("tree_svar_",splitVar,"sval_",splitVal,"_",tile/2,"_2")
        tile = tile+2
      }
    }
    l = l+(2*(ntile+1)) # each variable adds that many columns
}

# XX = cbind(XX,YY)
# XX = cbind(XX,df[,"y"])
# colnames(XX)[ncol(XX)]="y"
# 
# 
# XX.train = XX[trainIndex,]
# XX.test = XX[-trainIndex,]




YY = cbind(YY,df[,"y"])
colnames(YY)[ncol(YY)]="y"


YY.train = YY[trainIndex,]
YY.test = YY[-trainIndex,]

## lasso
library(glmnet)
cvfit.lasso<-cv.glmnet(YY.train[,-which(colnames(YY.train)=="y")], YY.train[,"y"], alpha=1,family="binomial") 
bestLambda.lasso = cvfit.lasso$lambda.min
fit.lasso = glmnet(YY.train[,-which(colnames(YY.train)=="y")], YY.train[,"y"], alpha=1,lambda=bestLambda.lasso,family="binomial")
coefs.lasso=as.matrix(coef(fit.lasso))

### take only non zero coefficients, remove intercept
nzlasso = coefs.lasso[coefs.lasso != 0]
nzlasso = as.matrix(nzlasso)
nzlasso=matrix(nzlasso,ncol=2,nrow=nrow(nzlasso))
nzlasso[,1]=rownames(coefs.lasso)[coefs.lasso!=0]
nzlasso=nzlasso[-1,]
###

### create a dataframe which describes the splits taken by lasso:
### one column for the variable index, second column for the split value
lassoVars = c()
for(i in 1:nrow(nzlasso)){
  
  v=nzlasso[i,1]
  
  lassoVars = rbind(lassoVars,c(as.numeric(strsplit(strsplit(v,"_")[[1]][3],"sval")[[1]][1]),as.numeric(strsplit(v,"_")[[1]][4]),as.numeric(nzlasso[i,2])))
}
colnames(lassoVars) = c("varIndex","splitValue","coeff")
colnames(zs) = c("varIndex","splitValue")



lbs_fun <- function(fit, ...) {
  L <- length(fit$lambda)
  x <- log(fit$lambda[L])
  y <- fit$beta[, L]
  labs <- names(y)
  text(x, y, labels=labs, ...)
}

plot(cvfit$glmnet.fit, "lambda",   label=TRUE)
lbs_fun(fit)




### ridge
cvfit.ridge<-cv.glmnet(YY.train[,-which(colnames(YY.train)=="y")], YY.train[,"y"], alpha=0,family="binomial") # lasso
bestLambda.ridge = cvfit.ridge$lambda.min
fit.ridge = glmnet(YY.train[,-which(colnames(YY.train)=="y")], YY.train[,"y"], alpha=0,lambda=bestLambda.ridge,family="binomial")
coefs.ridge=as.matrix(coef(fit.ridge))



# as.matrix(coef(mm))[as.matrix(coef(mm))[,1]>0]
# colnames(YY)[as.matrix(coef(mm))[,1]>0]
# lass=predict(mm,newx = YY.test[,-which(colnames(YY.test)=="y")],type="response",s=mm$lambda.min)
# 
# mean(abs(gbmpred-test[,"y"]))
# mean(abs(linp-YY.test[,"y"]))



mm2<-cv.glmnet(as.matrix(df[,-which(colnames(df)=="y")]), as.matrix(df[,"y"]), alpha=0) #ridge
lin2=predict(mm2,newx = as.matrix(df[,-which(colnames(df)=="y")]),type="response",s=mm$lambda.min)



plot(mm$glmnet.fit, "lambda", label=FALSE)