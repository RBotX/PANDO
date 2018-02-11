setwd("/home/dan/Thesis")
source("helper.R")
source("PANDO.r")
source("PANDO2.r")
source("pando3.r")
library(MASS)
library(grpreg)
library(xtable)
library(knitr)
library(R.matlab)


   
schools2 = readMat("/home/dan/Thesis/microarray/Split_Data_Microarray_Ratio=0.6.mat")
X=c()
for(i in 1:18){
  x = data.frame(schools2$DX)
  x = cbind(x,matrix(rep(paste0("fam",i),nrow(x)),ncol=1))
  x = cbind(x,schools2$DY[,i]) # take i-th family
  X = rbind(X,x)
}
schools2=X
colnames(schools2)=c(paste0("f",1:21),"Family","Label")
schools2=data.frame(schools2)


# alltests=c()
# NUM_TESTS=5
# for(l in 1:NUM_TESTS){
  #set.seed(l+1)
set.seed(l+10)
testidx = c()
validx = c()
trainidx = c()
for(fam in unique(schools2[,"Family"])){
  
  testidx = c(testidx, sample(which(schools2[,"Family"]==fam),floor(length(which(schools2[,"Family"]==fam))*0.2) ))
  famtrain = setdiff(which(schools2[,"Family"]==fam),testidx)
  famvalidx = sample(famtrain,0.2*length(famtrain) ) #val 10% of train
  validx = c(validx,famvalidx)
  famtrain = setdiff(famtrain,validx) # remove validation from train
  trainidx = c(trainidx,famtrain)
}



data=list()
data[["data"]]=schools2
data[["testidx"]]=testidx
data[["trainidx"]]=trainidx
data[["validx"]]=validx


iter=1000
rate=0.01
ridge.lambda=1  
#data=GenerateData(d=d,ntrain=ntrain,ntest=ntest,seed=i)

train = data$data[data$trainidx,]
test = data$data[data$testidx,]
val = data$data[data$validx,]

mshared=TunePando(TrainMultiTaskClassificationGradBoost3,train,val,fitTreeCoef="nocoef",fitLeafCoef="ridge",trainrate = 0.1,target="regression",maxdepths = c(10))

perTaskModels=list()
logitModels=list()
for(fam in unique(train[,"Family"])){
  
  cat("fam ",fam,"\n")
  tr = train[train[,"Family"]==fam,]
  tr.val = val[val[,"Family"]==fam,]
  #    m0 = TrainMultiTaskClassificationGradBoost(tr,valdata=tr.val,groups = matrix(fam,nrow=nrow(tr),ncol=1),iter=iter,v=0.01,
  #                                               controls=rpart.control(), ridge.lambda = ridge.lambda,target="binary") 
  m0=TunePando(vanillaboost,tr,tr.val,target="regression")
  # m0 = TrainMultiTaskClassificationGradBoost(tr,valdata=tr.val,groups = matrix(fam,nrow=nrow(tr),ncol=1),iter=iter,v=0.01,
  #                                            controls=rpart.control(), ridge.lambda = ridge.lambda,target="binary")  
  
  perTaskModels[[toString(fam)]]=m0
  logitModels[[toString(fam)]]= cv.glmnet(x=as.matrix(rbind(tr,tr.val)[,-which(colnames(tr) %in% c("Family","Label"))]),y=rbind(tr,tr.val)[,"Label"],family="gaussian",alpha=1,maxit=10000,nfolds=3, thresh=1E-4)
}

### train binary model, ignoring multi tasking:
binaryData = train
binaryData["Family"]="1"
binaryVal = val
binaryVal["Family"]="1"
mbinary=TunePando(vanillaboost,binaryData,binaryVal,target="regression")
mlogitbinary = cv.glmnet(x=as.matrix(rbind(binaryData,binaryVal)[,-which(colnames(tr) %in% c("Family","Label"))]),y=rbind(binaryData,binaryVal)[,"Label"],family="gaussian",alpha=1,maxit=10000,nfolds=4, thresh=1E-4)

gplassotraindata = CreateGroupLassoDesignMatrix(rbind(train,val))
gplassoX = (gplassotraindata$X)[,-ncol(gplassotraindata$X)]
gplassoy =  (gplassotraindata$X)[,ncol(gplassotraindata$X)]
#gplassoy[gplassoy==-1]=0

mgplasso = cv.grpreg(gplassoX, gplassoy, group=gplassotraindata$groups, nfolds=2, maxit=10000,seed=777,family="gaussian",trace=TRUE,penalty="grLasso")

gplassotestdata = CreateGroupLassoDesignMatrix(test)
gplassotestX = (gplassotestdata$X)[,-ncol(gplassotestdata$X)]
gplassotesty =  (gplassotestdata$X)[,ncol(gplassotestdata$X)]

gplassoPreds = predict(mgplasso,gplassotestX,type="response",lambda=mgplasso$lambda.min)

methods = c("PANDO","PTB","BB","PTLogit","BinaryLogit","GL")#,"PANDO3","GL")

rc=list()
tt=list()
compmat = c()
digitsfmt = matrix(-2,nrow=length(methods),ncol=length(methods))
xtables=list()
k=0
allpreds = matrix(nrow=nrow(test),ncol=length(methods)+2)
colnames(allpreds)=c(methods,"Label","testnum")
allpreds[,"Label"]=test[,"Label"]
allpreds[,"testnum"]=l
##################### test:
for(fam in unique(test[,"Family"])){
  k = k+1
  testidxs = which(test["Family"]==fam)
  compmatrix = matrix(NA,nrow=length(methods),ncol = length(methods))
  
  tr.test = test[test["Family"]==fam,]
  tr.test = tr.test[,-which(colnames(tr.test)=="Family")]
  
  bestIt=min(which(as.vector(mshared$log$vscore)==max(as.vector(mshared$log$vscore))))    
  tt[[methods[which(methods=="PANDO")]]]= predict(mshared[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE,bestIt=bestIt)
  tt[[methods[which(methods=="PTB")]]]= predict(perTaskModels[[toString(fam)]][[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE)
  
  bestIt=min(which(as.vector(mbinary$log$vscore)==max(as.vector(mbinary$log$vscore))))    
  tt[[methods[which(methods=="BB")]]] = predict(mbinary[[toString(1)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE,bestIt=bestIt)
  tt[[methods[which(methods=="PTLogit")]]] =predict(logitModels[[toString(fam)]],newx=as.matrix(tr.test[,-which(colnames(tr) %in% c("Family","Label"))]),type="response",s=logitModels[[toString(fam)]]$lambda.min)
  tt[[methods[which(methods=="BinaryLogit")]]] =predict(mlogitbinary,newx=as.matrix(tr.test[,-which(colnames(tr) %in% c("Family","Label"))]),type="response",s=mlogitbinary$lambda.min)
  
  bestIt=min(which(as.vector(mshared2$log$vscore)==max(as.vector(mshared2$log$vscore))))    
  tt[[methods[which(methods=="PANDO2")]]] = predict(mshared2[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE,bestIt=bestIt)
  #tt[[methods[7]]] = predict(mshared3[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE)
  tt[[methods[which(methods=="GL")]]] = gplassoPreds[test[,"Family"]==fam]
  
  for(i in 1:length(methods)){
    rc[[methods[i]]]=sqrt(mean((tr.test[,"Label"]-tt[[methods[i]]])**2))
    allpreds[testidxs,methods[i]] = as.matrix(tt[[methods[i]]],ncol=1)
  }
  
  
  
  for(i in 1:length(methods)){
    #cat("storing ",round(pROC::auc(rc[[methods[i]]])[1],4)," in ",i,i,"\n" )
    compmatrix[i,i]=rc[[methods[i]]] 
    digitsfmt[i,i]=3
    for(j in 1:length(methods)){
      if(i >=j ){
        next
      }
      #cat("storing ",signif(pROC::roc.test(rc[[methods[i]]],rc[[methods[j]]])$p.value, digits = 3)," now in ",i,j,"\n")
      compmatrix[i,j] = signif((rc[[methods[i]]]-rc[[methods[j]]]), digits = 3)
      cat("RMSE for ",fam," ", methods[i]," VS ",methods[j],": ",round((rc[[methods[i]]])[1],4)," ",round((rc[[methods[j]]])[1],4),"\n")
    }
  }
  compmat = rbind(compmat,compmatrix)
  cat("***********\n")
  
  
  dft=data.frame(compmatrix)
  colnames(dft)=methods
  rownames(dft)=methods
  xtables[[toString(k)]]=xtable(dft,digits=cbind(rep(1,nrow(digitsfmt)),digitsfmt))
  
}
#####

for(method in methods){
  score=mean(sqrt((tr.test[,"Label"]-tt[[method]])**2))
  cat(method," ",score,"\n")
}

alltests=rbind(alltests,allpreds)
  

cat("*********************************************\n")
cat("final results - RMSE:\n")
finalresults= matrix(nrow=NUM_TESTS,ncol=length(methods))
finalresults2= matrix(nrow=NUM_TESTS,ncol=length(methods))
colnames(finalresults)=methods
colnames(finalresults2)=methods
for(l in 1:NUM_TESTS){
  for(method in methods){
    
    score=mean(sqrt((alltests[alltests[,"testnum"]==l,"Label"]-alltests[alltests[,"testnum"]==l,method])**2))
    cat(method," ",score,"\n")
    finalresults[l,method]=score
  }
}
for(m in colnames(finalresults)){
  cat(m ,"mean:",mean(finalresults[,m]),"std:",sd(finalresults[,m]), "mean+std:",mean(finalresults[,m])+sd(finalresults[,m]),"\n")
}

###exaplined vairance
cat("final results - Explained Variance:\n")
for(l in 1:NUM_TESTS){
  for(method in methods){
    totalvar=var(schools2[,"Label"])
    score=(totalvar - var(alltests[alltests[,"testnum"]==l,"Label"]-alltests[alltests[,"testnum"]==l,method]))/(totalvar)
    
    finalresults2[l,method]=score
  }
}

for(m in colnames(finalresults2)){
  cat(m ,"mean:",mean(finalresults2[,m]),"std:",sd(finalresults2[,m]), "mean-std:",mean(finalresults2[,m])-sd(finalresults2[,m]),"\n")
}





