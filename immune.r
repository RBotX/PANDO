setwd("/home/dan/Thesis")
source("helper.R")
source("PANDO.r")
source("PANDO2.r")
library(MASS)
library(grpreg)
library(xtable)
library(knitr)

load("/home/dan/Downloads/ImmuneSystemDiseases.rda")

data = c()
for(i in 1:length(x1.list)){
  x=rbind(x1.list[[i]],x2.list[[i]])
  y=c(rep(1,nrow(x1.list[[i]])),rep(0,nrow(x2.list[[i]])))
  dat=cbind(x,y)
  dat = cbind(dat, matrix( paste0("fam",i), ncol=1))
  data = rbind(data,dat)
  #write.table(dat,file=paste(diseases[i],".xls",sep=""),row.names=F,col.names=F,sep="\t")
}
colnames(data)=c(paste0("f",1:1809),"Label","Family")



schools2=data

alltests=c()
NUM_TESTS=1
for(l in 1:NUM_TESTS){
  #set.seed(l+1)
  set.seed(l+11)
  testidx = c()
  validx = c()
  trainidx = c()
  
  for(fam in unique(schools2[,"Family"])){
    testidx = c(testidx, sample(which(schools2[,"Family"]==fam),floor(length(which(schools2[,"Family"]==fam))*0.3) ))
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
  
  controls=rpart.control(cp=0.0001,maxdepth = 5)
  iter=5000
  rate=0.01
  ridge.lambda=0.0001  
  #data=GenerateData(d=d,ntrain=ntrain,ntest=ntest,seed=i)
  
  
  # cat("starting PCA...\n")
  # X.pca = prcomp(data$data[,-which(colnames(data$data) %in% c("Family","Label"))],
  #                center = TRUE,
  #                scale = TRUE)
  # PoV <- X.pca$sdev^2/sum(X.pca$sdev^2)
  # PoV=cumsum(PoV)
  # npc = min(which(PoV > 0.80)) # keep all PCs which will get us to 90% of explained variance
  # data.pca = predict(X.pca)[,1:npc]
  # cat("done with PCA...\n")
  # 
  # data$data = cbind(data.pca,data$data[,c("Label","Family")])
  
  
  
  train = data$data[data$trainidx,]
  test = data$data[data$testidx,]

  equalTrainPerFamily=FALSE
  if(equalTrainPerFamily){
    perFamilyTrain = min(table(train[,"Family"]))
    for(fam in unique(train[,"Family"])){
      trainSamplesInFamily = table(train[,"Family"])[fam]
      if(trainSamplesInFamily==perFamilyTrain){
        next
      }
      famidxs = which(train[,"Family"]==fam)
      extractedidxs = sample(famidxs, trainSamplesInFamily-perFamilyTrain)
      test = rbind(test, train[extractedidxs,])
      train = train[-extractedidxs,]
    }
  }
  
  val = data$data[data$validx,]

  perTaskModels=list()
  logitModels=list()
  for(fam in unique(train[,"Family"])){
    
    cat("fam ",fam,"\n")
    tr = train[train[,"Family"]==fam,]
    tr.val = val[val[,"Family"]==fam,]
    # m0 = TrainMultiTaskClassificationGradBoost(tr,valdata=tr.val,groups = matrix(fam,nrow=nrow(tr),ncol=1),iter=iter,v=0.01,
    #                                            controls=controls, ridge.lambda = ridge.lambda,target="binary")  
    # perTaskModels[[toString(fam)]]=m0
    logitModels[[toString(fam)]]= cv.glmnet(x=as.matrix(tr[,-which(colnames(tr) %in% c("Family","Label"))]),y=tr[,"Label"],family="binomial",alpha=1,maxit=10000,nfolds=5, thresh=1E-4,type.measure="deviance")
  }
  
  ### train binary model, ignoring multi tasking:
  binaryData = train
  binaryData["Family"]="1"
  binaryVal = val
  binaryVal["Family"]="1"

  mlogitbinary = cv.glmnet(x=as.matrix(rbind(binaryData,binaryVal)[,-which(colnames(tr) %in% c("Family","Label"))]),y=rbind(binaryData,binaryVal)[,"Label"],family="binomial",alpha=1,maxit=10000,nfolds=5, thresh=1E-4,type.measure="deviance",nlambda=50)
  
  cat("create gplasso train data\n")
  gplassotraindata = CreateGroupLassoDesignMatrix(train,interceptGrouped=FALSE)
  gplassoX = (gplassotraindata$X)[,-ncol(gplassotraindata$X)]
  gplassoy =  (gplassotraindata$X)[,ncol(gplassotraindata$X)]
  gplassoy[gplassoy==-1]=0
  gplassotraindata$X=gplassotraindata$X[sample(nrow(gplassotraindata$X)),]
  mgplasso = cv.grpreg(gplassoX, gplassoy, group=gplassotraindata$groups, nfolds=5, nlambda=50,max.iter=10000,seed=777,family="binomial",trace=TRUE,penalty="grLasso")
  
  vgplassomodels = list()
  lambdas = (0:3)/3
  for(lambda in lambdas){
    vgplassotraindata = CreateVibratingGroupLassoDesignMatrix(train,interceptGrouped=FALSE,binaryVariablesLambda=lambda)
    vgplassoX = (vgplassotraindata$X)[,-ncol(vgplassotraindata$X)]
    vgplassoy =  (vgplassotraindata$X)[,ncol(vgplassotraindata$X)]
    vgplassoy[vgplassoy==-1]=0
    vgplassotraindata$X=vgplassotraindata$X[sample(nrow(vgplassotraindata$X)),]
    
    mvgplasso = cv.grpreg(vgplassoX, vgplassoy, group=vgplassotraindata$groups, nfolds=5, nlambda=50,max.iter=10000,seed=777,family="binomial",trace=TRUE,penalty="grLasso")    
    vgplassomodels[[toString(lambda)]]=mvgplasso
  }
  
  bestLambda=-1
  bestCVE = 1
  for(lambda in lambdas){
    newCVE =  min(vgplassomodels[[toString(lambda)]]$cve)
    cat(newCVE,"\n")
    if(newCVE < bestCVE){
      bestCVE=newCVE
      bestLambda=lambda
    }
  }
  
  rm(gplassotraindata) # remove gplasso train data, it can be very large
  gplassotestdata = CreateGroupLassoDesignMatrix(test)
  gplassotestX = (gplassotestdata$X)[,-ncol(gplassotestdata$X)]
  gplassotesty =  (gplassotestdata$X)[,ncol(gplassotestdata$X)]
  
  gplassoPreds = predict(mgplasso,gplassotestX,type="response",lambda=mgplasso$lambda.min)
  rm(gplassotestdata) # remove gplasso test data, it can be very large
  #methods = c("PANDO","PTB","BB","PTLogit","BinaryLogit","PANDO2","GL")#"PANDO3","PANDO4")#"PANDO3","GL","PANDO4")
  methods = c("PTLogit","BinaryLogit","GL")#"PANDO3","PANDO4")#"PANDO3","GL","PANDO4")
  
  rc=list()
  tt=list()
  compmat = c()
  digitsfmt = matrix(-2,nrow=length(methods),ncol=length(methods))
  
  
  allpreds = matrix(nrow=nrow(test),ncol=length(methods)+2)
  colnames(allpreds)=c(methods,"Label","testnum")
  allpreds[,"Label"]=test[,"Label"]
  allpreds[,"testnum"]=l
  
  
  k=0  
  xtables=list()
  ##################### test:
  for(fam in unique(test[,"Family"])){
    k = k+1
    testidxs = which(test["Family"]==fam)
    compmatrix = matrix(NA,nrow=length(methods),ncol = length(methods))
    
    tr.test = test[test["Family"]==fam,]
    tr.test = tr.test[,-which(colnames(tr.test)=="Family")]
    
    
    bestIt=NULL
    # tt[[methods[which(methods=="PANDO")]]]= predict(mshared[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE,bestIt=bestIt)
    # tt[[methods[which(methods=="PTB")]]]= predict(perTaskModels[[toString(fam)]][[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE)
    # 
    # 
    # tt[[methods[which(methods=="BB")]]] = predict(mbinary[[toString(1)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE,bestIt=bestIt)
    tt[[methods[which(methods=="PTLogit")]]] =predict(logitModels[[toString(fam)]],newx=as.matrix(tr.test[,-which(colnames(tr) %in% c("Family","Label"))]),type="response",s=logitModels[[toString(fam)]]$lambda.min)
    tt[[methods[which(methods=="BinaryLogit")]]] =predict(mlogitbinary,newx=as.matrix(tr.test[,-which(colnames(tr) %in% c("Family","Label"))]),type="response",s=mlogitbinary$lambda.min)
    
    
    # tt[[methods[which(methods=="PANDO2")]]] = predict(mshared2[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE,bestIt=bestIt)
    #tt[[methods[7]]] = predict(mshared3[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE)
    tt[[methods[which(methods=="GL")]]] = gplassoPreds[test[,"Family"]==fam]
    
    #bestIt=min(which(as.vector(mshared3$log$vscore)==max(as.vector(mshared3$log$vscore))))    
    #tt[[methods[which(methods=="PANDO3")]]] = predict(mshared3[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE,bestIt=bestIt)
    
    #bestIt=min(which(as.vector(mshared4$log$vscore)==max(as.vector(mshared4$log$vscore))))    
    #tt[[methods[which(methods=="PANDO4")]]] = predict(mshared4[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE,bestIt=bestIt)
    #tt[[methods[9]]] = predict(mshared4[[toString(fam)]],tr.test[,-which(colnames(tr.test) %in% c("Family","Label"))],calibrate=TRUE)
    
    
    
    
    for(i in 1:length(methods)){
      rc[[methods[i]]]=pROC::roc(as.factor(tr.test[,"Label"]),as.numeric(tt[[methods[i]]]))
      allpreds[testidxs,methods[i]] = as.matrix(tt[[methods[i]]],ncol=1)
    }
    
    for(i in 1:length(methods)){
      compmatrix[i,i]=rc[[methods[i]]]$auc[1]
      digitsfmt[i,i]=6
      for(j in 1:length(methods)){
        if(i >=j ){
          next
        }
        #cat("setting compmatrix",i," ",j,"\n")
        compmatrix[i,j] = pROC::roc.test(rc[[methods[i]]],rc[[methods[j]]],paired=TRUE)$p.value
        cat("auc  for ",fam," ", methods[i]," VS ",methods[j],": ",round((rc[[methods[i]]]$auc[1]),4)," ",round((rc[[methods[j]]]$auc[1]),4)," with pval: ",pROC::roc.test(rc[[methods[i]]],rc[[methods[j]]],paired=TRUE)$p.value,"\n")
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
  
  alltests=rbind(alltests,allpreds)
  cat("round ",l," summary:\n********************\n")
  for(method in methods){
    score=pROC::roc(as.factor(alltests[alltests[,"testnum"]==l,"Label"]),alltests[alltests[,"testnum"]==l,method])$auc[1]
    cat(method," ",score,"\n")
  }
  
  
  
}

cat("\n*********************************************\n")
cat("final results:\n")
finalresults= matrix(nrow=NUM_TESTS,ncol=length(methods))
#finalAUC= matrix(nrow=NUM_TESTS,ncol=length(methods))
colnames(finalresults)=methods
#colnames(finalAUC)=methods
for(l in 1:NUM_TESTS){
  for(method in methods){
    
    score=pROC::roc(as.factor(alltests[alltests[,"testnum"]==l,"Label"]),alltests[alltests[,"testnum"]==l,method])
    cat(method," ",score$auc[1],"\n")
    finalresults[l,method]=score$auc[1]
    #finalAUC[l,method]=score
  }
}
for(m in colnames(finalresults)){
  cat(m ,"mean:",mean(finalresults[,m]),"std:",sd(finalresults[,m]), "mean+std:",mean(finalresults[,m])-sd(finalresults[,m]),"\n")
}


compmat2=matrix(NA,nrow=length(unique(test[,"Family"])),ncol=length(methods))
rownames(compmat2)=unique(test[,"Family"])
colnames(compmat2)=methods
for(fam in unique(test[,"Family"])){
  testidxs = (test["Family"]==fam)
  for(method in methods){
    compmat2[which(rownames(compmat2)==fam),which(colnames(compmat2)==method)]=roc(as.factor(alltests[(alltests[,"testnum"]==l)&(testidxs),"Label"]),alltests[(alltests[,"testnum"]==l)&(testidxs),method])$auc[1]
    
  }
}


#### seed=1 or 2 got BB vs PANDO2 with 
# Z = -2.2629, p-value = 0.02364
# alternative hypothesis: true difference in AUC is not equal to 0
# sample estimates:
#   AUC of roc1 AUC of roc2 
# 0.9796171   0.9823319 

roc.test(roc(as.factor(alltests[alltests[,"testnum"]==l,"Label"]),alltests[alltests[,"testnum"]==l,"GL"]),
         roc(as.factor(alltests[alltests[,"testnum"]==l,"Label"]),alltests[alltests[,"testnum"]==l,"BinaryLogit"]))



#save.image(compress=TRUE)



for(m in methods){
  cat(m," ",mean(compmat2[,m]),"\n")
}





#################### interpretation gplasso
plot(mgplasso)
plot(mgplasso$lambda,mgplasso$cve,main="group lasso CV err as function of lambda")
numOfGroupsSelected = sum(mgplasso$fit$beta[,mgplasso$min][2:1812]!=0)
#mgplasso$fit$beta[,mgplasso$min]
#length(as.matrix(mgplasso$fit$beta[,mgplasso$min],ncol=1)[-1,])
d=ncol(train)-2 ### ignore label/family
#betas=mgplasso$fit$beta[,mgplasso$min][-1]
betas=mgplasso$fit$beta[,mgplasso$min]
taskBetas = matrix(NA,ncol=5,nrow=d)
taskBetas2 = c()
for(i in 1:5){
  ##betas for i-th task

  start = (d*(i-1))+1
  end = start+d-1  
  
  start = start+i
  end = end+i
  
  #cat("i is ",i,"\n","start ",start+( (i>1)+0), "end ",end+((i>1)+0), "\n*****************\n")
  cat("i is ",i,"\n","start ",start, "end ",end, "\n*****************\n")
  taskBetas[,i] = mgplasso$fit$beta[,mgplasso$min][start:end]
  taskBetas2 = rbind(taskBetas2,cbind(mgplasso$fit$beta[,mgplasso$min][start:end],1:d))
  #cat("length is ", length(mgplasso$fit$beta[,mgplasso$min][start:end]),"\n")
  #cat(which(is.na(mgplasso$fit$beta[,mgplasso$min][start:end])))
}
taskBetas2=matrix(taskBetas2,ncol=2)
taskBetas2=data.frame(taskBetas2)
colnames(taskBetas2)=c("varValue","varName")

if(length(which(is.na(taskBetas)))!=0){
  cat("ERROR! some coefficients from group lasso reported as NA\n")
}




require(ggplot2)
taskBetas2[,"color"]="1"
taskBetas2=taskBetas2[taskBetas2["varValue"]!=0,]
ggplot(data = taskBetas2, aes(x=factor(varName), y=varValue)) + geom_boxplot(aes(fill=color)) +
  labs(x="covariate index", y="coef value boxplot over tasks", title="gplasso coef boxplots")



bestBinaryLogitLambda=which(mlogitbinary$lambda==mlogitbinary$lambda.min)
binarybetas=mlogitbinary$glmnet.fit$beta[,bestBinaryLogitLambda]




taskBetas[taskBetas==0]=NA ## don't plot zero variables
binarybetas[binarybetas==0]=NA ## don't plot zero variables

ptsize=1
pchval=21
plot(1:d,taskBetas[,1],cex=ptsize,type="p",pch=pchval, bg="lightgreen",xlab="variable index",ylab="coefficient values",main="Binary Logit VS gplasso coefs")
lines(1:d,taskBetas[,2],cex=ptsize,type="p",pch=pchval, bg="lightgreen")
lines(1:d,taskBetas[,3],cex=ptsize,type="p",pch=pchval, bg="lightgreen")
lines(1:d,taskBetas[,4],cex=ptsize,type="p",pch=pchval, bg="lightgreen")
lines(1:d,taskBetas[,5],cex=ptsize,type="p",pch=pchval, bg="lightgreen")
lines(1:d,binarybetas,cex=ptsize,type="p",pch=pchval,bg="lightblue")

legend("bottomright", legend=c("group lasso", "binary logit"),
       col=c("lightgreen", "lightblue"),  cex=0.8,fill=c("lightgreen", "lightblue"))

