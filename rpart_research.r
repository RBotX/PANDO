# Classification Tree with rpart
library(rpart)

# grow tree 
fit <- rpart(Kyphosis ~ Age + Number + Start,
             method="class", data=kyphosis)

predict(fit)
newfit=fit
newframe = newfit$frame[,"yval2"]
newframe[,4]=1
newframe[,5]=0
newfit$frame[,"yval2"]=newframe
predict(newfit)


unique(fit$where)
dim(fit$frame)
rownames(fit$frame)



library(foreach)
library(doParallel)
library(registerDoMC)
cl <- makeCluster(3)
registerDoParallel(cl)

foreach(i=4:7) %dopar% sqrt(i)