library(ElemStatLearn)
library(glmnet)
data(prostate)
target="lpsa"
x=prostate[prostate$train==TRUE,-which(colnames(prostate) %in% c(target,"train"))]
y=prostate[prostate$train==TRUE,target]
glmmod <- glmnet(as.matrix(x), y=as.matrix(y), alpha=1)

# Plot variable coefficients vs. shrinkage parameter lambda.
plot(glmmod, xvar="lambda")
plot(glmmod) ## similar to thesis version
colors=c("red","blue","cyan","green","yellow","black","orange","purple")

plot(glmmod,col=colors)

beta=glmmod$beta
beta=t(beta)
colnames(beta)=colnames(x)
lambdas=glmmod$lambda
lambdas = log(lambdas)

plot(x=lambdas,y=beta[,1],xlim=rev(range(lambdas)),ylim=range(c(-0.2,beta[,1])),xlab="log(lambda)",ylab="coefficients",type="l",col=colors[1],cex=3.5)
lines(x=lambdas,y=beta[,2],xlim=rev(range(lambdas)),col=colors[2])
lines(x=lambdas,y=beta[,3],xlim=rev(range(lambdas)),col=colors[3])
lines(x=lambdas,y=beta[,4],xlim=rev(range(lambdas)),col=colors[4])
lines(x=lambdas,y=beta[,5],xlim=rev(range(lambdas)),col=colors[5])
lines(x=lambdas,y=beta[,6],xlim=rev(range(lambdas)),col=colors[6])
lines(x=lambdas,y=beta[,7],xlim=rev(range(lambdas)),col=colors[7])
lines(x=lambdas,y=beta[,8],xlim=rev(range(lambdas)),col=colors[8])


text(locator(), labels = colnames(beta))
