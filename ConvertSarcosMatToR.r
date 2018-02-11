require(R.matlab)
t=readMat("/home/dan/Downloads/sarcos_inv_test.mat")
t=data.frame(t)
m = c()
for(i in 1:7){
  taskdata=cbind(t[,c(1:21,21+i)],rep(paste0("fam",i),nrow(t)))
  colnames(taskdata)=c(paste0("X.",1:21),"Label","Family")
  m = rbind(m,taskdata)
}

write.csv(m,"sarcos_test.csv",row.names = FALSE)

t=readMat("/home/dan/Downloads/sarcos_inv.mat")
t=data.frame(t)
m = c()
for(i in 1:7){
  taskdata=cbind(t[,c(1:21,21+i)],rep(paste0("fam",i),nrow(t)))
  colnames(taskdata)=c(paste0("X.",1:21),"Label","Family")
  m = rbind(m,taskdata)
}

write.csv(m,"sarcos_train.csv",row.names = FALSE)