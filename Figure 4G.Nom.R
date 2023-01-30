library(rms)
library(ROCR)
setwd("C:\\Users\\HASEE\\Desktop\\BMs\\20.test_Nom")
data1=read.table("hub_gene.txt",sep="\t",header = F)
data2=read.table("input.txt",sep="\t",header = T,row.names = 1)
data3=data2[as.vector(data1[,1]),]
data3=t(data3)
data3=data.frame(id=rownames(data3),data3)
data4=read.table("clinical.txt",sep="\t",header = T)
data5=merge(data4,data3,by="id")
write.table(data5,"lg.txt",sep = "\t",quote = F,row.names = F)

inputfile2="lg.txt"
data1<- read.table(inputfile2,header=T,sep="\t",row.names = 1)

for (i in colnames(data1[,2:ncol(data1)])){
  exp=factor(ifelse(data1[,i]<median(data1[,i]),"0","1"))
  data1[,i]=exp
  data1[,i]=factor(data1[,i],labels=c("Low","High"))
}

head(data1)
median(as.numeric(data1[,"MMP7"]))

ddist <- datadist(data1)
options(datadist="ddist") 

mylog<- glm(status~.,family=binomial(link = "logit"),data = data1)  

summary(mylog)
coefficients(mylog)#回归系数
exp(coefficients(mylog))#OR
exp(confint(mylog))#95%CI


#nom
mylog<-lrm(status~.,data=data1,x=T,y=T)


mynom<- nomogram(mylog, fun=plogis,fun.at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),lp=F, funlabel="Risk of IPF")

pdf("Nom.pdf",10,8)
plot(mynom)
dev.off()
#C-index
Cindex <- rcorrcens(data1$status~predict(mylog))
Cindex
#Calibration
mycal<-calibrate(mylog,method="boot",B=1000)

pdf("Calibration.pdf")
plot(mycal,xlab="Nomogram-predicted probability of IPF ",ylab="Actual diagnosed IPF (proportion)",sub=F)
dev.off()

#AUC
pre_rate<-predict(mylog)
ROC1<- prediction(pre_rate,data1$status) 
ROC2<- performance(ROC1,"tpr","fpr")
AUC <- performance(ROC1,"auc")

print(AUC,max = 100000)

AUC<-0.917 #改为自己的AUC

pdf("ROC.pdf")
plot(ROC2,col="blue", xlab="False positive rate",ylab="True positive rate",lty=1,lwd=3,main=paste("AUC=",AUC))
abline(0,1,lty=2,lwd=3)
dev.off() 



