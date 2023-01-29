library(corrplot)                   
expFile="panGeneExp.txt"             
scoreFile="estimateScores.txt"      
scoreType="ESTIMATEScore"             
setwd("E:\\.生信自学网\\压缩包\\panCancer_multiGenes\\15.estimateCor")     


exp=read.table(expFile, header=T,sep="\t",row.names=1,check.names=F)
exp=exp[(exp[,"Type"]=="Tumor"),]


TME=read.table(scoreFile, header=T,sep="\t",row.names=1,check.names=F)


sameSample=intersect(row.names(TME),row.names(exp))
TME=TME[sameSample,]
exp=exp[sameSample,]


outTab=data.frame()
pTab=data.frame()

for(i in levels(factor(exp[,"CancerType"]))){
    exp1=exp[(exp[,"CancerType"]==i),]
    TME1=TME[(TME[,"CancerType"]==i),]
    x=as.numeric(TME1[,scoreType])
    pVector=data.frame(i)
    outVector=data.frame(i)
	
	genes=colnames(exp1)[1:(ncol(exp1)-2)]
	for(j in genes){
		y=as.numeric(exp1[,j])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pValue=corT$p.value
		pVector=cbind(pVector,pValue)
		outVector=cbind(outVector,cor)
	}
	pTab=rbind(pTab,pVector)
	outTab=rbind(outTab,outVector)
}

colNames=c("CancerType",colnames(exp1)[1:(ncol(exp1)-2)])
colnames(outTab)=colNames
write.table(outTab,file="estimateCor.cor.txt",sep="\t",row.names=F,quote=F)

colnames(pTab)=colNames
write.table(pTab,file="estimateCor.pval.txt",sep="\t",row.names=F,quote=F)


pdf(file=paste0(scoreType,".pdf"),height=2.6,width=8)
row.names(outTab)=outTab[,1]
outTab=outTab[,-1]
corrplot(corr=as.matrix(t(outTab)),
		 title=paste0("\n\n\n\n",scoreType),
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()