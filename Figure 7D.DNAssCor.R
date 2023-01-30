library(limma)
library(corrplot)

expFile="panGeneExp.txt"                            
scoreFile="StemnessScores_DNAmeth_20170210.tsv"     
scoreType="DNAss"                                  
setwd("C:\Users\HASEE\Desktop\BMs\39.DNAssCor")     


exp=read.table(expFile, header=T,sep="\t",row.names=1,check.names=F)
exp=exp[(exp[,"Type"]=="Tumor"),]


STEM=read.table(scoreFile, header=T,sep="\t",row.names=1,check.names=F)
STEM=t(STEM)


outTab=data.frame()
pTab=data.frame()

for(i in levels(factor(exp[,"CancerType"]))){
    exp1=exp[(exp[,"CancerType"]==i),]
    exp1=as.matrix(exp1[,1:(ncol(exp1)-2)])
    row.names(exp1)=gsub(".$","",row.names(exp1))
    exp1=avereps(exp1)

	
	sameSample=intersect(row.names(STEM),row.names(exp1))
	STEM1=STEM[sameSample,]
	exp1=exp1[sameSample,]
	
    x=as.numeric(STEM1[,scoreType])
    pVector=data.frame(i)
    outVector=data.frame(i)
	
	genes=colnames(exp1)
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

colNames=c("CancerType",colnames(exp)[1:(ncol(exp)-2)])
colnames(outTab)=colNames
write.table(outTab,file="DNAssCor.cor.txt",sep="\t",row.names=F,quote=F)

colnames(pTab)=colNames
write.table(pTab,file="DNAssCor.pval.txt",sep="\t",row.names=F,quote=F)


pdf("DNAssCor.pdf",height=2.6,width=8)
row.names(outTab)=outTab[,1]
outTab=outTab[,-1]
corrplot(corr=as.matrix(t(outTab)),
		 title=paste0("\n\n\n\n",scoreType),
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()


