library(pheatmap)            
inputFile="panGeneExp.txt"   
setwd("C:\\Users\\HASEE\\Desktop\\BMs\\30.pheatmap")              
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)   


logFcTab=data.frame()
pTab=data.frame()
genelist=colnames(rt[,(1:(ncol(rt)-2))])
for(cancerType in levels(factor(rt[,"CancerType"]))){
	data=rt[(rt[,"CancerType"]==cancerType),]
	normal=data[(data[,"Type"]=="Normal"),]
	tumor=data[(data[,"Type"]=="Tumor"),]
	if(nrow(normal)>=5){
		logFcVector=data.frame(cancerType)
		pVector=data.frame(cancerType)
		normal=normal[,(1:(ncol(normal)-2))]
		tumor=tumor[,(1:(ncol(tumor)-2))]
		for(gene in colnames(tumor)){
			logFC=mean(tumor[,gene])-mean(normal[,gene])
			test=wilcox.test(tumor[,gene],normal[,gene])
			pVector=cbind(pVector,test$p.value)
			logFcVector=cbind(logFcVector,logFC)
		}
		logFcTab=rbind(logFcTab,logFcVector)
		pTab=rbind(pTab,pVector)
	}
}

colnames(logFcTab)=c("CancerType",genelist)
write.table(logFcTab,file="logFC.txt",sep="\t",row.names=F,quote=F)

colnames(pTab)=c("CancerType",genelist)
write.table(pTab,file="pvalue.txt",sep="\t",row.names=F,quote=F)


row.names(logFcTab)=logFcTab[,1]
logFcTab=logFcTab[,-1]
pdf(file="heatmap.pdf",height=6,width=5)
pheatmap(logFcTab, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         fontsize = 8,
         fontsize_row=8,
         fontsize_col=8)
dev.off()
