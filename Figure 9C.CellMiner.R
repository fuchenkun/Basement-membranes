library(impute)
library(limma)
library(ggplot2)
library(ggpubr)

drugFile="drug.txt"        
expFile="geneExp.txt"      
geneFile="gene.txt"        
corPlotNum=20             
setwd("C:\Users\HASEE\Desktop\BMs\42.CellMiner")    


rt=read.table(drugFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
drug=rt[,2:ncol(rt)]
dimnames=list(rownames(drug),colnames(drug))
data=matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)


mat=impute.knn(data)
drug=mat$data
drug=avereps(drug)


rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp=avereps(data)


gene=read.table(geneFile,sep="\t",header=F,check.names=F)
genelist=as.vector(gene[,1])
genelist=gsub(" ","",genelist)
genelist=intersect(genelist,row.names(exp))
exp=exp[genelist,]


outTab=data.frame()
for(Gene in row.names(exp)){
	x=as.numeric(exp[Gene,])
	#??ҩ??ѭ??
	for(Drug in row.names(drug)){
		y=as.numeric(drug[Drug,])
		corT=cor.test(x,y,method="pearson")
		cor=corT$estimate
		pvalue=corT$p.value
		if(pvalue<0.05){
			outVector=cbind(Gene,Drug,cor,pvalue)
			outTab=rbind(outTab,outVector)
		}
	}
}

outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="drugCor.txt",sep="\t",row.names=F,quote=F)


plotList=list()
if(nrow(outTab)<corPlotNum){
	corPlotNum=nrow(outTab)
}
for(i in 1:corPlotNum){
	Gene=outTab[i,1]
	Drug=outTab[i,2]
	x=as.numeric(exp[Gene,])
	y=as.numeric(drug[Drug,])
	cor=sprintf("%.03f",as.numeric(outTab[i,3]))
	pvalue=0
	if(as.numeric(outTab[i,4])<0.001){
		pvalue="p<0.001"
	}else{
		pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
	}
	df1=as.data.frame(cbind(x,y))
	p1=ggplot(data = df1, aes(x = x, y = y))+
		geom_point(size=1)+
		stat_smooth(method="lm",se=FALSE, formula=y~x)+
		labs(x="",y="",title = paste0(Gene,", ",Drug),subtitle = paste0("Cor=",cor,", ",pvalue))+
		theme(axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank())+
		theme_bw()
	plotList[[i]]=p1
}


nrow=ceiling(sqrt(corPlotNum))
ncol=ceiling(corPlotNum/nrow)
pdf(file="drugCor.pdf", width = 13,height =9)
ggarrange(plotlist=plotList,nrow=nrow,ncol=ncol)
dev.off()

