library(limma)
library(ggplot2)
library(reshape2)

expFile="panGeneExp.txt"                        
subtypeFile="Subtype_Immune_Model_Based.txt"    
setwd("D:\\biowolf\\panGenes\\13.immuneType")     


exp=read.table(expFile, header=T,sep="\t",row.names=1,check.names=F)
exp=exp[(exp[,"Type"]=="Tumor"),]
exp=as.matrix(exp[,1:(ncol(exp)-2)])
row.names(exp)=gsub(".$","",row.names(exp))
exp=avereps(exp)


subtype=read.table(subtypeFile, header=T,sep="\t",row.names=1,check.names=F)


sameSample=intersect(row.names(subtype),row.names(exp))
subtype=subtype[sameSample,]
subtype=gsub(".+Immune |\\)","",subtype)
exp=exp[sameSample,]
exp=cbind(as.data.frame(exp),subtype)


outTab=data.frame()
geneSig=c()
for(gene in colnames(exp)[1:(ncol(exp)-1)]){
	rt1=exp[,c(gene,"subtype")]
	colnames(rt1)=c("expression","subtype")
	ks=kruskal.test(expression ~ subtype, data = rt1)
	p=ks$p.value
	outTab=rbind(outTab,cbind(gene,pvalue=p))
	Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
	geneSig=c(geneSig,Sig)
}
geneSig=c(geneSig,"")
colnames(exp)=paste0(colnames(exp),geneSig)
write.table(outTab,file="immuneType.xls",sep="\t",row.names=F,quote=F) 


data=melt(exp)
colnames(data)=c("Subtype","Gene","Expression")


p1=ggplot(data,aes(x=Subtype,
                y=Expression,
                fill=Subtype))+
    guides(fill=guide_legend(title="Immune Subtype"))+
    labs(x = "Immune Subtype", y = "Gene expression")+
	geom_boxplot()+ facet_wrap(~Gene,nrow =1)+ theme_bw()
pdf(file="immuneType.pdf",width=9,height=5)
print(p1)
dev.off()
