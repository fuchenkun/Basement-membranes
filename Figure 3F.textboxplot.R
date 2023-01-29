


#???ð?
library(limma)
library(ggpubr)

expFile="diffGeneExp.txt"        #?????????ļ?
geneFile="intersectGenes.txt"      #?????б??ļ?
setwd("C:\\Users\\HASEE\\Desktop\\BMs\\17.test_boxplot")     #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#??ȡ?????б??ļ?,??ȡ?????????????ı???��
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=F]

#??ȡ??Ʒ???ͣ????ñȽ???
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
Type=ifelse(Type=="con", "Con", "Treat")
my_comparisons=list()
my_comparisons[[1]]=levels(factor(Type))

#?Ի???????ѭ???????Ʋ???????ͼ
for(i in row.names(data)){
	#data[i,][data[i,]>quantile(data[i,], 0.99)]=quantile(data[i,], 0.99)
	rt1=data.frame(expression=data[i,], Type=Type)
	
	#?Բ??????????п??ӻ???????????ͼ
	boxplot=ggboxplot(rt1, x="Type", y="expression", color="Type",
				      xlab="",
				      ylab=paste(i, "expression"),
				      legend.title="",
				      palette = c("blue", "red"),
				      add = "jitter")+ 
		stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
		#stat_compare_means(comparisons = my_comparisons)
			
	#????ͼƬ
	pdf(file=paste0("boxplot.",i,".pdf"), width=5, h4.eit=4.5)
	print(boxplot)
	dev.off()
}
		

###