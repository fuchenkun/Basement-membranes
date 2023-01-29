library(ggplot2)       
setwd("C:\\Users\\HASEE\\Desktop\\BMs\\04.PCA")   

bioPCA=function(inputFile=null, outFile=null){
  #读取输入文件,提取数据
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(rt)
	Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
	rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))
	
	#PCA分析
	data.pca=prcomp(data)
	pcaPredict=predict(data.pca)
	PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)
	
	#绘制PCA图
	pdf(file=outFile, height=5, width=6)       
	p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
			 scale_colour_manual(name="",  values=c("blue", "red"))+
		     theme_bw()+
		     theme(plot.margin=unit(rep(1.5,4),'lines'))+
		     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
}

#批次矫正前的PCA图
bioPCA(inputFile="merge.preNorm.txt", outFile="PCA.preNorm.pdf")
#批次矫正后的PCA图
bioPCA(inputFile="merge.normalzie.txt", outFile="PCA.normalzie.pdf")

