library(venn)
diff1File="LASSO.gene.txt" 
diff3File="SVM-RFE.gene.txt"
diff2File="rfGenes.txt" 

setwd("C:\\Users\\HASEE\\Desktop\\BMs\\15.venn")   
geneList=list()


rt=read.table(diff1File, header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
uniqGene=unique(geneNames)
geneList[["LASSO"]]=uniqGene


rt=read.table(diff2File, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
uniqGene=unique(geneNames)
geneList[["RandomForest"]]=uniqGene

rt=read.table(diff3File, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
uniqGene=unique(geneNames)
geneList[["SVM-RFE"]]=uniqGene



mycol=c("#fbb4ae","#b3cde3","#ccebc5")


# mycol=c("#b3cde3","#ccebc5")
pdf(file="venn.pdf", width=5, height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()


intersectGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt", intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)


