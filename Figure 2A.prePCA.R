library(limma)
library(sva)
mergeFile="merge.preNorm.txt"            
normalizeFile="merge.normalzie.txt"         
setwd("C:\\Users\\HASEE\\Desktop\\BMs\\04.PCA")        
files=c("GSE32537.txt", "GSE53845.txt")       


geneList=list()
for(i in 1:length(files)){
	fileName=files[i]
    rt=read.table(fileName, header=T, sep="\t", check.names=F)
    header=unlist(strsplit(fileName, "\\.|\\-"))
    geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)


allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
	fileName=files[i]
    header=unlist(strsplit(fileName, "\\.|\\-"))
    
    rt=read.table(fileName, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))
    
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    rt=normalizeBetweenArrays(rt)
    
    if(i==1){
    	allTab=rt[intersectGenes,]
    }else{
    	allTab=cbind(allTab, rt[intersectGenes,])
    }
    batchType=c(batchType, rep(header[1],ncol(rt)))
}


allTabOut=rbind(geneNames=colnames(allTab), allTab)
write.table(allTabOut, file=mergeFile, sep="\t", quote=F, col.names=F)


normalizeTab=ComBat(allTab, batchType, par.prior=TRUE)
normalizeTab=rbind(geneNames=colnames(normalizeTab), normalizeTab)
write.table(normalizeTab, file=normalizeFile, sep="\t", quote=F, col.names=F)