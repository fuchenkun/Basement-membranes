inputFile="panGeneExp.txt"   
setwd("C:\\Users\\HASEE\\Desktop\\BMs\\28.boxplot")              
rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)    
rt=rt[(rt[,"Type"]=="Tumor"),]    
data=rt[,(1:(ncol(rt)-2))]         
geneNum=ncol(data)                 


pdf(file="boxplot.pdf",width=6,height=5)    
par(mar=c(5, 5, 2, 2))                      
boxplot(data,ylab ="Gene expression",col = rainbow(geneNum),xaxt = "n",outline = FALSE)    #????boxplot
text(1:geneNum, par("usr")[3]-0.25, srt=45, adj=1, labels=colnames(data), xpd=TRUE, cex=0.8)  #x???̶?
dev.off()

