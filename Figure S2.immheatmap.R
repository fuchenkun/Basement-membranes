library(pheatmap)
setwd("C:\\Users\\HASEE\\Desktop\\BMs\\22.heatmap")
data1=read.table("clinical.txt",header = T,sep = "\t")
data2=read.table("ssgsea_result.txt",header = T,sep = "\t",check.names = F)
data3=merge(data1,data2,by="id")
write.table(data3,"data3.txt",sep = "\t",quote = F,row.names = F)

ht=read.table("data3.txt",header = T,sep = "\t")
typedf=read.table("type.txt",header = T,sep = "\t")
ndf=merge(typedf,ht,by="id")
write.table(ndf,"ndf.txt",sep = "\t",row.names = F,quote = F)

ht=read.table("ndf.txt",header = T,sep = "\t",row.names = 1)
Group=c(rep("Con",58),rep("IPF",207))    
names(Group)=colnames(ht)
Group=as.data.frame(Group)
mytype=read.table("type.txt",header = T,sep = "\t",row.names = 1)
pdf("heatmap.pdf",10,8)
pheatmap(ht, annotation_col = Group, annotation_row = mytype,cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("blue", "black", "red"))(50),border=F,show_colnames = F
         )
dev.off()
