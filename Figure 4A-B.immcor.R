install.packages("corrplot")

library(corrplot)

setwd("C:\\Users\\HASEE\\Desktop\\new\\5.Nontumor_BMs_gene\\7_step7")
data1 = read.table("ssgsea_result.txt",header = T,sep = "\t",check.names = F,row.names = 1)

data2 = read.table("type.txt",header = T,sep = "\t",check.names = F)
t1 = data2[data2$Type=="imc",]
t2= data2[data2$Type=="imf",]
rownames(t1)=t1$id
rownames(t2)=t2$id
imc=data1[rownames(t1),]
imf=data1[rownames(t2),]
imc=t(imc)
imf=t(imf)

pdf("corrplot_imc.pdf",10,8)              
corrplot(corr=cor(imc),method="color",order = "alphabet",tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()


pdf("corrplot_imf.pdf",10,8)              
corrplot(corr=cor(imf),method="color",order = "alphabet",tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

