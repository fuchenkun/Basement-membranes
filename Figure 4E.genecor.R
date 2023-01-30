install.packages("psych")
install.packages("ggcorrplot")

library(psych)
library(ggcorrplot)

setwd("C:\\Users\\HASEE\\Desktop\\new\\5.Nontumor_BMs_gene\\9_step9")
data1=read.table("ssgsea_result.txt",header = T,sep = "\t",check.names =F)
data2=read.table("hub.txt",header = F,sep = "\t",check.names =F)
data3=read.table("BMs_gene1.txt",header = T,sep = "\t",check.names =F,row.names = 1)
data3=data3[as.vector(data2[,1]),]
data3=t(data3)
data1=t(data1)
colnames(data1)=data1[1,]
data1=data1[-1,]
data1=apply(data1,2,as.numeric)
data3=apply(data3,2,as.numeric)


mycor<- corr.test(data3,data1,method = 'spearman')
rvalue <- mycor$r
pvalue <- mycor$p.adj
rnewdf=data.frame(rvalue)
write.table(rnewdf,"cor.txt",sep = "\t",quote = F)
pnewdf=data.frame(pvalue)
write.table(pnewdf,"p.txt",sep = "\t",quote = F)

pdf("geneimc.pdf",10,8)
ggcorrplot(t(rvalue), show.legend = T, p.mat = t(mycor$p.adj), digits = 2,  
          sig.level = 0.05,
           insig = 'blank',lab_size = 3,
           lab = T)
dev.off()







