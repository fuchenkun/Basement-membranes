library(ggpubr)
library(reshape2)
setwd("C:\\Users\\HASEE\\Desktop\\Example\\8_step8")
data1=read.table("ssgsea_result.txt",header = T,sep = "\t",check.names =F)
data1=t(data1)
write.table(data1,"data1.txt",row.names = T,col.names = F,sep = "\t",quote = F)
data1=read.table("data1.txt",header = T,sep = "\t",check.names =F)
data2=read.table("clinical.txt",header = T,sep = "\t",check.names = F)
data3=merge(data2,data1,by="id")

imc=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
      "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
      "Tfh","Th1_cells","Th2_cells","TIL","Treg")
imcd=data3[,c(imc,"mygroup")]
newdata=melt(imcd,id.vars=c("mygroup"))
colnames(newdata)=c("Group","Type","Score")
newdata$Group=factor(newdata$Group, levels=c("Con","IPF"))
a=ggboxplot(newdata, x="Type", y="Score", color = "Group",
            ylab="Score",add = "none",xlab="",palette = c("green","red") )
a=a+rotate_x_text(51)
pdf("imc.pdf",8,8)            
a+stat_compare_means(aes(group=newdata$Group),symnum.args=list(cutpoints = c(0,0.001,0.01,0.05, 1), symbols = c("***", "**","*", "ns")),label = "p.signif")
dev.off()


imf=c("APC_co_inhibition","APC_co_stimulation","CCR",
      "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
      "MHC_class_I","Parainflammation","T_cell_co-inhibition",
      "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
imfd=data3[,c(imf,"mygroup")]
newdata=melt(imfd,id.vars=c("mygroup"))
colnames(newdata)=c("Group","Type","Score")
newdata$Group=factor(newdata$Group, levels=c("Con","IPF"))
a=ggboxplot(newdata, x="Type", y="Score", color = "Group",
            ylab="Score",add = "none",xlab="",palette = c("blue","red") )
a=a+rotate_x_text(51)
pdf("imf.pdf",8,8)            
a+stat_compare_means(aes(group=newdata$Group),symnum.args=list(cutpoints = c(0,0.001,0.01,0.05, 1), symbols = c("***", "**","*", "ns")),label = "p.signif")
dev.off()
