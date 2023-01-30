library(ggplot2)          
logFCfilter=1             
adj.P.Val.Filter=0.05    
inputFile="all.txt"        
setwd("C:\\Users\\HASEE\\Desktop\\BMs\\08.volcano")   

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
rt=cbind(rt, Sig=Sig)


p=ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("blue", "gray", "red"))+
    xlim(-5,5)+
    labs(title = " ")+
    geom_vline(xintercept=c(-logFCfilter,logFCfilter), col="black", cex=0.5, linetype=2)+
    geom_hline(yintercept= -log10(adj.P.Val.Filter), col="black", cex=0.5, linetype=2)+
    theme(plot.title=element_text(size=16, hjust=0.5, face="bold"))
p=p+theme_bw()


pdf(file="volcano.pdf", width=6, height=5.1)
print(p)
dev.off()
