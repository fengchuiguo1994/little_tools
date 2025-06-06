library(ggtern)
mydat<-read.table("enrich.compare.forsuanyuantu.txt",header=F)
rownames(mydat)<-mydat$V1
colnames(mydat)<-c("gene","A549.CTCF","A549.Pol2","GM12878.CTCF")
pdf("enrich.compare.sanyuantu.pdf",width=12)
ggtern(data=mydat,aes(A549.CTCF,A549.Pol2,GM12878.CTCF))+geom_point(aes(fill=gene,group=gene),color="red",alpha=0.25)+guides(fill=FALSE)+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
dev.off()
