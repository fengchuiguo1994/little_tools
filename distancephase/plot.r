library(ggplot2)
mydat <- read.table("pet.txt",header=F)
names(mydat) <- c("type","count")
mydat$normal <- log10(mydat$count)
pdf("dis.pdf")
ggplot(mydat,aes(x=normal,group=type,colour=factor(type)))+geom_line(stat = 'density')+theme_bw()+theme(axis.title.y=element_text(size=14,face="bold"),axis.text.x=element_text(size=12,face="bold",vjust=1, hjust=1,angle=45),axis.text.y=element_text(size=8),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()


