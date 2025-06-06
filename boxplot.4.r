library(ggpubr)
library(reshape2)
mydat = read.table("boxplot.4.r.result",header=T)
md = melt(mydat,id=c("gene"),measure=c("promoterpearson","promoterspearman","X100Kpearson","X100Kspearman"))

meandata = aggregate(value ~ variable, data = md, mean)

pdf("k562.cor.box.pdf")
# ggboxplot(md,x="variable",fill="variable",y="value",color = "black",width=0.6,outlier.colour=NA)+
ggboxplot(md,x="variable",fill="variable",y="value",color = "black",width=0.6,outlier.colour=NA)+
  geom_hline(yintercept = median(md$value), linetype=2) +
  geom_point(data=meandata,mapping=aes(x=variable,y=value),size=5) +
  theme_bw() +
  theme(panel.grid=element_blank(),axis.line=element_line(size=1,colour="black"))

ggviolin(md,x="variable",fill="variable",y="value",color = "black",width=0.6,outlier.colour=NA)+
  geom_hline(yintercept = median(md$value), linetype=2) +
  geom_point(data=meandata,mapping=aes(x=variable,y=value),size=5) +
  theme_bw() +
  theme(panel.grid=element_blank(),axis.line=element_line(size=1,colour="black"))
dev.off()

