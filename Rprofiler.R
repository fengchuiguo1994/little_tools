library(ggplot2)
b=read.table("input.txt")
step<-0.025
length<-(ceiling(1/step)-1)*3 #117
num<-c(0:(length-1))
lengthp<-length+1
#b<-read.table(Infile,header=F,sep="\t")
data<-rbind(data.frame(num,meth=as.double(b[1,2:lengthp]),Context=rep("CG",length),sub_genome=rep("small",length)),data.frame(num,meth=as.double(b[2,2:lengthp]),Context=rep("CHG",length),sub_genome=rep("small",length)),data.frame(num,meth=as.double(b[3,2:lengthp]),Context=rep("CHH",length),sub_genome=rep("small",length)),data.frame(num,meth=as.double(b[4,2:lengthp]),Context=rep("CG",length),sub_genome=rep("big",length)),data.frame(num,meth=as.double(b[5,2:lengthp]),Context=rep("CHG",length),sub_genome=rep("big",length)),data.frame(num,meth=as.double(b[6,2:lengthp]),Context=rep("CHH",length),sub_genome=rep("big",length)) )

pdf("output.pdf",width=10,height=6)
max<-max(data$meth)
max<-max*1.1
p <- ggplot(data,aes(x=num,y=meth,color=Context,linetype=sub_genome))
p + geom_line(lwd=1.5)+xlab("")+theme_bw()+theme(panel.grid=element_blank())+scale_x_continuous("",limits=c(0, length-1),breaks=c(0,length/3,length/3*2,length-1),labels=c("up","TSS","TES","down"))+geom_vline(xintercept=c(length/3,length/3*2),linetype="dotted",color="blue",size=1)+
  theme(legend.key=element_rect(linetype='dashed',color="white"),axis.text.y = element_text(size=13),axis.text.x = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=15),legend.key.height=unit(1.2,'cm')) +ylim(0,max)
dev.off()