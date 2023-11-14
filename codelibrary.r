## 命令行参数
args=commandArgs(T) 
# args[1] args[2] args[3]...

## scale
rowscale = function(x) {
  return = (x - mean(x))/sd(x) * 1000
}
mydatscale = apply(mydat, 2, rowscale)

## 标准化
```
scalemat = scale(mat, center = TRUE, scale = TRUE) # 默认按列
```

## 安装 rstudio server 服务
```
yum -y install epel-release
yum install R
wget https://download2.rstudio.org/server/centos7/x86_64/rstudio-server-rhel-2023.09.1-494-x86_64.rpm
yum install rstudio-server-rhel-2023.09.1-494-x86_64.rpm
adduser rstudio
passwd rstudio # rstudio123456
usermod -g rstudio-server rstudio
# firewall-cmd --zone=public --list-ports
firewall-cmd --zone=public --add-port=8787/tcp --permanent
firewall-cmd --reload

# ip:8787/
```

## 画布分屏
```
aa = par(cex=0.7, mai=c(0.1,0.1,0.1,0.1))
bb = par(fig=c(0,0.7,0,0.7))
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)
plot(dose, drugA, type="b")

cc = par(cex=0.7, mai=c(0.1,0.1,0.1,0.1))
dd = par(fig=c(0.7,1,0.7,1))
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)
plot(dose, drugA, type="b")

dev.off()

bb
boxplot(dose, drugA)


par(no.readonly = FALSE)
par(cex=0.7, mai=c(0.1,0.1,0.1,0.1), no.readonly = FALSE)
par(fig=c(0,0.7,0,0.7))
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)
aa <- plot(dose, drugA, type="b")

par(fig=c(0,0.7,0.7,1), new=TRUE)
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)
gg = plot(dose, drugA, type="b")

par(fig=c(0,0.7,0,0.7), new=TRUE)
bb = boxplot(dose, drugA)

par(fig=c(0,0.7,0.7,1))
boxplot(dose, drugA)

par(fig=c(0,0.7,0,0.7))


par(cex=0.7, mai=c(0.1,0.1,0.1,0.1), no.readonly = FALSE)
par(fig=c(0,0.7,0,0.7))
plot(dose, drugA, type="b")
par(fig=c(0,0.7,0.7,1), new=TRUE)
plot(dose, drugA, type="b")


par(fig=c(0,0.7,0.7,1), new=TRUE)
boxplot(dose, drugA)
```

## 热图
```
library(factoextra)
library(pheatmap)
set.seed(1234)
mat = matrix(rnorm(6*36), nrow = 6, ncol = 36)
rownames(mat) = seq(1,6)
colnames(mat) = seq(1,36)

scalemat = scale(mat, center = TRUE, scale = TRUE) # 默认按列

p1=pheatmap(matrix,scale = "row",
            display_numbers = T,
            cluster_rows = F,
            cluster_cols = F)
p2=pheatmap(matrix,scale = "column",
            display_numbers = T,
            cluster_rows = F,
            cluster_cols = F)
p3 = pheatmap(mat,scale = "column",
            display_numbers = T,
            cluster_rows = F,
            cluster_cols = T,
            kmeans_k = 4)

cl <- kmeans(mat, 3, nstart = 24)

cl$size
cc = cl$cluster

cc1 = cc[cc == 1]
fltdat1 = fltdat[,names(cc1)]
fltdat1$mean = rowMeans(fltdat1)

cc2 = cc[cc == 2]
fltdat2 = fltdat[,names(cc2)]
fltdat2$mean = rowMeans(fltdat2)

cc3 = cc[cc == 3]
fltdat3 = fltdat[,names(cc3)]
fltdat3$mean = rowMeans(fltdat3)

cc4 = cc[cc == 4]
fltdat4 = fltdat[,names(cc4)]
fltdat4$mean = rowMeans(fltdat4)

m1 = data.frame(name=rownames(fltdat1),mm=fltdat1$mean,type="m1")
m2 = data.frame(name=rownames(fltdat2),mm=fltdat2$mean,type="m2")
m3 = data.frame(name=rownames(fltdat3),mm=fltdat3$mean,type="m3")
m4 = data.frame(name=rownames(fltdat4),mm=fltdat4$mean,type="m4")
alldat = rbind(m1,m2,m3,m4)
alldat$name = factor(alldat$name,levels = unique(alldat$name))
library(ggplot2)
ggplot(data = alldat, mapping = aes(x = name, y = mm, colour = type, group = type)) +
  geom_line(size=1) +
  ylim(c(0,1))
```
```
library(RColorBrewer)
library(dplyr)
aa = c("A","B","C","D","E","F","G","H")
bb = c(3,2,3,4,2,6,7,8)
cc = data.frame(aa=aa,bb=bb)
cc$index = as.numeric(rownames(cc))
sp=spline(cc$index,cc$bb,n=10)
cc$xleft = -(ceiling(max(bb)/2))
cc$xright = 0
cc$ybottom = seq(0,nrow(cc)-1)
cc$ytop = seq(1,nrow(cc))
col = colorRampPalette(brewer.pal(8, "Blues"))(8)
colframe = data.frame(bb=seq(0,length(col)),col=c('white',col))
final = left_join(x=cc, y=colframe, by="bb")

col = "red"
par(fig=c(0,0.7,0,1))
# plot(c(1,1),type='n')
# barplot(cc$bb,col=col,horiz=T,border=col,width=1,space=0,names.arg="cell",legend.text="GEX",main="geneGEX",sub="GEXbarplot",xlab="UMI",ylab="cell")
barplot(cc$bb,col=col,horiz=T,border=col,width=1,space=0,names.arg="cell",legend.text="GEX",main="geneGEX",sub="GEXbarplot",xlab="UMI")
rect(xleft=cc$xleft, ybottom=cc$ybottom, xright=cc$xright, ytop=cc$ytop,col=final$col,border=final$col,lwd=0.1)
# rect(xleft=c(1,2,3), ybottom=c(0,0,0), xright=c(4,5,6), ytop=c(1,1,1),col=c('red','green','blue'),border=c('blue','green','red'),lwd=2)
# lines(sp,col=”green”,type=”l”,xlim=c(0,6),ylim=c(0,30),lwd=2,xlab=”WEEK”,ylab=”STUDENT”,main=”lesson”,sub=”class”,col.main=”green”,font.main=2)
lines(sp$y,sp$x-0.5,col="green",type="l",lwd=2,xlab="WEEK",ylab="STUDENT",main="lesson",sub="class",col.main="green",font.main=2)
points(y=cc$index-0.5,x=cc$bb)
# points(x=cc$index,y=cc$bb)
axis(2,at=cc$index-0.5,labels=cc$aa)


ff = data.frame(aa=c(1,2,3),bb=c(4,5,6),cc=c(7,8,9))
rownames(ff) 
ff["2","aa"]

df <- scale(mtcars)
heatmap(df, scale = "none")
heatmap(df, Colv = NA, Rowv = NA, scale="column")
```

## 设置色阶
```
library(RColorBrewer)

args=commandArgs(T)

cols<-brewer.pal(9, "YlOrRd") # cols<-brewer.pal(3, "YlOrRd") 
pal<-colorRampPalette(cols)
mycolors<-pal(16)

pdf(args[2])
mydat<-read.table(args[1])
names(mydat)<-c("dis",'pet')
plot(c(2.5,9),c(0,1), type = "n", ylab = "distance",xlab = "",xaxt="n",main = "complex=x")
for(i in 2:15){
    tmp<-mydat[(mydat$pet==i),]
    dd<-density(log10(tmp$dis+1))
    lines(dd, col = mycolors[i], lwd = 2)
}
axis(1,at=c(1,2,3,4,5,6,7,8,9),labels = c("10","100","1k","10k","100k","1M","10M",'100M','1000M'))
tmp<-mydat[(mydat$pet > 15),]
dd<-density(log10(tmp$dis+1))
lines(dd, col = mycolors[16], lwd = 2)
id<-paste("=",2:15,sep="")
id[15]<-'>15'
colo<-mycolors[2:16]
legend("topleft",id,lty=array(1,c(2,2)),col=colo,lwd=2)


plot(c(3,9),c(0,1), type = "n", ylab = "distance",xlab = "",xaxt="n",main = "complex>=x")
for(i in 2:15){
    tmp<-mydat[(mydat$pet>=i),]
    dd<-density(log10(tmp$dis+1))
    lines(dd, col = mycolors[i], lwd = 2)
}
axis(1,at=c(1,2,3,4,5,6,7,8,9),labels = c("10","100","1k","10k","100k","1M","10M",'100M','1000M'))
tmp<-mydat[(mydat$pet > 15),]
dd<-density(log10(tmp$dis+1))
lines(dd, col = mycolors[16], lwd = 2)
id<-paste(">=",2:15,sep="")
id[15]<-'>15'
colo<-mycolors[2:16]
legend("topleft",id,lty=array(1,c(2,2)),col=colo,lwd=2)

dev.off()

```

## 多维数据的复合相关性分析
```
cor_taylor <- function(X){
    tryCatch(
        { n <- ncol(X)
          return((1/sqrt(n))*sd(eigen(cor(X))$values)) },
        warning = function(w) {return(4)},
        error = function(e) { return(5)}
    )
}

cor_taylorspearman <- function(X){
    tryCatch(
        { n <- ncol(X)
          return((1/sqrt(n))*sd(eigen(cor(X,method="spearman"))$values)) },
        warning = function(w) {return(4)},
        error = function(e) { return(5)}
    )
}

```