## 安装R
```
wget https://sourceforge.net/projects/pcre/files/pcre/8.45/pcre-8.45.tar.gz/download
tar zxf download && cd pcre-8.45
./configure --prefix=/usr/local/pcre-7.8 --libdir=/usr/local/lib/pcre --includedir=/usr/local/include/pcre
make && make install

wget https://curl.se/download/curl-7.55.0.tar.xz
tar xf curl-7.55.0.tar.xz
cd curl-7.55.0/
./configure
make && make install

wget http://mirrors.ctan.org/fonts/inconsolata.zip
unzip inconsolata.zip 
cp -Rfp inconsolata/* /usr/share/texmf
mktexlsr

yum-builddep R
yum install readline readline-devel readline-static libX11-devel libXt-devel libcurl-devel 
yum install libpng libpng-devel libtiff libtiff-devel libjpeg-turbo libjpeg-turbo-devel
wget https://mirror.nju.edu.cn/CRAN/src/base/R-4/R-4.2.3.tar.gz
tar xf R-4.2.3.tar.gz && cd R-4.2.3
./configure --enable-R-shlib --with-x --with-cairo --with-libpng --with-jpeglib --prefix=/DIR/
make -j12 && make install

wget https://ftp.gnu.org/gnu/texinfo/texinfo-6.8.tar.gz
tar zxf texinfo-6.8.tar.gz
cd texinfo-6.8/
./configure --prefix=/data/home/ruanlab/huangxingyu/Tools/texinfo-6.8
make && make install
# add path

# texlive
wget https://mirrors.tuna.tsinghua.edu.cn/CTAN/systems/texlive/tlnet/install-tl-unx.tar.gz
tar zxf install-tl-unx.tar.gz
cd install-tl-20240602/

# 自定义安装路径
export TEXLIVE_INSTALL_PREFIX=~/.local/texlive
export TEXLIVE_INSTALL_TEXDIR=~/.local/texlive/2019

./install-tl


关闭selinux
wget https://download2.rstudio.org/server/rhel8/x86_64/rstudio-server-rhel-2024.04.1-748-x86_64.rpm
yum install rstudio-server-rhel-2024.04.1-748-x86_64.rpm
# vim /etc/rstudio/rserver.conf
rstudio-server start
rstudio-server status


yum install perl-IPC-Cmd perl-CPAN
wget https://www.openssl.org/source/openssl-1.1.1u.tar.gz
tar zxf openssl-1.1.1u.tar.gz
cd openssl-1.1.1u
./config  --prefix=/usr/local/openssl
make
make install
ln -s /usr/local/openssl/lib64/libssl.so.1.1 /usr/lib/libssl.so.1.1
ln -s /usr/local/openssl/lib64/libcrypto.so.1.1 /usr/lib/libcrypto.so.1.1
```

## 安装包
```
install.packages("pak")
library(pak)
repo_get() # 查看镜像
pak::pak("tibble") # 安装包
pak::pak(c("BiocNeighbors", "ComplexHeatmap", "circlize", "NMF")) # 安装包
pak::pak("tibble", lib = "PATH") # 安装包到指定路径
pak::pak("tidyverse/tibble") # 从GitHub安装
pak::pak("tidyverse/tibble@tag") # 从GitHub安装
pak::pak("url::https://cran.r-project.org/src/contrib/Archive/tibble/tibble_3.1.7.tar.gz") # 在线安装
pak::local_install("CytoTRACE") # 本地安装
pak::pak("local::./CytoTRACE_0.3.3.tar.gz") # 本地安装
pak::pak("tibble", upgrade = TRUE) # 更新包的所有依赖，默认不更新依赖
pkg_remove("tibble") # 卸载
pak::pkg_deps("tibble") # 查看包的依赖
pak::pkg_deps_tree("tibble") # 查看包的依赖
pak::pkg_deps_tree("tidyverse/tibble") # 查看github包依赖

# 网络不是太好，可以使用国内的cran、Bioconductor镜像。操作为：在 $HOME/.Rprofile 中添加如下两行
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

R CMD INSTALL Matrix_1.6-5.tar.gz


options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))  # 更换默认镜像
install.packages('Seurat')
```
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

## 分组计数
```
library(dplyr)
mtcars %>% group_by(cyl,am) %>% summarise(total=n())
```

## 将每行合并为字符串
```
df$combined <- apply(df, 1, function(row) paste(row, collapse=""))
```

## 合并图/拼图
```
library(patchwork)
```

## 字符串分割/分隔
```
library(stringr)

df <- data.frame(player=c('John_Wall', 'Dirk_Nowitzki', 'Steve_Nash'), dots=c(22, 29, 18), assists=c(8, 4, 15))
df[c('First', 'Last')] <- str_split_fixed(df$player, '_', 2)
```

## 颜色读取，转换
```
is_color <- function(color_str) {
  # 判断是否为颜色字符串
  # black
  color_names <- colors()
  if (tolower(color_str) %in% tolower(color_names)) {
    return(TRUE)
  }
  
  # 判断是否为 16 进制颜色代码
  # #234567
  hex_pattern <- "^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$"
  if (grepl(hex_pattern, color_str)) {
    return(TRUE)
  }
  # 234567
  hex_pattern <- "^([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$"
  if (grepl(hex_pattern, color_str)) {
    return(TRUE)
  }
  
  # 判断是否为 RGB 字符串
  # rgb(1,2,3) rgb(125,242,32)
  rgb_pattern <- "^rgb\\(\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*\\)$"
  if (grepl(rgb_pattern, color_str, perl = TRUE)) {
    return(TRUE)
  }
  # 125,242,32
  rgb_pattern <- "^\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*$"
  if (grepl(rgb_pattern, color_str, perl = TRUE)) {
    return(TRUE)
  }
  
  return(FALSE)
}

colorts <- function(color_str) {
  # 判断是否为颜色字符串
  # black
  color_names <- colors()
  if (tolower(color_str) %in% tolower(color_names)) {
    # library(gplots)
    # color_str = col2hex(color_str)
    return(color_str)
  }
  
  # 判断是否为 16 进制颜色代码
  # #234567
  hex_pattern <- "^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$"
  if (grepl(hex_pattern, color_str)) {
    return(color_str)
  }
  # 234567
  hex_pattern <- "^([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$"
  if (grepl(hex_pattern, color_str)) {
    return(paste0("#",color_str))
  }
  
  # 判断是否为 RGB 字符串
  # rgb(1,2,3) rgb(125,242,32)
  rgb_pattern <- "^rgb\\(\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*\\)$"
  if (grepl(rgb_pattern, color_str, perl = TRUE)) {
    # color_nums <- as.numeric(gsub("rgb\\(|\\)", "", color_str, perl=TRUE, fixed=FALSE))
    color_nums <- as.numeric(strsplit(gsub("rgb\\(|\\)", "", color_str), ",")[[1]])
    return (rgb(color_nums[1], color_nums[2], color_nums[3], maxColorValue = 255))
  }
  # 125,242,32
  rgb_pattern <- "^\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*$"
  if (grepl(rgb_pattern, color_str, perl = TRUE)) {
    color_nums <- as.numeric(strsplit(color_str, ",")[[1]])
    return (rgb(color_nums[1], color_nums[2], color_nums[3], maxColorValue = 255))
  }
}

color_str <- "125,242,32"
is_color(color_str)
colorts(color_str)
```

## 累计曲线
```
library(ggplot2)
library(dplyr)

# 创建示例数据框
data <- data.frame(
  x = c(1, 2, 3, 4, 5),
  y = c(1, 3, 2, 5, 4)
)

# 计算累计数量
data <- data %>%
  arrange(x) %>%
  mutate(cumulative_count = cumsum(y))

# 绘制累计数量的散点图
p <- ggplot(data, aes(x = x, y = cumulative_count)) +
  geom_point() +
  geom_line() +  # 可选：连接点
  labs(title = "累计数量的散点图",
       x = "X 轴",
       y = "累计数量") +
  theme_minimal()

# 显示图形
print(p)
```