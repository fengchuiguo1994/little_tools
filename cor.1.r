# https://mp.weixin.qq.com/s?__biz=MzkzMjQyOTA5Nw==&mid=2247484573&idx=1&sn=76df0bf03a9cdf82bc130454b215081f&chksm=c25aa8dbf52d21cdb7a329d0e233febb56b34fba3e4cdd0f32bc43e27dd3ddfe82dd89660eab&cur_album_id=3181722267144339456&scene=190#rd
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggExtra)

p <- ggplot(iris,aes_string(x = 'Sepal.Length',y = 'Sepal.Width')) +
  geom_point(size = 2,color = '#EC0101',alpha = 0.5) +
  theme_bw() +
  # 主题细节调整
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
  ) +
  # 添加回归线
  geom_smooth(method = 'lm',se = T,color = '#F9B208',size = 1.5,fill = '#FEA82F') +
  # 添加相关性系数及p值
  stat_cor(method = "pearson",digits = 3,size=6)

# 添加边际柱形密度图
p1 = ggMarginal(p,type = "densigram",
           xparams = list(binwidth = 0.1, fill = "#B3E283",size = .7),
           yparams = list(binwidth = 0.1, fill = "#8AB6D6",size = .7))
print(p1)




library(ggpmisc)
library(ggpubr)
data("iris")
df <- iris
head(df)
p2 = ggscatter(df, 
          x = "Sepal.Length",
          y = "Petal.Length", 
          palette = "jco", 
          add = "reg.line",                    
          conf.int = TRUE, 
          color = "Species"
) +
  labs(title = "test")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 15),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
  )+
  stat_cor(aes(color = Species), label.x = 3) +
  #添加回归方程
  stat_poly_eq(aes(color = Species,label =  paste(..rr.label..,..p.value.label.., sep = "~~~~")))
print(p2)




# 加载包
library(ggplot2)
library(ggpubr)
# 加载iris数据集
data(iris)
aa<-iris
### 出图1
p3 = ggplot(aa,aes(x=aa$Sepal.Length,y=aa$Petal.Width))+#x,y参数
  geom_point(size=4,aes(color=aa$Species))+#点大小
  scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+#颜色
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.02)+
  labs(y="Sepal Width",x="Sepa Length")+#坐标轴
  theme(axis.title = element_text(color='black',size=9),#主题
    axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
    axis.line = element_line(colour = "black"), 
    axis.title.x=element_text(colour='black', size=11),
    axis.title.y=element_text(colour='black', size=11),
    axis.text=element_text(colour='black',size=9),
    legend.title=element_blank(),
    legend.text=element_text(size=9),
    legend.key=element_blank(),
    legend.background = element_rect(colour = "White"))
print(p3)

### 出图2 

p4 = ggplot(aa,aes(x=aa$Sepal.Length,y=aa$Petal.Length))+
  geom_point(size=4,aes(color=aa$Species))+
  scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.02)+
  theme_classic()+
  labs(y="Sepal Width",x="Sepa Length")+
  theme(axis.title = element_text(color='black',size=9),
    axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
    axis.line = element_line(colour = "black"), 
    axis.title.x=element_text(colour='black', size=11),
    axis.title.y=element_text(colour='black', size=11),
    axis.text=element_text(colour='black',size=9),
    legend.title=element_blank(),
    legend.text=element_text(size=9),
    legend.key=element_blank(),
    legend.background = element_rect(colour = "White"))
print(p4)




# install.packages("devtools")
# # devtools::install_github("Hy4m/linkET", force = TRUE)
# packageVersion("linkET")

library(linkET)
## matrix_data
matrix_data(list(mtcars = mtcars))
## md_tbl
matrix_data(list(mtcars = mtcars)) %>% 
  as_md_tbl()
as_matrix_data(mtcars)
as_md_tbl(mtcars)
correlate(mtcars) %>% 
  as_matrix_data()
correlate(mtcars) %>% 
  as_md_tbl()

correlate(mtcars) %>% 
  as_md_tbl() %>% 
  qcorrplot() +
  geom_square()

library(vegan)
#> 载入需要的程辑包：permute
#> 载入需要的程辑包：lattice
#> This is vegan 2.6-4
data("varespec")
data("varechem")
correlate(varespec[1:30], varechem) %>% 
  qcorrplot() +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

qcorrplot(varespec[1:30], type = "lower") +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))
#> The input data is not a correlation matrix,
#> you can override this behavior by setting the `is_corr` parameter.

## you can set your style
set_corrplot_style()
qcorrplot(mtcars) + geom_square()
#> The input data is not a correlation matrix,
#> you can override this behavior by setting the `is_corr` parameter.


library(dplyr)
#> 
#> 载入程辑包：'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
data("varechem", package = "vegan")
data("varespec", package = "vegan")

mantel <- mantel_test(varespec, varechem,
                      spec_select = list(Spec01 = 1:7,
                                         Spec02 = 8:18,
                                         Spec03 = 19:37,
                                         Spec04 = 38:44)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#> `mantel_test()` using 'bray' dist method for 'spec'.
#> `mantel_test()` using 'euclidean' dist method for 'env'.

qcorrplot(correlate(varechem), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

qpairs(iris) + geom_pairs()





library(ggplot2)
library(RColorBrewer)
library(ggsci)

df=iris
##Spearman's correlation
colnames(df)
df$Sample3<-"Spearman's correlation"

p5 = ggplot(df, aes(x=Sepal.Length, y=Sepal.Width,color=as.factor(Species)))+
  geom_point(df,mapping =  aes(x=Sepal.Length, y=Sepal.Width),size=2)+
  geom_smooth(aes(x=Sepal.Length, y=Sepal.Width),method = 'lm',level=0.95,se = F, size=1.6)+
  stat_cor(method = "spearman",label.x.npc ="left",label.y.npc = 0.97)+#或pearson
  theme_bw()+scale_color_npg()+scale_fill_npg()+
  theme(axis.text=element_text(colour='black',size=9))+
  labs(x="Sepal.Length", y="Sepal.Width", color = "",fill = "")+
  facet_grid( ~Sample3, drop=TRUE,scale="free", space="free_x")+
  scale_y_continuous(expand = c(0, 0), limit = c(2, 4.5))+
  scale_x_continuous(expand = c(0, 0), limit = c(4, 8))+
  annotate("rect", xmin = -1.5, xmax =3,  ymin = -1.5, ymax =3, alpha = 0.1,fill="#FAE3AD") +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(strip.background = element_rect(fill=c("#FAE3AD")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "#B8733A"))+
  theme(axis.text=element_text(colour='black',size=11))
print(p5)
##Pearson's correlation

#将pH转化为数值


df$Sample3<-"Pearson's correlation"

p6 = ggplot(df, aes(x=Sepal.Length, y=Sepal.Width,color=as.factor(Species)))+
  geom_point(df,mapping =  aes(x=Sepal.Length, y=Sepal.Width),size=2)+
  geom_smooth(aes(x=Sepal.Length, y=Sepal.Width),method = 'lm',level=0.95,se = F, size=1.6)+
  stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.97)+#或pearson
  theme_bw()+scale_color_npg()+scale_fill_npg()+
  theme(axis.text=element_text(colour='black',size=9))+
  labs(x="Sepal.Length", y="Sepal.Width", color = "",fill = "")+
  facet_grid( ~Sample3, drop=TRUE,scale="free", space="free_x")+
  scale_y_continuous(expand = c(0, 0), limit = c(2, 4.5))+
  scale_x_continuous(expand = c(0, 0), limit = c(4, 8))+
  annotate("rect", xmin = -1.5, xmax =3,  ymin = -1.5, ymax =3, alpha = 0.1,fill="#CCE3FA") +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(strip.background = element_rect(fill=c("#CCE3FA")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "#328BBF"))+
  theme(axis.text=element_text(colour='black',size=11))
print(p6)





library(ggplot2)
data=iris
colnames(data)
data$Species
pp1 <- ggplot(data)+
  geom_point(aes(Sepal.Length, Sepal.Width, color = Species))+
  scale_color_manual(values = c("setosa" = "#ffbe5d",
                                "versicolor" = "#82a0c3",
                                "virginica" = "#f28c88"))+
  theme_classic()+
  theme(legend.position = c(0.99,0.99),
        legend.justification = c(1,1))+
  guides(colour = guide_legend(""))


# 上方直方图：对应X轴
pp2 <- ggplot(data)+
  geom_histogram(aes(Sepal.Length), 
                 binwidth = 0.2, 
                 fill = "#9db1c1",
                 color = "#ffffff",
                 size = 1)+
  theme_classic()+
  xlab("")+
  ylab("# of vOTUs")


# 右侧直方图：对应Y轴
pp3 <- ggplot(data)+
  geom_histogram(aes(Sepal.Width), 
                 binwidth = 0.2, 
                 fill = "#9db1c1",
                 color = "#ffffff",
                 size = 1)+
  theme_classic()+
  xlab("")+
  ylab("# of vOTUs")+
  coord_flip()

# 拼图：
library(cowplot)

p7 = ggdraw() + 
  draw_plot_label("a", size = 20)+
  # 四个数值中前两个代表左下角位置：
  # 后两个代表图形宽和高：
  draw_plot(pp1, 0.05, 0, .75, .8)+
  draw_plot(pp2, 0.05, 0.8, .75, .2)+
  draw_plot(pp3, 0.8, 0, .2, .8)
print(p7)




library(ggplot2)
library(scales)
p8 = ggplot(iris, aes())+
  geom_jitter(aes(x=iris$Sepal.Length,y=iris$Sepal.Width),position=position_jitter(0.17),size=3, alpha=1,color="#925F8C")+
  labs(x="Sepal.Length", y="Sepal.Width")+
  geom_smooth(aes(x=iris$Sepal.Length, y=iris$Sepal.Width),method=lm,level=0.95,color="#925F8C", formula = y~poly(x, 2))+#拟合线
  geom_jitter(aes(x=iris$Sepal.Length,y=iris$Petal.Length),position=position_jitter(0.17),size=3, alpha=1,color="#537F88", formula = y~poly(x, 2))+
  geom_smooth(aes(x=iris$Sepal.Length,y=iris$Petal.Length),method=lm,level=0.95,color="#537F88")+
  theme_classic()+
  scale_y_continuous(sec.axis = sec_axis(~.*0.5, name = 'SOC (g/kg)'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=9))
print(p8)




library(ggplot2)
library(scales)
colnames(iris)
p9 = ggplot(iris,aes(x = Sepal.Length,y = Sepal.Width))+
  geom_point(size = 3)+
  geom_smooth(method = "lm",color = "#9f0000")+
  stat_cor(method = "pearson",#相关系数计算方法
           label.sep = '\n',#r和p换行
           p.accuracy = 0.0001,#P值精确度
           r.digits = 4,#相关性系数小数点位数
           label.x = 7,size = 7.5) +
  scale_x_continuous(limits = c(4,8.5),expand = c(0,0)) +
  scale_y_continuous(limits = c(1.5,5),expand = c(0,0)) +
  labs(x = "MAPK-PROGENy GES",y = "IL-17 GES") +
  ggtitle("PROGENy versus IL-17") +
  theme_classic() +
  theme(axis.text = element_text(color = "black",size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5))
print(p9)

print(p1+p2+p3+p4+p5+p6+p7+p8+p9)