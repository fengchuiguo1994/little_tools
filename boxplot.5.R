# https://github.com/IndrajeetPatil/ggstatsplot/issues/796
# https://mp.weixin.qq.com/s/ucw5OI6EdhQp-Jf1xV3H2A
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(cowplot)
# Simulating data 'fm'
set.seed(123)
fm <- matrix(rnorm(100), nrow=200) %>% as.data.frame()
colnames(fm) <- c(paste0("Gene_", 1))
rownames(fm) <- paste0("Sample_", 1:200)
fm$alt_subtype <- rep(c('type1','type2','type3','type4'),time=c(50,50,50,50))
fm$group=NULL

p1 = ggplot(fm, aes(x = alt_subtype, y = Gene_1, fill = alt_subtype)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = "grey20",size = 1) +
  ylim(c(min(fm$Gene_1)-(max(fm$Gene_1)-min(fm$Gene_1))*0.1,
    max(fm$Gene_1)+(max(fm$Gene_1)-min(fm$Gene_1))*0.1)) +
  scale_fill_manual(values = brewer.pal(6, "Paired")[c(6,5,4,3)]) +
  ylab("Cytolytic score (log2)") +
  xlab("") +
  ggtitle('Gene_1') +
  theme_cowplot() +
  theme(legend.position="none",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0)) +
  geom_signif(comparisons = list(c(1, 2)),map_signif_level = T, textsize = 5, test = wilcox.test,step_increase = 0.2)+
  geom_signif(comparisons = list(c(3, 4)),map_signif_level = T, textsize = 5, test = wilcox.test,step_increase = 0.2)+
  theme(plot.margin = unit(c(0.5,1,0,1), "cm"))
print(p1)




library(tidyverse)
library(ggpubr)
library(magrittr)
library(ggsci)
library(ggsignif)
text <- fm %>%
  group_by(alt_subtype) %>%
  top_n(1,Gene_1) %>%
  column_to_rownames('alt_subtype') %>%
  set_colnames('num') %>%
  mutate(num=round(num,digits = 1)) %>%
  mutate(alt_subtype=rownames(.))  

p2 = ggplot(fm, aes(alt_subtype, Gene_1)) +  
  stat_boxplot(geom="errorbar", position=position_dodge(width=0.2), width=0.1) + # 使用 stat_boxplot 添加误差条
  geom_boxplot(position=position_dodge(width=0.2), width=0.4) + # 添加箱线图
  geom_point(aes(fill=alt_subtype, group=alt_subtype, color=alt_subtype), pch=21, position=position_dodge(0.2)) +
  geom_text(data=text, aes(label=num, y=num + 0.1), size=4, color="black", hjust=0.5, vjust=0.5) + # 添加文本标签
  geom_signif(comparisons=list(c("type1", "type2")), map_signif_level=TRUE, textsize=6, test="wilcox.test", step_increase=0.2) + # 为特定比较添加显著性标签
  scale_size_continuous(range=c(1, 3)) + # 调整大小比例
  stat_cor(label.y=4, aes(label=paste(..rr.label.., ..p.label.., sep="~`,`~"), group=1), color="black", label.x.npc="middle") + # 添加相关性标签
  stat_regline_equation(label.y=3.75, aes(group=1), color="black", label.x.npc="middle") +  scale_fill_simpsons(alpha=0.7) +
  scale_color_simpsons(alpha=0.7) +  labs(x=NULL, y=NULL) + # 删除x和y轴标签
  theme(    plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), units="cm"),  # 调整图形边距
    axis.line=element_line(color="black", size=0.4), # 自定义坐标轴线    
    panel.grid.minor=element_blank(),  # 删除次要网格线
    panel.grid.major=element_line(size=0.2, color="#e5e5e5"), # 自定义主要网格线
    panel.background=element_blank(),  # 删除面板背景
    axis.text.y=element_text(color="black", size=10, face="bold"), # 自定义y轴上的文本
    axis.text.x=element_text(color="black", size=10, hjust=1, face="bold", angle=45), # 自定义x轴上的文本
    axis.line.x.top=element_line(color="black"),  # 自定义顶部x轴线    
    axis.text.x.top=element_blank(), # 删除顶部x轴上的文本
    axis.ticks.y.right=element_blank(),  # 删除右侧y轴刻度
    axis.text.y.right=element_blank(),  # 删除右侧y轴文本
    axis.ticks.x.top=element_blank(), # 删除顶部x轴刻度
    panel.spacing.x=unit(0, "cm"), # 删除x轴间距
    panel.border=element_blank(), # 删除面板边框
    legend.position="none", # 删除图例
    panel.spacing=unit(0, "lines")) +  # 删除面板之间的间距
  guides(x.sec="axis", y.sec="axis")
print(p2)




set.seed(123)
fm <- matrix(rnorm(100), nrow=200,ncol = 8) %>% as.data.frame()
colnames(fm) <- c(paste0("Gene_", 1:8))
rownames(fm) <- paste0("Sample_", 1:200)
fm$alt_subtype <- rep(c('type1','type2','type3','type4'),time=c(50,50,50,50))
fm$group <- rep(c('group1','group2'),time=c(100,100))

fm$alt_subtype %>% unique()
# 为每个不同的组设置颜色。
col <- c("#5CB85C", "#337AB7", "#F0AD4E", "#D9534F")
comparisons <- list(c("type1", "type2"),
                    c("type1", "type3"),
                    c("type1", "type4"),
                    c("type2", "type3"),
                    c("type2", "type4"),
                    c("type3", "type4"))
plist <- list() # 创建一个空的列表 'plist'，用于存储后续循环中生成的图形
colnames(fm)

for (i in 1:8) {
  bar_tmp <- fm[, c(colnames(fm)[i], "alt_subtype")]  # 从 'Exp_plot' 中提取当前基因的表达信息和样本组
  colnames(bar_tmp) <- c("Expression", "group") 
  pb1 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                   x = "group", # X-axis is for groups.
                   y = "Expression", # Y-axis is for expression levels.
                   color = "group", # Fill by sample group.
                   fill = NULL,
                   add = "jitter", # Add jitter points.
                   bxp.errorbar.width = 0.8,
                   width = 0.5,
                   size = 0.1,
                   font.label = list(size = 20), 
                   palette = col) +
    theme(panel.background = element_blank())
  pb1 <- pb1 + theme(axis.line = element_line(colour = "black")) + 
    theme(axis.title.x = element_blank()) # 调整坐标轴
  pb1 <- pb1 + theme(axis.title.y = element_blank()) + 
    theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1)) 
  pb1 <- pb1 + theme(axis.text.y = element_text(size = 15)) + 
    ggtitle(colnames(fm)[i]) + 
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) 
  pb1 <- pb1 + theme(legend.position = "NA")  # 删除图例（因为样本类型已经显示在横轴上）
  pb1 <- pb1 + stat_compare_means(method = "t.test", hide.ns = FALSE, 
                                  comparisons = comparisons,label = "p.signif",vjust=0.02,bracket.size=0.6) # 执行显著性测试，使用 t 检验，并添加不同组别之间的比较
  plist[[i]] <- pb1 # 将生成的图形存储在 'plist' 中
}

# Align and arrange the plots into a grid.
p3 = plot_grid(plist[[1]], plist[[2]], plist[[3]],
          plist[[4]], plist[[5]], plist[[6]],
          plist[[7]], plist[[8]],ncol = 4) # ncol = 4 indicates the number of columns in the grid.
print(p3)




# 像这种多组之间两两比较，需要提前设置一个list，里面包含的是两两比较的对象
list = list(c("setosa","versicolor"),c("virginica","versicolor"),c("setosa","virginica")) 

# 把数据转成纵性数据
df=reshape2::melt(iris, # 需要转换的数据
                  id.vars=5, # 第5列Species不动
                  measure.vars=1:4) # 第1列到第4列全部转换成纵性数据

p4 = ggplot(df, # 用来画图的数据集名称
       aes(Species, # 数据集的Species列作为横坐标
           value,  # 数据集的value列作为纵坐标
           color = Species, # 根据数据集的Species列设置箱线图的边框颜色
           fill = Species))+ # 根据数据集的Species列设置散点图和箱线图的填充颜色
  geom_jitter(width=.15)+ # 将散点图的宽度设置为0.15,缩窄一些
  geom_boxplot(width=.4, alpha=.2)+ # 将箱线图的宽度设置为0.4，缩窄了一些，不透明度设置为0.2，填充颜色变得很浅
  scale_color_manual(values = c("#BF5960","#6F99AD","#F9A363"))+ # 自定义边框颜色
  scale_fill_manual(values = c("#BF5960","#6F99AD","#F9A363"))+ # 自定义填充颜色，可以和边框颜色不一样，我觉得配套的比较好看，就用了一样的颜色
  theme_test()+ # 把背景颜色去掉并添加边框，可以修改为其他背景，比如theme_classic
  labs(y="Length of petal", x="")+ # 修改横纵坐标名称
  stat_compare_means(comparisons = list, method = "wilcox.test", #像这种非配对的非正态分布的数据，可以用非参数检验方法，这里选的是非参的其中一种，叫秩和检验
                     size=4, # 将p值字体改小一些
                     label = "p.signif")+ # 把p值改为星号
  scale_y_continuous(expand = c(0.1,0.1))+ # 将y轴上下拉宽一点，这样被遮住的p值就可以呈现出来啦~
  facet_wrap(.~variable, # 同时做几个图，这几个图分别对应数据转换前的第1列到第4列
             nrow=2, ncol=2, # 图排成两行两列
             scale="free")+ # 横纵坐标不要统一成一样的，各个图根据自己情况来自动设置横纵坐标的范围
  theme(strip.background = element_blank(), #把头上的灰框去掉
        axis.text = element_text(color="black", size=10),
        axis.line = element_blank()) # 把横纵坐标的字体颜色改为黑色，看起来更沉稳
print(p4)




set.seed(123)
fm <- iris[,c(1,5)]
fm$Species %>% unique()
# fm$group <- rep(c('group1','group2'),time=c(100,100))
df2=fm
#统计平均值，以Group为基础统计，分别为AA，BB，CC，DD。
mean_df<-aggregate(df2[,1],by=list(df2[,2]),FUN=mean)
rownames(mean_df)<-mean_df[,1]
#统计标准差，以Group为基础统计，分别为AA，BB，CC，DD。
sd_df<-aggregate(df2[,1],by=list(df2[,2]),FUN=sd)
rownames(sd_df)<-sd_df[,1]
sd_df<-sd_df[,-1]
#计算标准误，样品数量为4，因此标准差为除以2，获得se值。
sd_df<-sd_df/2 
se<-as.data.frame(sd_df)
#合并数据框，以及标准误。
df1<-cbind(mean_df,se)
colnames(df1)<-c("group","mean","se")
colnames(df2)<-c("Gene_Abundance","group")
#绘图，geom_bar添加平均值柱子，geom_jitter用df2数据框加散点，geom_errorbar添加标准误误差棒，geom_hline添加虚线。
p5 = ggplot()+geom_bar(data=df1,mapping=aes(x=group,y=mean),fill = "white",size = 2, color = c("#4169B2","#479E9B","#BB2BA0"),
                  position="dodge", stat="identity",width = 0.65)+  
  geom_jitter(data=df2,mapping=aes(x=group,y=Gene_Abundance,fill = group,colour = group),
              size = 2,height = 0.02,width = 0.1)+ 
  ggsignif::geom_signif(comparisons = list(c('setosa','versicolor'),
                                           c('setosa','virginica'),
                                           c('versicolor','virginica')))+
  scale_color_manual(values = c("#4169B2","#479E9B","#BB2BA0"))+ 
  geom_errorbar(data=df1,mapping=aes(x = group,ymin = mean-se, ymax = mean+se), width = 0.3, 
                color = c("#4169B2","#479E9B","#BB2BA0"),size=1)+
  labs(y="Nitrogen metabolism", x="")+
  geom_hline(aes(yintercept =3),linetype="dashed", size=1.2, colour="gray56")+
  theme_classic(base_line_size = 1.05,base_rect_size =1.05)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=12))
print(p5)




library(ggh4x)
set.seed(123)
dfb <- iris[,c(1,5)]

dfb %>% 
  pivot_longer(!Species) %>% 
  select(Species,value) %>% 
  na.omit() -> dfb.1
library(ggplot2)
ebtop<-function(x){
  return(mean(x)+sd(x))
}
ebbottom<-function(x){
  return(mean(x)-sd(x))
}
p6 = ggplot(data=dfb.1,aes(x=Species,y=value))+
  stat_summary(geom = "bar",
               fun = mean,
               fill="#c6c3c3")+
  stat_summary(geom = "errorbar",
               fun.min = ebbottom,
               fun.max = ebtop,
               width=0.2)+
  geom_jitter(width = 0.3)+
  geom_signif(comparisons = list(c('setosa','versicolor'),
                                 c('versicolor','virginica'),
                                 c('setosa','virginica')
  ),
  test = t.test,
  test.args = list(var.equal=T,
                   alternative="two.side"),
  y_position = c(8,8.5,9.2),#自助修改
  annotations = c(""),
  parse = T)+
  annotate(geom = "text",
           x=1.5,y=8.65,
           label=expression(italic(P)~'='~2.22%*%10^-16))+
  annotate(geom = "text",
           x=2.5,y=9.15,
           label=expression(italic(P)~'='~1.7%*%10^-7))+
  annotate(geom = "text",
           x=2,y=9.8,
           label=expression(italic(P)~'='~2.22%*%10^-16))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,10),
                     breaks = seq(0,10,2))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(hjust=8,
                                    size=15),
        axis.text.x = element_text(angle = 30,
                                   hjust = 1,
                                   size=10))+
  guides(y=guide_axis_truncated(trunc_lower = 0,
                                trunc_upper = 10))+
  labs(x=NULL,y="Survival Rate")
print(p6)




library(ggplot2)
library(ggpubr)
library(RColorBrewer)
set.seed(123)
df <- iris[,c(1,5)] %>% subset(Species %in% c('setosa','versicolor'))
df$pair=rep(paste0('pair',seq(1:50)),time=2)
df <- set_colnames(df,c('Expression','Group1','Group2'))
pair <- c("A" = "red", "B" = "blue", "C" = "green")
p7 = ggplot(df,aes(x=Group1,y=Expression,color=Group1))+
  geom_boxplot(aes(color=Group1),width=0.8,lwd=1,outlier.shape=NA)+#箱线图
  scale_color_manual(limits=c("setosa","versicolor"),values=c("#993365","#666632"))+ #设置颜色
  geom_jitter(size=3.5,shape=21,alpha=0.8, #加散点，设置散点大小、形状、透明度
              aes(fill=Group2),
              position=position_dodge(0.5))+
  #scale_fill_brewer(palette='Set3')+
  stat_compare_means(method="t.test",paired=TRUE,#配对t检验
                     comparisons=list(c("setosa","versicolor")))+##按分组进行统计检验
  geom_line(aes(group=Group2),
            color='grey40',lwd=0.6,#线条粗细
            position=position_dodge(0.5))+
  labs(x='',y='Expression',title='',subtitle='')+
  ##以下是ggplot2的主题设置，修改边框、背景、标题、字体等
  theme_classic()+#主题设置
  theme(axis.line=element_line(colour='black',size=1),#坐标轴线条颜色、粗细 
        axis.text=element_text(size=14,color='black'),#坐标轴字体大小、颜色
        axis.text.x=element_text(angle=45,hjust=1),legend.position = '',
        axis.title=element_text(size=20,color='black'))
print(p7)




library(tidyverse)
library(ggtext)
library(ggprism)
library(ggsignif)
library(rstatix)
library(ggpubr)
df <- iris %>% pivot_longer(-Species) %>% 
  dplyr::filter(Species=="setosa") %>% select(-Species)
result.aov <- aov(value ~ name, data = df)
result.tukey <- TukeyHSD(result.aov)

# 转换p值
aov_pvalue <- result.tukey$name %>% as.data.frame() %>% 
  rownames_to_column(var="group") %>% 
  dplyr::select(1,`p adj`) %>% 
  separate(`group`, into=c("group2", "group1"), sep="-", convert = TRUE) %>% 
  select(2,1,3) %>% 
  select(-1,-2) %>% 
  mutate(p_signif=symnum(`p adj`, corr = FALSE, na = FALSE,  
                         cutpoints = c(0, 0.01, 0.05,1), 
                         symbols = c("**", "*", "ns")))

df_pvalue <- df %>% 
  wilcox_test(value ~ name) %>% 
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="name") %>% select(-p.adj) %>% 
  bind_cols(aov_pvalue)

 
p8= ggplot(df, aes(name,value))+
  stat_summary(fun = mean,geom = "errorbar", width=.2,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "crossbar",width = 0.4,color = "black",size=0.5) +
  geom_jitter(aes(fill=name,color=name,shape=name),width = 0.1, height = 0)+
  stat_pvalue_manual(df_pvalue,label="p_signif",label.size=5,hide.ns=T)+
  scale_shape_manual(values = c(21,22,23,24)) +
  scale_fill_manual(values=c("#679289","#ee2e31","#c9cba3","#f4c095"))+
  scale_color_manual(values=c("#679289","#ee2e31","#c9cba3","#f4c095"))+
  scale_y_continuous(guide = "prism_minor",
                     limits = c(0, 10),
                     expand = c(0, 0))+
  labs(x=NULL,y="setosa (10<sup>3</sup>/g)")+
  theme_prism()+
  theme(strip.background = element_blank(),
        legend.background = element_rect(color=NA),
        legend.key = element_blank(),
        legend.spacing.x = unit(-0.09,"in"),
        legend.spacing.y = unit(-0.09,"in"),
        legend.text = element_text(color="black",size=6,face="bold"),
        legend.position = "top",
        axis.text.x=element_blank(),
        axis.text.y=element_text(color="black",size=8),
        axis.title.y = element_markdown(color="black",face="bold",size=10),
        axis.ticks.x = element_blank())+
  guides(shape = guide_legend(override.aes = list(size=3),
                              direction = "horizontal",
                              nrow=3, byrow=TRUE))
print(p8)




#宽数据转成长数据
dd2 <-iris %>% 
  pivot_longer(cols=1:4,
               names_to= "Immune",
               values_to = "Fraction")

p9 = ggplot(dd2,aes(x=Immune, y=Fraction, 
               fill = Species, color = Species)) + 
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c('#CD534C','#0073C2','#EFC000')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")
print(p9)




p10 = ggplot(data = iris[,c(1,5)]) + 
  geom_boxplot(mapping = aes(x = Species, y = Sepal.Length, colour = Species),  # 箱线图的映射
               #alpha = 0.5,  # 箱线图的透明度
               size = 1.5,  # 箱线图边界的粗细
               width = 0.6) +  # 箱线图的宽度
  geom_jitter(mapping = aes(x = Species, y = Sepal.Length, colour =Species),  # 散点图的映射
              size = 1.5,alpha = 0.3) +  # 散点的大小和透明度
  scale_color_manual(limits = c("setosa", "versicolor", "virginica"),  # 手动设置颜色的范围
                     values = c("#85B22E", "#5F80B4", "#E29827")) +  # 不同样本的颜色
  geom_signif(mapping = aes(x = Species, y = Sepal.Length),  # 显著性线的映射
              comparisons = list(c("setosa", "versicolor"), c("setosa", "virginica"), c("versicolor", "virginica")),  # 两两比较
              map_signif_level = TRUE,  # 映射显著性水平
              tip_length = c(0, 0, 0),  # 显著性线的长度
              y_position = c(15, 15.8, 16.6),  # 显著性线的纵坐标位置
              size = 1,  # 显著性线的粗细
              textsize = 4,  # 显著性标记的大小
              test = "t.test") +  # 统计检验的类型
  theme_classic(base_line_size = 1) +  # 经典主题，轴线较粗
  labs(title = "Gene Expression by Sepal.Length", x = "Species", y = "Expression level") +  # 设置标题和轴标签
  theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5),  # 标题的外观
        axis.title.y = element_text(size = 15, color = "black", face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),  # 纵坐标轴标题的外观
        legend.title = element_text(color = "black", size = 15, face = "bold"),  # 图例标题的外观
        legend.text = element_text(color = "black", size = 10, face = "bold"),  # 图例标签的外观
        axis.text.x = element_text(size = 13, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0),  # 横坐标轴标签的外观
        axis.text.y = element_text(size = 13, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))  # 纵坐标轴标签的外观
print(p10)




p11 = ggplot(iris[,c(1,5)], aes(x=Species,y=Sepal.Length))+
  geom_violin(aes(fill = Species,color=Species), trim = FALSE)+
  geom_boxplot(aes(fill=Species),notch = F,width=0.3)+
  labs(x="", y="Sequencing depth")+
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  scale_fill_manual(values=c("#3EB185","#CCE3FA","#3F88C5","#F3AC66","#DD8CBC"))+ 
  scale_color_manual(values=c("#3EB185","#CCE3FA","#3F88C5","#F3AC66","#DD8CBC"))+ 
  annotate('text', label = 'Bacteria', x =3.5, y =8, angle=0, size =5,color="#DC9E05")+
  annotate("rect", xmin = 0, xmax =4,  ymin = 2.5, ymax = 9, alpha = 0.1,fill="#FAE3AD") +#"#CCE3FA"
  theme(axis.text=element_text(colour='black',size=10))+
  scale_y_continuous(expand = c(0, 0), limit = c(2.5, 9))
print(p11)

library(tidyverse)
library(ggbeeswarm)
library(scales)
library(ggfun)
set.seed(1115)

df <- tibble(
  `Day 0-Old` = abs(rnorm(100, mean = 10, sd = 10)),
  `Day 0-Young` =  abs(rnorm(100, mean = 30, sd = 30)),
  `Day 21-Old` =  abs(rnorm(100, mean = 60, sd = 60)),
  `Day 21-Young` = abs(rnorm(100, mean = 500, sd = 500)),
  `Day 35-Old` = abs(rnorm(100, mean = 1000, sd = 1000)),
  `Day 35-Young` = abs(rnorm(100, mean = 5000, sd = 5000))
)

# data cleaning
df_clean <- df %>%
  tidyr::gather(key = "key", value = "value") %>%
  dplyr::mutate(group = str_split(string = key, pattern="-", simplify = T)[,2])
# data statistics
df_clean_mean <- df_clean %>%
  group_by(key) %>%
  summarise(mean = mean(value),
            mean_scale = log10(mean))

df_clean_mean

# data signigicant
table(df_clean$key)

signif_out <- c()

for (i in c("Day 0","Day 21","Day 35")) {
  out <- t.test(df_clean %>% dplyr::filter(key == paste(i,"Old",sep = "-")) %>% pull(value),
                df_clean %>% dplyr::filter(key == paste(i,"Young",sep = "-")) %>% pull(value))
  signif_out <- c(signif_out, out$p.value)
}
# scale_y_log10
breaks_log10 <- function(x) {
  low <- floor(log10(min(x)))
  high <- ceiling(log10(max(x)))
  10^(seq.int(low, high))
}

p12 = ggplot(data = df_clean) + 
  geom_quasirandom(aes(x = key, y = value, shape = group, fill = group),
                   method = "pseudorandom", size = 3, alpha = 0.7) + 
  scale_shape_manual(values = c(21, 22)) + 
  scale_fill_manual(values = c("#80b1d3", "#fdb462")) + 
  annotate(geom = "segment", x = 0.6, xend = 1.4, y = 11.4, yend = 11.4, linewidth = 1) +
  annotate(geom = "segment", x = 1.6, xend = 2.4, y = 34.8, yend = 34.8, linewidth = 1) +
  annotate(geom = "segment", x = 2.6, xend = 3.4, y = 68.4, yend = 68.4, linewidth = 1) +
  annotate(geom = "segment", x = 3.6, xend = 4.4, y = 529, yend = 529, linewidth = 1) +
  annotate(geom = "segment", x = 4.6, xend = 5.4, y = 1179, yend = 1179, linewidth = 1) + 
  annotate(geom = "segment", x = 5.6, xend = 6.4, y = 6217, yend = 6217, linewidth = 1) +
  scale_y_log10(breaks = breaks_log10,
                labels = trans_format(log10, math_format(10^.x))) +
  annotation_logticks(sides = "l", outside = TRUE) + 
  coord_cartesian(clip = "off") + 
  scale_x_discrete(labels = c("Day 0", "Day 0", "Day 21", "Day 21", "Day 35", "Day 35")) +
  labs(x = "Sample", y = "Data") + 
  geom_hline(yintercept = 100, linetype = "dashed") + 
  # Day 0
  annotate(geom = "segment", x = 1, xend = 2, y = 500, yend = 500) + 
  annotate(geom = "segment", x = 1, xend = 1, y = 500, yend = 300) + 
  annotate(geom = "segment", x = 2, xend = 2, y = 500, yend = 300) + 
  annotate(geom = "text", x = 1.5, y = 800, label =  bquote(italic("p")~"< 0.01"), size = 5) +
  # Day 21
  annotate(geom = "segment", x = 3, xend = 4, y = 10000, yend = 10000) + 
  annotate(geom = "segment", x = 3, xend = 3, y = 10000, yend = 2000) + 
  annotate(geom = "segment", x = 4, xend = 4, y = 10000, yend = 8000) + 
  annotate(geom = "text", x = 3.5, y = 15000, label =  bquote(italic("p")~"< 0.01"), size = 5) +
  # Day 35
  annotate(geom = "segment", x = 5, xend = 6, y = 80000, yend = 80000) + 
  annotate(geom = "segment", x = 5, xend = 5, y = 80000, yend = 10000) + 
  annotate(geom = "segment", x = 6, xend = 6, y = 80000, yend = 60000) + 
  annotate(geom = "text", x = 5.5, y = 120000, label =  bquote(italic("p")~"< 0.01"), size = 5) +
  theme_classic() + 
  theme(axis.text = element_text(size = 15),
        axis.text.y.left = element_text(margin = margin(r = 10)),
        legend.background = element_roundrect(color = "#808080", linetype = 1))
print(p12)




library(tidyverse)
library(ggrain)
library(ggpp)
library(ggprism)
iris %>%
  dplyr::group_by(Species) %>%
  summarise(mean = round(mean(Sepal.Length),2),
            sd = round(sd(Sepal.Length),2),
            se = round(sd/sqrt(n()),2)) -> stat_df
stat_df_pos <- tibble(x = 0.5, y = 1, tb = list(stat_df))
# insert table
p13 = ggplot(data = iris, 
       aes(x = Species, 
           y = Sepal.Width, 
           fill = Species, 
           color = Species)
) + 
  geom_rain(alpha = 0.6,
            rain.side = "l") + 
  geom_table(data = stat_df_pos, aes(x = x, y = y, label = tb)) + 
  theme_bw()
print(p13)




ggplot(iris, aes(x = Sepal.Length, Sepal.Width, color = Species)) + 
  geom_point(size = 3) + 
  theme_bw() + 
  theme(legend.position = "none")-> p

plot_pos <- tibble(x = 3.5, y = 0.2, plot = list(p))

p14 = ggplot(data = iris, 
       aes(x = Species, y = Sepal.Width, fill = Species, color = Species)) + 
  geom_rain(alpha = 0.6, rain.side = "l") + 
  geom_table(data = stat_df_pos, aes(x = x, y = y, label = tb)) + 
  geom_plot(data = plot_pos, aes(x = x, y = y, label = plot)) +
  theme_bw()
print(p14)




library(ggbeeswarm)
df <- iris %>% pivot_longer(-Species)
p15 = ggplot(data=df, aes(x=name, y=value,color=name)) +
  geom_quasirandom(dodge.width = 1, alpha=0.8, size=0.6) +
  ggsci::scale_color_aaas(guide=FALSE) + 
  labs(x='Antibiotics', y='log2(% + 1)')+
  theme(axis.text.x = element_text(angle=40, hjust=1, size=8),
        strip.text = element_text(size=7, face='bold.italic')) +
  facet_wrap(~Species,nrow=2)+theme_classic()
print(p15)




library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
mytheme <- theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
                 axis.line = element_line(color = "black",size = 0.4),
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
                 axis.text.y = element_text(color="black",size=10),
                 axis.text.x = element_text(angle = 45, hjust = 1 ,color="black",size=10),
                 legend.position = "none",
                 panel.spacing = unit(0,"lines"))
fill.color = c(setosa="#56B4E9",versicolor= "#E69F00")
library(ggpubr)
library(RColorBrewer)
set.seed(123)
df <- iris %>% subset(Species %in% c('setosa','versicolor'))
df$pair=rep(paste0('pair',seq(1:50)),time=2)
df <- df %>% pivot_longer(-c(Species,pair))

p16 = ggplot(df, aes(x = Species,y = value,fill = Species)) +
  geom_line(aes(group=pair),position = position_dodge(0.2),color="grey80") +
  geom_point(aes(group=pair,size=value),
             pch=21,
             size=3,
             # position = position_dodge(0),
             position = position_jitter(w = 0.1))+
  stat_summary(fun.data = 'mean_se', geom = "errorbar", color = "red", 
               width = 0.25, size = 0.8, position = position_dodge( .9)) +
  # scale_size_continuous(range=c(1,3)) +
  facet_wrap(.~name,scales = "free_y",ncol = 11) +
  scale_fill_manual(values = fill.color) +
  # scale_y_continuous(limits = c(-4,4),minor_breaks = seq(0,90,1)) +
  labs(x= NULL,y="Gene expression")+
  stat_compare_means(aes(group =  Species),
                     paired = T,
                     # label.y = 3.5, 
                     label.x = 1.5,
                     method = "t.test",
                     label = "p.signif",
                     size = 4) + 
  theme_bw() + mytheme
print(p16)

#上面的P值是用stat_compare_means计算的，其实多组间的两两比较还可以考虑用校正后的P值，可以使用rstatix包进行计算：
stat.test<- df %>%
  group_by(name) %>%
  t_test(value~Species) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>%
  add_xy_position(x="Group",dodge = 0.8)
stat.test

####科学计数法显示保留前两位小数
# stat.test$p.adj <- format(stat.test$p.adj, scientific = TRUE) %>% as.numeric()
# stat.test$p.adj <- sprintf("%.2e", stat.test$p.adj)
# stat.test$p.adj = ifelse(as.numeric(stat.test$p.adj) > 0.05,"NS",stat.test$p.adj)
p17 = ggplot(df, aes(x = Species,y = value)) + 
  geom_line(aes(group=pair),position = position_dodge(0.2),color="grey80") +
  geom_point(aes(group=pair,size=value,
                 # alpha=Expression,
                 fill= Species),pch=21,size=3,
             position = position_jitter(w = 0.1))+
  stat_summary(fun.data = 'mean_se', geom = "errorbar", color = "red", 
               width = 0.25, size = 0.8, position = position_dodge( .9)) +
  facet_wrap(~name ,scale="free_y",nrow = 2) + 
  scale_fill_manual(values= fill.color) +
  scale_x_discrete( labels=NULL ) + 
  stat_pvalue_manual(stat.test,label = 'p.adj.signif', #p.adj.signif
                     tip.length = 0.05,
                     #y.position = 6,
                     remove.bracket = T,hjust=1)+
  labs(x= NULL,y="Gene expression")+
  theme_bw() + mytheme
print(p17)




library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(ggsignif) # 

df <- iris[,c(1,5)] %>% subset(Species %in% c('setosa','versicolor'))
colnames(df)=c('value','group')
df %>% 
  group_by(group) %>% 
  summarise(mean_value=mean(value)) %>% 
  bind_cols(x=c(1,2))-> df1

p18 = ggplot(df,aes(group,value))+
  #散点图
  geom_jitter(size=5,aes(fill=group),color="black",width = 0.08,shape=21)+
  #均值
  geom_segment(data=df1,aes(x=x-0.25,xend=x+0.25,y=mean_value,yend=mean_value),
               color="black",linewidth=1.5)+
  #显著性
  geom_signif(comparisons = list(c("setosa","versicolor")),
              map_signif_level=T, 
              tip_length=0, 
              y_position = 8, 
              size=1, textsize = 7, 
              test = "t.test")+
  #主题相关设置
  theme_classic()+
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 1),
        axis.text.y = element_text(color="black",size = 15),
        axis.text.x = element_text(color="black",size = 16,angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(color="black",size = 18),
        axis.ticks.y = element_line(size=1),
        axis.ticks.x = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(color="black",size = 20))+
  #颜色
  scale_fill_manual(values = c("#00c000","#a0ffa0"))+
  #标题
  labs(x=NULL,y="Percentage in arrest",title = "Arrest coefficier")+
  #轴范围
  scale_y_continuous(limits = c(0,10),breaks = seq(0,10, len = 6))
print(p18)




library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(ggsignif) # 

setosa<- "#636775"
versicolor <-"#F44686"
df <- iris %>% subset(Species %in% c('setosa','versicolor'))
df= df %>% pivot_longer(-c(Species))
colnames(df)
p19 = ggplot(df,aes(x=Species, y= value,fill=Species))+
  #geom_violin(aes(fill=condition), scale = "area", trim=F, size=0.5)+
  geom_violin(scale = "area", trim=F, size=0.5,color="white",cex=1,alpha=1)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="white")+
  #facet_grid(.~cell_types_3_groups, labeller = labeller(cell_types_3_groups=cell_types_labels))+
  facet_grid(.~name)+
  xlab("")+
  ylab("MG score")+
  scale_fill_manual(values=c(setosa, versicolor))+
  theme_classic(base_size=14)+
  ggsignif::geom_signif(comparisons = list(c("setosa", "versicolor")), # 第2组与第3组的比较
                        map_signif_level = F)
print(p19)




library(ggplot2) #绘图
library(ggsignif) #添加统计检验
library(ggdist) #云雨图
iris$Species
Vec1 <- c("setosa", "versicolor", "virginica")
comb_list <- list()
for(i in 1:(length(Vec1)-1)) {
  for(j in (i+1):length(Vec1)) {
    comb <- combn(c(Vec1[i], Vec1[j]), 2)
    if(!any(comb[1,] == comb[2,])) {
      comb_list[length(comb_list)+1] <- list(comb)
    }
  }
}
#选择一个颜色运行,如果修改数据，请根据你的分组数量，设置对应数量的颜色
Custom.color <- c("#d3838a","#efa48a","#cbd27e")
Custom.color <- c("#4a9a5b","#8ab522","#d1d628")
Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#ae3f51","#d6b55a","#79aec1")

p20 = ggplot(iris, aes(x = Species, y = Sepal.Length,fill=Species)) +
  geom_jitter(mapping = aes(color=Species),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Species),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  # ylim(11,35)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("Histological origins of cancer cells") +  #设置X轴标题
  ylab("Degree of cellular infiltration") +   #设置Y轴标题
  ggtitle("A visual case")+  #设置主标题
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  geom_signif(comparisons = comb_list,step_increase = .1,map_signif_level = TRUE,vjust = 0.5,hjust= 0)
print(p20)


library(patchwork)
pdf("test.pdf", width=20, height=30)
# print((p1+p2+p3+p4+p5)/(p6+p7+p8+p9+p10)/(p11+p12+p13+p14+p15)/(p16+p17+p18+p19+p20))
print(p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19+p20)
dev.off()
