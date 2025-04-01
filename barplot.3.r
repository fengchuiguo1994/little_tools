library(ggplot2)
library(patchwork)
data <- data.frame(Experiment=c(rep('Control',6),rep('Exp1',6),rep('Exp2',6)), Type=factor(c(rep(1:6,3))), Value=round(rnorm(18,10,4),digits = 1))
knitr::kable(head(data,10))
p1 = ggplot(data,aes(x=Experiment,y=Value, fill=Experiment))+ 
  geom_bar(stat='identity') + # identity: 加总 # sum(data[data$Experiment=="Control",]$Value)
  labs(x=NULL)+             #修改坐标轴标题
  theme_bw(base_size = 18)+ #自定义主题等
  theme(axis.text = element_text(colour = 'black'), legend.position = 'none')

p2 = ggplot(data,aes(x=Experiment,y=Value,fill=Type))+  #根据不同Type填充颜色
  geom_bar(stat = 'identity')+ 
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+ 
  theme(axis.text = element_text(colour = 'black'))

p3 = ggplot(data,aes(x=Experiment,y=Value,fill=Type))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label=Value), size=5, position = position_stack(vjust = 0.5))+ #添加Value的文本注释 #调整文本位置
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))

p4 = ggplot(data,aes(x=Experiment,y=Value,fill=Type))+
  geom_bar(stat = 'identity', position = 'fill')+  #标准化成每组总量为1，也可以表示各类别所占百分比，修改一下y轴刻度即可
  geom_text(aes(label=Value), size=5,position = position_fill(vjust = 0.5))+
  labs(x=NULL,y='Proprotin of Value')+ 
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))

p5 = ggplot(data,aes(x=Experiment,y=Value,fill=Type))+
  geom_bar(stat = 'identity', position = 'dodge', width = 0.8, color='black')+ #设置柱子宽度,使变量之间分开 #柱状图位置并排: #使用position=position_dodge(width=0.9),可使组内柱子间隔,自行试一下。
  geom_text(aes(label=Value),size=4, position = position_dodge(width = 0.8), vjust=-0.3)+    #调节注释高度 #相应的注释宽度也调整
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))

p6 = ggplot(data,aes(x=Experiment,y=Value,fill=Type))+
  geom_bar(stat = 'identity',position = 'fill')+
  geom_text(aes(label=Value), size=5, position = position_fill(vjust = 0.5))+
  labs(x=NULL,y='Proprotin of Value')+ 
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  #红绿色盲人群友好的颜色设置：
  scale_fill_manual(values = rep(c('#D55E00','#0072B2','#F0E442','#009E73','#56B4E9','#CC79A7'),3)) #rep()函数重复作用，这里重复这组颜色3次

p7 = ggplot(data,aes(x=Experiment,y=Value,fill=Type))+
  geom_bar(stat = 'identity', position = 'dodge', width = 0.8,color='black')+        
  geom_text(aes(label=Value),size=4,position = position_dodge(width = 0.8), vjust=-0.3)+ 
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#D55E00','#0072B2','#F0E442','#009E73','#56B4E9','#CC79A7'),3))


design <- "AB
           CD
           EF
           G#"
# design <- "AB
# CD
# EF
# GG"
pdf("test.test.pdf", width=14, height = 12)
wrap_plots(A = p1, B = p2, C = p3, D = p4, E = p5, F = p6, G = p7, design = design)
dev.off()
