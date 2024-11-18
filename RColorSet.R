library(smallstuff)
library(RColorBrewer)

cols = brewer.pal(11,'RdYlBu')
cols = brewer.pal(9, "YlOrRd") # cols<-brewer.pal(3, "YlOrRd") 
pal<-colorRampPalette(cols)
mycolors<-pal(16)
plotCol(mycolors)
color = colorRampPalette(c("#4DBBD5FF","#FF6F00FF","#E64B35FF"))(101)


library(ggsci)
library(scales)
cols = pal_npg()(10)
show_col(cols)


library(scales)
hex <- hue_pal()(3)
show_col(hex)


library(viridis)
scale_color_viridis() + # 配合ggplot2


c1a='#e0f3db' # 单色色系的热图配色
c2a='#a8ddb5'
c3a='#43a2ca'
c1='#edf8fb'
c2='#b2e2e2'
c3='#66c2a4'
c4='#2ca25f'
c5='#006d2c'
my_palette <- colorRampPalette(c(c1a, c1, c2, c4, c5))(n = 499)
col_breaks = c(seq(0.5,0.59,length=100),  # for c1
  seq(0.60,0.69,length=100),          # for c2
  seq(0.70,0.79,length=100),          # for c3
  seq(0.80,.89,length=100),           # for c4
  seq(0.90,1,length=100))             # for c5
