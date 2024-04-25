library(smallstuff)
library(RColorBrewer)

cols = brewer.pal(11,'RdYlBu')
cols = brewer.pal(9, "YlOrRd") # cols<-brewer.pal(3, "YlOrRd") 
pal<-colorRampPalette(cols)
mycolors<-pal(16)
plotCol(mycolors)
color = colorRampPalette(c("#4DBBD5FF","#FF6F00FF","#E64B35FF"))(101)


library(ggsci)
cols = pal_npg()(10)


library(scales)
hex <- hue_pal()(3)


library(viridis)
scale_color_viridis() + # 配合ggplot2