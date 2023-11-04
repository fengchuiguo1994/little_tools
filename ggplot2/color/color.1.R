library(ggplot2)
library(RColorBrewer)
library(patchwork)

par(cex=0.5, mai=c(0.1,0.1,0.1,0.1))
par(fig=c(0,0.5,0,0.5))
display.brewer.all(type = "seq")

par(fig=c(0.5,1,0,0.5), new=TRUE)
display.brewer.all(type = "div")

par(fig=c(0.5,1,0.5,1), new=TRUE)
p3 = display.brewer.all(type = "qual")

par(fig=c(0,0.5,0.5,0.6), new=TRUE)
display.brewer.pal(9,"Reds")

# display.brewer.pal(所取颜色的个数,"调色板名称")
# brewer.pal(9,"Set1")[3:8]) # 使用颜色
par(fig=c(0,0.5,0.5,0.6), new=TRUE)
display.brewer.pal(9,"Reds")
par(fig=c(0,0.5,0.6,0.7), new=TRUE)
display.brewer.pal(3,"Reds") 
par(fig=c(0,0.5,0.7,1), new=TRUE)
barplot(c(1:6),col = brewer.pal(9,"Set1")[3:8])