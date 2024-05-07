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
par(fig=c(0,0.5,0.7,0.8), new=TRUE)
barplot(c(1:6),col = brewer.pal(9,"Set1")[3:8])
par(fig=c(0,0.1,0.8,0.9), new=TRUE)
pie(rep(1, 20), col = rainbow(20), main = "rainbow20")
par(fig=c(0.1,0.2,0.8,0.9), new=TRUE)
pie(rep(1, 1000), labels = NA, col = rainbow(1000), border = rainbow(1000), main = "rainbow1000")

par(fig=c(0.2,0.3,0.8,0.9), new=TRUE)
cols = colors()
pie(rep(1, 657), labels = NA, col = cols, border = cols, main = "colors")

par(fig=c(0.3,0.5,0.8,0.9), new=TRUE)
r <- seq(7, 255, 8)
col <- rgb(r, 0, 0, alpha = 255, maxColorValue = 255)
image(x = 1:8, y = 1:4, z = matrix(1:32, ncol = 4), col = col, axes = F, ann = F)


## ggsci
```
library("ggsci")
library("ggplot2")
library("gridExtra")

scale_color_palname() # scale_color_npg() scale_color_aaas()
scale_fill_palname() # scale_fill_npg() scale_fill_aaas()

mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal <- pal_aaas("default", alpha = 0.7)(9)
library("scales")
show_col(mypal)
```