library(Seurat)
data("pbmc_small")
library(ggplot2)
library(ggmagnify)
library(magrittr)
dt <- pbmc_small@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(ident=Idents(pbmc_small))
#绘制常规图形；
p1 <- ggplot(dt, aes(tSNE_1, tSNE_2, color = ident)) +
  geom_point(alpha = 0.8, size = 1) +
  labs(x = "tSNE 1",y = "tSNE 2")
p1

#自定义配色；
color1 <- c("#00a8e1","#99cc00","#e30039","#fcd300","#800080","#00994e",
            "#ff6600","#808000","#db00c2","#008080","#0000ff","#c8cc00")
color2 <- c("#fa2c7b","#ff38e0", "#ffa235","#04c5f3","#0066fe","#8932a5",
            "#c90444", "#cb9bff","#434348","#90ed7d","#f7a35c","#8085e9")
#使用自定义颜色1,并调整图例图形（点）的大小；
p2 <- p1+scale_color_manual(values = color1)+
  guides(color=guide_legend(override.aes = list(size=5)))+
  theme_classic()
p2

p3 <- p1+scale_color_manual(values = color2)+
  guides(color=guide_legend(override.aes = list(size=5)))+
  theme_classic()
p3


#设置放大区域（左下、右上两个点的坐标）；
target = c(-40, -25, -25, -10)
#设置插入区域（左下、右上两个点的坐标）；
insert = c(0, 20, 10, 20)

#放大局部区域；
p4 <- p3 + geom_magnify(from = target, to = insert)
p4


#放大局部区域（添加坐标轴）,并自定义轮廓线颜色；
p5 <- p2 + geom_magnify(from = target,
                        to = insert,
                        proj ="facing",
                        colour="#0066fe",
                        linewidth = 0.5,
                        axes = "xy")
p5
