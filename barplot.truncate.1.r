library(ggplot2)
library(ggprism)
library(ggbreak) 

df <- data.frame(x = c('a','b','c','d','e','f','g','h','i','j'), y = c(rnorm(3) + 20, rnorm(3) + 10, rnorm(4) + 50) )
p1 <- ggplot(df,aes(x,y))+
  geom_col(aes(fill=x))+
  theme_prism(palette = "flames",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16,
              base_line_size = 0.8,
              axis_text_angle = 45)+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0))
p1

# 截断一次
p2<-p1+scale_y_break(c(30,40), # 截断位置及范围
                space = 0.3,   # 间距大小
                scales = 1.5)  # 上下显示比例，大于1上面比例大，小于1下面比例大
p2

# 截断两次
p3<-p1+scale_y_break(c(5,8),scales = 1.5,space = 0.3)+
  scale_y_break(c(40,45),scales = 1.5,space = 0.3)
p3

# 旋转图形并进行截断
p4<-p1+coord_flip() +
  scale_y_break(c(40,45),scales = 1.8,space = 0.3)
p4

(p1+p2)/(p3+p4)