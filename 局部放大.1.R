# https://mp.weixin.qq.com/s?__biz=MzkzMjQyOTA5Nw==&mid=2247484675&idx=1&sn=7f5c09a7bc1a4efbd00d0867cb49adc0&chksm=c25aa945f52d2053c855f81ad0190623199703d43f507e5c98dd52fd6799a23bf9ab6525fda0&cur_album_id=3181722267144339456&scene=190#rd

# remotes::install_github("hughjonesd/ggmagnify")
# library(devtools)
# devtools::install_local("ggmagnify-master.zip")
# install.packages("lifecycle")
# 从https://github.com/hughjonesd/ggmagnify网站下载ggmagnify软件包
# 本地安装

library(ggmagnify)
library(lifecycle)
library(ggplot2)
library(ggfx)

ggpi <- ggplot(iris, aes(Sepal.Width, Sepal.Length, colour = Species)) +
  geom_point() + xlim(c(1.5, 6))

ggpi + geom_magnify(aes(from = Species == "setosa" & Sepal.Length < 5), 
                    to = c(4, 6, 6, 7.5))

ggpi +  facet_wrap(vars(Species)) +
  geom_magnify(aes(from = Sepal.Length > 5 & Sepal.Length < 6.5), 
               to = c(4.5, 6, 6, 7.5),
               shadow = TRUE)
