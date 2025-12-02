# 设置布局：3行3列
par(mfrow = c(3, 3))
slices <- c(12545193, 192564, 2444290)
lbls <- c("<=1kb FR", "<=1kb RF+Tandem", ">1kb")
pie(slices, labels = lbls, main="各国市场份额")

# 计算百分比
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # 添加百分比到标签
lbls <- paste(lbls, "%", sep="") # 添加%符号
# 绘制带百分比的饼图
pie(slices, labels = lbls, col=rainbow(length(lbls)), main="各国市场份额(百分比)")


library(plotrix)
pie3D(slices, labels=lbls, explode=0.1, main="3D饼图示例") # 绘制3D饼图
# my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
pie(slices, labels = lbls, col = my_colors, main="自定义颜色的饼图")


slices <- c(494438575, 2876199, 97989989)
lbls <- c("<=1kb FR", "<=1kb RF+Tandem", ">1kb")
# 计算百分比
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # 添加百分比到标签
lbls <- paste(lbls, "%", sep="") # 添加%符号
# my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
pie(slices, labels = lbls, col = my_colors, main="TestTime3/scHi-C/1280M")

slices <- c(170636156, 3971, 173658)
lbls <- c("<=1kb FR", "<=1kb RF+Tandem", ">1kb")
# 计算百分比
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # 添加百分比到标签
lbls <- paste(lbls, "%", sep="") # 添加%符号
# my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
pie(slices, labels = lbls, col = my_colors, main="TestTime3/scATAC/320M/")

slices <- c(189845565, 1692125, 14090721)
lbls <- c("<=1kb FR", "<=1kb RF+Tandem", ">1kb")
# 计算百分比
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # 添加百分比到标签
lbls <- paste(lbls, "%", sep="") # 添加%符号
# my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
pie(slices, labels = lbls, col = my_colors, main="TestTime3/scChIATAC/640M/")
