library(reshape2)
library(ggplot2)
args=commandArgs(T)

mydat=read.table(args[1],header=F)
names(mydat)=c("sample","variable","value")
mydat$sample = factor(mydat$sample, levels=unique(mydat$sample))
pdf(args[2])
ggplot(data = mydat, mapping = aes(x = sample, y = value, fill = variable)) + # 解析数据
    geom_bar(stat = 'identity', position = 'fill',colour = 'black') + # 画图
    # scale_fill_manual(values = col) +
    labs(x = "",y = "Fraction", title = args[1]) +
    theme_bw() + 
    theme(panel.grid=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) # 横轴倾斜45度
dev.off()
