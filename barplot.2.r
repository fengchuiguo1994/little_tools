"
sample	unmap	uniqmap	lowmap	multimap
rep1-DNA	11388802	140042857	4712919	12827873
rep2-DNA	10718863	136639337	4694716	12310375
"

"
Rscript barplot.2.r barplot.2.r.txt barplot.2.r.txt.pdf
"
library(reshape2)
library(ggplot2)
args=commandArgs(T)
mydat=read.table(args[1],header=T)
md <- melt(mydat,id=c("sample"),measure=c("unmap","uniqmap","lowmap","multimap"))
col = c("#f8b62b","#00a0e9","#009944","#601986")
pdf(args[2],width=10)
ggplot(data = md, mapping = aes(x = sample, y = value, fill = variable)) + # 解析数据
    geom_bar(stat = 'identity', position = 'fill',colour = 'black') + # 画图
    scale_fill_manual(values = col) +
    ylab('Fraction') + 
    theme_bw() + 
    theme(panel.grid=element_blank())
dev.off()