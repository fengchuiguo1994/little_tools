library(ggplot2)
args=commandArgs(T)

mydat = read.table(args[1],header=F,stringsAsFactors=F)
# names(mydat) = c("sample",'count')
names(mydat) = c("sample",'total','chimeric')
mydat$count = mydat$chimeric/mydat$total*100
mydat$sample = factor(mydat$sample, levels=mydat$sample)

pdf(args[2],width=7)
ggplot(mydat, aes(x = sample,y = count,fill="#ffb548",color = "#ffb548",))+
    #####这部分的position_dodge(width=0.8)大于宽width = 0.6点，可以使得分组内柱子之间有缝隙，而不是贴合。
    # geom_bar(stat ="identity",width = 0.6,position = position_dodge(width=0.8))+        
    geom_bar(stat = "identity", ,position='dodge',width = 0.6,size=0.25,alpha=1) +
    labs(x = "",y = "Chimeric reads percentage (%)", title = args[1])+        ############坐标标签和图片title
    guides(fill = guide_legend(reverse = F))+                           ##############图例顺序反转
    theme_bw() +
    theme(panel.grid=element_blank()) +
    # geom_text(aes(label = mydat$count),position=position_dodge(width = 0.9),size = 3,vjust = -0.25)+    ###########设置柱子上的标签文字，文字的position_dodge(width=0.5)设置，保证分隔宽度。
    geom_text(aes(y = count + 0.05,label=round(count,digits=2)), position = position_dodge(0.9), vjust = -0.5) +
    theme(# legend.title = element_blank(),                 ##########图例名称为空
        # legend.text = element_text(size = 18, face = "bold"),   ##########图例文字大小
        # legend.position = 'right',      ############图例位置
        # legend.key.size=unit(0.8,'cm'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) # 横轴倾斜45度
dev.off()
