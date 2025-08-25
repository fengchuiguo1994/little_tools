library(ggplot2)
args=commandArgs(T) # args=c("all.c2n.cf2.fc2.featureCount.DESeq2.txt", "all.c2n.cf2.fc2.featureCount.DESeq2.txt.pdf")
# mydat = read.table("all.c2n.cf2.fc2.featureCount.DESeq2.txt") # all.c2n.cf2.fc2.featureCount.DESeq2.txt
mydat = read.table(args[1])

# 设置p_value和logFC的阈值
cut_off_qvalue = 0.05  #统计显著性
cut_off_logFC = 1      #差异倍数值

mydat$change = ifelse(mydat$padj <= cut_off_qvalue & abs(mydat$log2FoldChange) >= cut_off_logFC, ifelse(mydat$log2FoldChange > cut_off_logFC ,'Up','Down'), 'Stable')

# pdf("all.c2n.cf2.fc2.featureCount.DESeq2.txt.pdf")
pdf(args[2])
ggplot(
  mydat, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-cut_off_logFC, cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_qvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)", y="-log10 (padj)")+
  theme_bw()+
  # 图例
  theme(panel.grid=element_blank(), plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank())
dev.off()