library(ggplot2)
args=commandArgs(T)

mydat = read.table(args[1],header=F)
names(mydat) = c("sample","pos","percentage")
mydat$per = mydat$percentage*100

pdf(args[2],height=5)
# ggplot(mydat, aes(x = pos, y=percentage, group = sample, fill = sample, colour = sample)) +
ggplot(mydat, aes(x = pos, y=per, group = sample, fill = sample, colour = sample)) +
  geom_line(size=0.75) +
  # scale_color_manual( values = c("#da1c6f","#da1c6f","#f8b62b","#f8b62b","#00a0e9","#00a0e9","#009944","#009944")) +
  # scale_linetype_manual(values = c("solid","dashed","solid","dashed","solid","dashed","solid","dashed")) +
  scale_x_continuous(expand = c(0,0),limits = c(-0.1,42), breaks = seq(1,41,by=5),labels=c("-20","-15","-10","-5","Junction","5","10","15","20")) + 
  scale_y_continuous(expand = c(0,0),limits = c(-0.5,100.5), breaks = seq(0,100, by=20)) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "#555555") +
  labs(title=args[1],x="position (bp)",y="Cytosine (%)") + 
  geom_vline(xintercept=21, linetype="dotted") +
  theme_bw() + 
  theme(panel.grid=element_blank())

dev.off()
