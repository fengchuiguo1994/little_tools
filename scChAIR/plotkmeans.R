library(factoextra)
library(TreeDist)
library(caret)
args=commandArgs(T)
pdf(args[2])
mydat = read.table(args[1],header=T,check.names=F) # SCG0093.raw.genome.txt
rownames(mydat) = mydat[,1]
mydat = mydat[,-1]
# mydat[1:15,1:15]

normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

percentage <- function(x) {
  return (x / sum(x))
}

scalemydat = scale(mydat, center = TRUE, scale = TRUE) # Standard 标准化 # 中心化center：x-mean(x) # 标准化：中心化后除以标准差 (x-mean(x))/sd(x)
normalizemydat = sapply(mydat, normalize) # Normalize 归一化 0-1 (test <- preProcess(mydat, method = 'range')    testmydat <- predict(test, mydat)     table(testmydat == normalizemydat))
percentagemydat = sapply(mydat, percentage) # %
# scalemydat[1:15,1:15]

# (log2(50000) - 10) / 0.125 # bin44
# (log2(8000000) - 10) / 0.125 # bin104

scalemydatflt = scalemydat[44:138,]
# p1 = fviz_nbclust(t(scalemydatflt), kmeans, method = "wss", print.summary=T, k.max=40) + geom_vline(xintercept = 4, linetype = 2) + ggtitle("scale bin44-138 spearman kmeans")
p1 = fviz_nbclust(t(scalemydatflt), KMeansPP, method = "wss", print.summary=T, k.max=40) + ggtitle("scale bin44-138 KMeansPP wss")
print(p1)
p1 = fviz_nbclust(t(scalemydatflt), KMeansPP, print.summary=T, k.max=40) + ggtitle("scale bin44-138 KMeansPP")
print(p1)
p1 = fviz_nbclust(t(scalemydatflt), kmeans, method = "wss", print.summary=T, k.max=40) + ggtitle("scale bin44-138 kmeans wss")
print(p1)
p1 = fviz_nbclust(t(scalemydatflt), kmeans, print.summary=T, k.max=40) + ggtitle("scale bin44-138 kmeans")
print(p1)

scalemydatflt = scalemydat[44:143,]
p1 = fviz_nbclust(t(scalemydatflt), KMeansPP, method = "wss", print.summary=T, k.max=40) + ggtitle("scale bin44-143 KMeansPP wss")
print(p1)
p1 = fviz_nbclust(t(scalemydatflt), KMeansPP, print.summary=T, k.max=40) + ggtitle("scale bin44-143 KMeansPP")
print(p1)
p1 = fviz_nbclust(t(scalemydatflt), kmeans, method = "wss", print.summary=T, k.max=40) + ggtitle("scale bin44-143 kmeans wss")
print(p1)
p1 = fviz_nbclust(t(scalemydatflt), kmeans, print.summary=T, k.max=40) + ggtitle("scale bin44-143 kmeans")
print(p1)

normalizemydatflt = normalizemydat[44:138,]
p1 = fviz_nbclust(t(normalizemydatflt), KMeansPP, method = "wss", print.summary=T, k.max=40) + ggtitle("normalize bin44-138 KMeansPP wss")
print(p1)
p1 = fviz_nbclust(t(normalizemydatflt), KMeansPP, print.summary=T, k.max=40) + ggtitle("normalize bin44-138 KMeansPP")
print(p1)
p1 = fviz_nbclust(t(normalizemydatflt), kmeans, method = "wss", print.summary=T, k.max=40) + ggtitle("normalize bin44-138 kmeans wss")
print(p1)
p1 = fviz_nbclust(t(normalizemydatflt), kmeans, print.summary=T, k.max=40) + ggtitle("normalize bin44-138 kmeans")
print(p1)

normalizemydatflt = normalizemydat[44:143,]
p1 = fviz_nbclust(t(normalizemydatflt), KMeansPP, method = "wss", print.summary=T, k.max=40) + ggtitle("normalize bin44-143 KMeansPP wss")
print(p1)
p1 = fviz_nbclust(t(normalizemydatflt), KMeansPP, print.summary=T, k.max=40) + ggtitle("normalize bin44-143 KMeansPP")
print(p1)
p1 = fviz_nbclust(t(normalizemydatflt), kmeans, method = "wss", print.summary=T, k.max=40) + ggtitle("normalize bin44-143 kmeans wss")
print(p1)
p1 = fviz_nbclust(t(normalizemydatflt), kmeans, print.summary=T, k.max=40) + ggtitle("normalize bin44-143 kmeans")
print(p1)

percentagemydatflt = percentagemydat[44:138,]
p1 = fviz_nbclust(t(percentagemydatflt), KMeansPP, method = "wss", print.summary=T, k.max=40) + ggtitle("percentage bin44-138 KMeansPP wss")
print(p1)
p1 = fviz_nbclust(t(percentagemydatflt), KMeansPP, print.summary=T, k.max=40) + ggtitle("percentage bin44-138 KMeansPP")
print(p1)
p1 = fviz_nbclust(t(percentagemydatflt), kmeans, method = "wss", print.summary=T, k.max=40) + ggtitle("percentage bin44-138 kmeans wss")
print(p1)
p1 = fviz_nbclust(t(percentagemydatflt), kmeans, print.summary=T, k.max=40) + ggtitle("percentage bin44-138 kmeans")
print(p1)

percentagemydatflt = percentagemydat[44:143,]
p1 = fviz_nbclust(t(percentagemydatflt), KMeansPP, method = "wss", print.summary=T, k.max=40) + ggtitle("percentage bin44-143 KMeansPP wss")
print(p1)
p1 = fviz_nbclust(t(percentagemydatflt), KMeansPP, print.summary=T, k.max=40) + ggtitle("percentage bin44-143 KMeansPP")
print(p1)
p1 = fviz_nbclust(t(percentagemydatflt), kmeans, method = "wss", print.summary=T, k.max=40) + ggtitle("percentage bin44-143 kmeans wss")
print(p1)
p1 = fviz_nbclust(t(percentagemydatflt), kmeans, print.summary=T, k.max=40) + ggtitle("percentage bin44-143 kmeans")
print(p1)

dev.off()
