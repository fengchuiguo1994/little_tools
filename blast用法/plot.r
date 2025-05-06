library(ggplot2)
mydat <- read.table("EUSNP4.map.tongji",header=F)
names(mydat) <- c("count","species")
ggplot(data=mydat,aes(x=species,fill=species))+geom_bar(stat="count")
