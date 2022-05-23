args = commandArgs(T)
mydat = read.table(args[1])
pdf(args[2])
if (length(args) < 2){
    stop("the commans is like:\n\t\t\t\tRscript file.in file.in.pdf [height]\n")
} else if (length(args) == 2){
    height = 2
} else {
    height = as.numeric(args[3])
}
names(mydat) = c("gene",'seedling','panicle','leaf','entropy')

plot(c(0,2),c(0,height), type = "n", xlab = "Shannon Entropy",ylab = "denistity",xaxt="n")
dd = density(mydat$entropy)
lines(dd, col = rainbow(15)[2], lwd = 2)
axis(1,at=c(0,0.5,0.6,1,1.5,2),labels = c("0","0.5","0.6","1","1.5","2"))
abline(v=c(0.6),col="grey",lty=2)

dev.off()
