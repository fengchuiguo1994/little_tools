aa = c(1,2,3,4,5,6,7,8,9,10,11)
bb = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0)
mydat = data.frame(count=aa,per=bb)
plot(x=NULL, y=NULL,
     type=n,
     main=paste(Insert Size Histogram for, leve, nin file),
     xlab=Insert Size,
     ylab=Count,
     xlim=range(0, 11),
     ylim=range(0, 1))
colors - c()
labels - c()
# lines(mydat$count, mydat$per,  type=h, col=red)
lines(mydat$count, mydat$per, col=red)
colors - c(colors, red)
labels - c(labels, FR)
legend(topright, labels, fill=colors, col=colors, cex=0.7)
