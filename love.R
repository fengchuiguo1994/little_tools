plot(1:10,col="white",axes=FALSE, xlab = " ", ylab = " ")
for(i in 1:180){arrows(5,7,5,7.1,col="red",length =i/100, angle =180-i)}
arrows(2,8,7,4,col="red",lwd=2)
k=-0.8
for(i in 1:10){
  len=i/10
  arrows(2,8,2+len,8-0.8*len,col="red")
}

