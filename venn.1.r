library(eulerr)
dat <- c("GroupA" = 25, 
         "GroupB" = 5, 
         "GroupC" = 5,
         "GroupA&GroupB" = 5, 
         "GroupA&GroupC" = 5, 
         "GroupB&GroupC" = 3,
         "GroupA&GroupB&GroupC" = 3)
p1 = plot(euler(dat))

p2 = plot(euler(dat),
     fills = list(fill=c("red","blue",
                         "green","darkgreen",
                         "white","black",
                         "purple"),
                  alpha=0.5))
print(p2)

p3 = plot(euler(dat),
          fills = list(fill=c("red","blue",
                              "green","darkgreen",
                              "white","black",
                              "purple"),
                       alpha=0.5),
          quantities = c(25,5,5,
                         5,5,3,3))
print(p3)

p4 = plot(euler(dat),
          fills = list(fill=c("red","blue",
                              "green","darkgreen",
                              "white","black",
                              "purple"),
                       alpha=0.5),
          quantities = list(c(25,5,5,
                              5,5,3,3),
                            col="black",
                            cex=2),
          labels = list(col="white",font=3,cex=2))
print(p4)

p5 = plot(euler(dat),
          fills = list(fill=c("red","blue",
                              "green","darkgreen",
                              "white","black",
                              "purple"),
                       alpha=0.5),
          quantities = list(c(25,5,5,
                              5,5,3,3),
                            col="black",
                            cex=4),
          labels = list(col="white",font=3,cex=2),
          edges = list(col="white",alpha=0))
print(p5)

p6 = plot(euler(dat),
          fills = list(fill=c("red","blue",
                              "green","darkgreen",
                              "white","black",
                              "purple"),
                       alpha=0.5),
          quantities = list(c(25,5,5,
                              1,1,1,1),
                            col="black",
                            cex=2),
          labels = list(col="white",font=3,cex=1),
          edges = list(col="darkgreen",lwd=5,
                       lty=1:3))
print(p6)

p7 = plot(euler(dat),
          fills = list(fill=c("red","blue",
                              "green","darkgreen",
                              "white","black",
                              "purple"),
                       alpha=0.5),
          quantities = list(c(25,5,5,
                              1,1,1,1),
                            col="black",
                            cex=2),
          labels = list(col="white",font=3,cex=1),
          edges = list(col="darkgreen",lwd=5,
                       lty=1:3),
          main = list(label=c("ABC"),cex=1),
          legend = list(labels=c("GroupA","GroupB","GroupC"),
                        cex=2))
print(p7)

plot(euler(dat),
     fills = list(fill=c("red","blue",
                         "green","darkgreen",
                         "white","black",
                         "purple"),
                  alpha=0.5),
     quantities = c(25,5,5,
                    5,5,3,3)) -> p1
plot(euler(dat),
     fills = list(fill=c("red","blue",
                         "green","darkgreen",
                         "white","black",
                         "purple"),
                  alpha=0.5),
     quantities = list(c(25,5,5,
                         5,5,3,3),
                       col="black",
                       cex=1),
     labels = list(col="white",font=3,cex=1),
     edges = list(col="darkgreen",lwd=5,
                  lty=1:3),
     main = list(label=c("ABC"),cex=1),
     legend = list(labels=c("GroupA","GroupB","GroupC"),cex=2)) -> p2
p8 = gridExtra::grid.arrange(p1,p2,ncol=2)
print(p8)

gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=3)
