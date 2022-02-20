bgcolor <- "white"

pdf(file="../Figs/permtest.pdf", width=10, height=5, pointsize=12, onefile=TRUE)
par(fg="black",col="black",col.axis="black",col.lab="black",
    bg=bgcolor, bty="n")


par(mar=c(1.1, 1.1, 1.1, 1.1), cex=1.5)

plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(2,103),ylim=c(0,9))
u <- par("usr")
x0 <- 5; x1 <- 25; y1 <- 8; xm <- mean(c(x0,x1))
segments(c(x0,x0,x0,x1),c(0,0,y1,0),c(x1,x0,x1,x1),c(0,y1,y1,y1))
text(xm,y1/2,"genotype\ndata")
y2 <- mean(c(u[4],y1))
text(xm,y2,"markers")
text(x0/2,y1/2,"individuals",srt=90)
x2 <- 27.5; x3 <- 32.5
segments(c(x2,x2,x2,x3),c(0,0,y1,0),c(x3,x2,x3,x3),c(0,y1,y1,y1))
text(mean(c(x2,x3)),y1/2,"phenotypes",srt=90)
x4 <- 50; x5 <- 70; y3 <- 3; y4 <- 5; xm <- mean(c(x4,x5)); ym <- mean(c(y3,y4))
segments(c(x4,x4,x4,x5),c(y3,y3,y4,y3),c(x5,x4,x5,x5),c(y3,y4,y4,y4))
text(xm,ym,"LOD scores")
arrows(x3+3, ym, x4-3, ym, length=0.1)
x6 <- 87.5; x7 <- 97.5
arrows(x5+3, ym, x6-3, ym, length=0.1)
text(x7,ym,"maximum\nLOD score")

dev.off()
