# illustration of hidden markov model

ocolor <- "slateblue"
gcolor <- "green4"
arrowcolor <- "violetred"
xgap <- 5

library(here)

pdf(here("Figs/hmm.pdf"), height=2.5, width=11, pointsize=24)
par(mar=rep(0.1, 4), bty="n", cex=1.3)
plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n",
     xlim=c(-xgap, 14+xgap), ylim=c(8, 92), xaxs="i")

gx_pts <- c(0, 2, 6, 8, 12, 14)
gx_arr <- c(1, 3, 5, 7, 9, 11, 13)
gx_dots <- c(4, 10)
gy <- 25
gy_txt <- 18

u <- par("usr")
oy <- mean(u[3:4])+gy
oy_txt <- oy + (gy-gy_txt)

points(gx_pts, rep(gy, length(gx_pts)), pch=21, col=gcolor, bg=gcolor, lwd=2)
points(gx_pts, rep(oy, length(gx_pts)), pch=21, col=ocolor, bg=ocolor)

arrow_width <- 0.75

arrows(gx_arr - arrow_width/2, rep(gy, length(gx_arr)),
       gx_arr + arrow_width/2, rep(gy, length(gx_arr)), len=0.1, col=arrowcolor, lwd=2)

arrows(gx_pts, rep(gy + (oy-gy)*0.2, length(gx_arr)),
       gx_pts, rep(oy - (oy-gy)*0.2, length(gx_arr)),
       len=0.1, col=arrowcolor, lwd=2)

text(gx_dots, rep(gy, length(gx_dots)), "...", col=gcolor)

text(14.5+xgap/2, gy, "(unobserved)\ndiplotypes", col=gcolor, cex=0.7)
text(-xgap/2-0.5, oy, "SNP\ngenotypes", col=ocolor, cex=0.7)

text(gx_pts-0.15, rep(gy_txt, length(gx_pts)),
     c(expression(G[1]), expression(G[2]),
       expression(G[i]), expression(G[i+1]),
       expression(G[m-1]), expression(G[m])),
     adj=c(0, 1),
     col=gcolor, cex=0.5)

text(gx_pts-0.22, rep(oy_txt, length(gx_pts)),
     c(expression(O[1]), expression(O[2]),
       expression(O[i]), expression(O[i+1]),
       expression(O[m-1]), expression(O[m])),
     adj=c(0, 0),
     col=ocolor, cex=0.5)


dev.off()
