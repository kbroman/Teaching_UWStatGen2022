# plot of genome scan, of snps and of haplotypes

library(here)
library(qtl2)

dir <- here("QTLresults")
pmap <- readRDS(file.path(dir, "pmap.rds"))
out_add <- readRDS(file.path(dir, "out_add.rds"))
out_snps <- readRDS(file.path(dir, "out_snps.rds"))
operm_add <- readRDS(file.path(dir, "operm_add.rds"))
operm_snps <- readRDS(file.path(dir, "operm_snps.rds"))
out_full <- readRDS(file.path(dir, "out_full.rds"))
operm_full <- readRDS(file.path(dir, "operm_full.rds"))

# line colors
altcolor <- "green4"
linecolor <- "violetred"


thr_full <- summary(operm_full)
thr_add <- summary(operm_add)
thr_snps <- summary(operm_snps)

ymx_add <- maxlod(out_add)*1.04
ymx_snps <- thr_snps$A/thr_add$A*maxlod(out_add)*1.04
ymx <- max(c(ymx_add, ymx_snps))
ymx2 <- max(unlist(thr_full))*1.04


res <- 255


endA <- xpos_scan1(pmap, thechr=19, thepos=max(pmap[[19]]))+25/2

png(here("Figs/qtl_scan.png"), height=6.5*res, width=13*res, pointsize=10, res=res)
par(mar=c(2.1,4.1,2.1,1.1))
layout(c(1,2), height=c(0.5, 0.5))

plot(out_snps$lod, out_snps$snpinfo, altcol=altcolor, xlab="",
     ylim=c(0, ymx))
mtext(side=3, adj=1, "SNP scan", cex=1.5, line=0.5, col="darkslateblue")
u <- par("usr")
segments(u[1], thr_snps$A, endA, thr_snps$A, col=linecolor, lty=2)
segments(endA, thr_snps$X, u[2], thr_snps$X, col=linecolor, lty=2)

plot(out_add, pmap, altcol=altcolor)
mtext(side=3, adj=1, "additive founder alleles", cex=1.5, line=0.5, col="darkslateblue")
segments(u[1], thr_add$A, endA, thr_add$A, col=linecolor, lty=2)
segments(endA, thr_add$X, u[2], thr_add$X, col=linecolor, lty=2)
dev.off()


png(here("Figs/qtl_scan_full.png"), height=6.5*res, width=13*res, pointsize=10, res=res)
par(mar=c(2.1,4.1,2.1,1.1))
layout(c(1,2), height=c(0.5, 0.5))

plot(out_full, pmap, altcol=altcolor, ylim=c(0, ymx2))
mtext(side=3, adj=1, "full genotype model", cex=1.5, line=0.5, col="darkslateblue")
segments(u[1], thr_full$A, endA, thr_full$A, col=linecolor, lty=2)
segments(endA, thr_full$X, u[2], thr_full$X, col=linecolor, lty=2)

plot(out_add, pmap, altcol=altcolor)
mtext(side=3, adj=1, "additive founder alleles", cex=1.5, line=0.5, col="darkslateblue")
segments(u[1], thr_add$A, endA, thr_add$A, col=linecolor, lty=2)
segments(endA, thr_add$X, u[2], thr_add$X, col=linecolor, lty=2)
dev.off()
