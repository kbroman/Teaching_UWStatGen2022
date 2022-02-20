# compare results with and without loco kinship

library(here)
library(qtl2)

do <- readRDS(here("QTLresults/do.rds"))
apr <- readRDS(here("R/_cache/aprobs.rds"))
k_loco <- readRDS(here("QTLresults/kinship.rds"))
k_noloco <- readRDS(here("QTLresults/kinship_noloco.rds"))

phe <- log(do$pheno[,2])
covar <- cbind(sex=(do$is_female*1),
               wbc=log10(do$pheno[,1]))

out_add <- readRDS(here("QTLresults/out_add.rds"))
out_noloco <- scan1(apr, phe, addcovar=covar, kinship=k_noloco, cores=0)
out_nok <- scan1(apr, phe, addcovar=covar, cores=0)

pmap <- readRDS(here("QTLresults/pmap.rds"))

pdf(here("Figs/qtl_scan_loco.pdf"), height=5.5, width=10, pointsize=10)
par(mfrow=c(2,1), mar=c(3.1, 4.1, 0.6, 0.6))

plot(out_nok, pmap, chr=c(1,4,7,10,14), col="green3", mgp.x=c(1.8, 0.2, 0))
plot(out_add, pmap, chr=c(1,4,7,10,14), add=TRUE, col="slateblue")
plot(out_noloco, pmap, chr=c(1,4,7,10,14), add=TRUE, col="violetred")
legend("topleft", c("no kinship", "loco", "single kinship"),
       lwd=2, col=c("green3", "slateblue", "violetred"),
       bg="gray90")

plot(out_nok - out_add, pmap, chr=c(1,4,7,10,14), col="green3", ylim=c(-2, 2),
     ylab="LOD difference", mgp.x=c(1.8, 0.2, 0))
plot(out_noloco - out_add, pmap, chr=c(1,4,7,10,14), add=TRUE, col="violetred",
     ylim=c(-2, 2))
legend("topleft", c("no kinship - loco", "single kinship - loco"),
       lwd=2, col=c("green3", "violetred"),
       bg="gray90")

dev.off()
