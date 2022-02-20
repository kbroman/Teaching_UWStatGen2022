library(qtl)

bgcolor <- "white"
color <- c("violetred", "slateblue")

pdf(file="../Figs/lodsupint.pdf", width=9, height=6, pointsize=12, onefile=TRUE)
par(fg="black",col="black",col.axis="black",col.lab="black", #bty="n",
    bg=bgcolor)#, cex.axis=1.5, cex.lab=1.5)

file <- "_cache/insulin_lod.RData"
if(file.exists(file)) {
  load(file)
} else {
  attach("~/Projects/Attie/GoldStandard/FinalData/aligned_geno_with_pmap.RData")
  attach("~/Projects/Attie/GoldStandard/FinalData/lipomics_final_rev2.RData")
  phe <- "INSULIN (ng/ml) 10 wk"

  ids <- findCommonID(f2g$pheno$MouseNum, lipomics$MouseNum)

  f2g <- f2g[,ids$first]
  f2g$pheno$insulin <- log10(lipomics[ids$second,phe])

  f2g <- calc.genoprob(f2g, step=1, err=0.002, map.function="c-f")

  sex <- as.numeric(f2g$pheno$Sex)-1
  out <- scanone(f2g, phe="insulin", method="hk", addcovar=sex)

  operm <- scanone(f2g, phe="insulin", method="hk", addcovar=sex, n.perm=1000, n.cluster=24)

  g <- pull.geno(fill.geno(f2g[2,], err=0.002, map.function="c-f"))[,rownames(max(out))]
  gnames <- getgenonames(class(f2g)[1], "A", cross.attr=attributes(f2g))
  g <- factor(gnames[g], levels=gnames)
  y <- f2g$pheno$insulin

  me <- tapply(y, g, mean, na.rm=TRUE)
  ci <- matrix(unlist(tapply(y, g, function(a) t.test(a)$conf.int)), nrow=2)
  dimnames(ci) <- list(c("lo", "hi"),gnames)

  save(out, g, y, me, ci, operm, file=file)
}

par(mar=c(4.1,4.1,0.1,1.1),yaxs="i")
plot(out,ylim=c(0,8.5), chr=2, ylab="LOD score", col="black")
u <- par("usr")

li <- lodint(out,"2")
segments(li[2,2],li[2,3],40,li[2,3],lwd=2,lty=2, col=color[2])
segments(li[1,2],li[1,3],40,li[1,3],lwd=2,lty=2, col=color[2])
arrows(40,li[2,3]-0.2,40,li[1,3]+0.2,len=0.1,lwd=2, col=color[2])
text(42,mean(c(li[2,3],li[1,3])),"1.5", adj=c(0,0.5), col=color[2])

u <- par("usr")
segments(li[1,2],u[3]+diff(u[3:4])*0.09,
         li[3,2],u[3]+diff(u[3:4])*0.09,lwd=2, col=color[1])
for(i in c(1,3))
  segments(li[i,2],u[3]+diff(u[3:4])*0.07,
           li[i,2],u[3]+diff(u[3:4])*0.11,lwd=2, col=color[1])
text(36,u[3]+diff(u[3:4])*0.09, "1.5-LOD support interval", adj=c(0,0.5),
     col=color[1])


dev.off()
