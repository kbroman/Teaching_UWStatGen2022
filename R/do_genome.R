# reconstruction of a DO genome

library(qtl2)
library(here)

load(here("Data/do359_genome.RData"))

pdf(here("Figs/do_genome.pdf"), width=9.75, height=5, pointsize=16)
par(mar=c(3.1,4.1,0.6,0.6))
plot_onegeno(ph, pmap, mgp.x=c(1.8,0.3,0),
             ylab="Position (Mbp)")
legend(9.5, 163, names(CCcolors), pt.bg=CCcolors, pch=22, bg="gray90", ncol=4)
dev.off()
