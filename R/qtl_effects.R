library(here)
library(broman)
library(qtl2)

# locus on chr 1 (rs235841505, A,F,G vs others)
# locus on chr 7 (rs32459628, A vs others)
# locus on chr 5 (SV_5_66529849_66529851, A,C,F
# locus on chr 17 (rs50608335, F,G)

# all 36 effects; 8 alleles; 3 snp genotypes

chr1_max <- structure(list(chr = "1", pos = 129.809947, pheno1 = 10.018838147179), class = "data.frame", row.names = "backupJAX00266606")
chr7_max <- structure(list(chr = "7", pos = 20.027566, pheno1 = 7.30785166348042), class = "data.frame", row.names = "UNC070931469")
chr5_max <- structure(list(chr = "5", pos = 71.4184343333333, pheno1 = 17.3450713492845), class = "data.frame", row.names = "c5.loc38")
chr17_max <- structure(list(chr = "17", pos = 83.3622237692308, pheno1 = 17.4338652119796), class = "data.frame", row.names = "c17.loc52.7")


out_snps <- readRDS(here("QTLresults/out_snps.rds"))
do <- readRDS(here("QTLresults/do.rds"))
pr_sub <- readRDS(here("QTLresults/pr_sub.rds"))
apr_sub <- readRDS(here("QTLresults/apr_sub.rds"))
kinship <- readRDS(here("QTLresults/kinship.rds"))

phe <- log10(do$pheno[,2])
covar <- cbind(sex=(do$is_female*1),
               wbc=log10(do$pheno[,1]))


fit1_f <- lapply(seq_along(pr_sub), function(i)
    fit1(pr_sub[[i]], phe, addcovar=covar, kinship=kinship[[c("1","7","5","17")[i]]]))
fit1_a <- lapply(seq_along(pr_sub), function(i)
    fit1(apr_sub[[i]], phe, addcovar=covar, kinship=kinship[[c("1","7","5","17")[i]]]))

groups <- list(rs235841505=c("B","F","G"),
               rs32459628="A",
               SV_5_66529849_66529851=c("A","C","F"),
               rs50608335=c("F","G"))


pr_to_snp <-
    function(p, alleles=c("B","F","G"))
    {
        g <- colnames(p)
        g <- strsplit(g, "")
        gn <- sapply(g, function(x) sum(x %in% alleles))
        x <- matrix(0,ncol=3, nrow=ncol(p))
        for(i in 1:3) x[gn==i-1,i] <- 1
        p %*% x
    }

snp_pr_sub <- lapply(seq_along(pr_sub), function(i)
                     pr_to_snp(pr_sub[[i]], groups[[i]]))
names(snp_pr_sub) <- names(groups)

fit1_snp <- lapply(seq_along(pr_sub), function(i)
    fit1(snp_pr_sub[[i]], phe, addcovar=covar, kinship=kinship[[c("1","7","5","17")[i]]]))


for(i in 1:4) {
    file <- here("Figs/qtl_effects")
    if(i > 1) file <- paste0(file, "_", i)
    file <- paste0(file, ".pdf")

    pdf(file, height=6, width=11, pointsize=14)

    layout(rbind(c(1,1),c(2,3)), width=c(2,1))
    par(mar=c(3.1,4.1,0.6,0.6))

    yl <- range(c(fit1_f[[i]]$coef[1:36]-2*fit1_f[[i]]$SE[1:36],
                  fit1_f[[i]]$coef[1:36]+2*fit1_f[[i]]$SE[1:36],
                  fit1_a[[i]]$coef[1:8]-2*fit1_a[[i]]$SE[1:8],
                  fit1_a[[i]]$coef[1:8]+2*fit1_a[[i]]$SE[1:8],
                  fit1_snp[[i]]$coef[1:3]-2*fit1_snp[[i]]$SE[1:3],
                  fit1_snp[[i]]$coef[1:3]+2*fit1_snp[[i]]$SE[1:3]))
    yl <- c(-max(abs(yl)), max(abs(yl)))


    ciplot(fit1_f[[i]]$coef[1:36], fit1_f[[i]]$SE[1:36], las=2,
           ylab="QTL effect", xlab="Genotype", ylim=yl, mgp.x=c(1.8, 0.2, 0))
    ciplot(fit1_a[[i]]$coef[1:8], fit1_a[[i]]$SE[1:8],
           ylab="QTL effect", xlab="Allele", ylim=yl, mgp.x=c(1.8, 0.2, 0))

    ciplot(fit1_snp[[i]]$coef[1:3], fit1_snp[[i]]$SE[1:3],
           labels=c("AA","AB","BB"),
           ylab="QTL effect", xlab="SNP genotype", ylim=yl,
           mgp.x=c(1.8, 0.2, 0))

    dev.off()
}
