library(qtl2)
library(here)

do <- readRDS(here("QTLresults/do.rds"))

gmap <- insert_pseudomarkers(do$gmap, stepwidth="max", step=0.2)
pmap <- interp_map(gmap, do$gmap, do$pmap)
phe <- log10(do$pheno[,2])
covar <- cbind(sex=(do$is_female*1),
               wbc=log10(do$pheno[,1]))

kinship <- readRDS(here("QTLresults/kinship.rds"))


file <- here("R/_cache/probs.rds")
if(file.exists(file)) {
#    pr <- readRDS(file)
} else {
    pr <- calc_genoprob(do, gmap, error_prob=0.002, map_function="c-f", cores=0)
    saveRDS(pr, file)
}

# calculate allele dosages
file <- here("R/_cache/aprobs.rds")
if(file.exists(file)) {
#    apr <- readRDS(file)
} else {
    apr <- genoprob_to_alleleprob(pr, cores=0)
    saveRDS(apr, file)
}

# calculate kinship_noloco
file <- here("QTLresults/kinship_noloco.rds")
if(file.exists(file)) {
    k <- readRDS(file)
} else {
    k <- calc_kinship(apr, cores=0)
    saveRDS(k, file)
}


chr1_max <- structure(list(chr = "1", pos = 129.809947, pheno1 = 10.018838147179), class = "data.frame", row.names = "backupJAX00266606")
chr7_max <- structure(list(chr = "7", pos = 20.027566, pheno1 = 7.30785166348042), class = "data.frame", row.names = "UNC070931469")
chr5_max <- structure(list(chr = "5", pos = 71.4184343333333, pheno1 = 17.3450713492845), class = "data.frame", row.names = "c5.loc38")
chr17_max <- structure(list(chr = "17", pos = 83.3622237692308, pheno1 = 17.4338652119796), class = "data.frame", row.names = "c17.loc52.7")

pr_sub <- list(chr1=pull_genoprobpos(pr, marker=rownames(chr1_max)),
               chr7=pull_genoprobpos(pr, marker=rownames(chr7_max)),
               chr5=pull_genoprobpos(pr, marker=rownames(chr5_max)),
               chr17=pull_genoprobpos(pr, marker=rownames(chr17_max)))
apr_sub <- list(chr1=pull_genoprobpos(apr, marker=rownames(chr1_max)),
                chr7=pull_genoprobpos(apr, marker=rownames(chr7_max)),
                chr5=pull_genoprobpos(apr, marker=rownames(chr5_max)),
                chr17=pull_genoprobpos(apr, marker=rownames(chr17_max)))
saveRDS(pr_sub, here("QTLresults/pr_sub.rds"))
saveRDS(apr_sub, here("QTLresults/apr_sub.rds"))
