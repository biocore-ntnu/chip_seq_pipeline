library(edgeR)

design = read.table(snakemake@input[["design"]], sep=" ", header=1, row.names=1)

df = read.table(snakemake@input[["e_values"]], sep=" ", header=1, row.names=1)
weights = read.table(snakemake@input[["weigths"]], sep=" ", header=1, row.names=1)

lm = lmFit(as.matrix(df), design, weigths=weights)
fit = eBayes(lm)

tt = topTable(fit, sort.by="none", n=Inf, coef=c(snakemake@wildcards[["coefficient"]]))
write.table(tt, snakemake@output[[1]], sep=" ")
