library(edgeR)

df = read.table(gzfile(snakemake@input[["counts"]]), header=T, row.names=1)
design = read.table(snakemake@input[["design"]], header=T, row.names=1)

dge = DGEList(counts=as.matrix(df))
dge <- calcNormFactors(dge, method=snakemake@params["voom_normalization_dge"])

pdf(snakemake@output[["plot"]])
y = voom(dge, design, plot=TRUE, normalize.method=snakemake@params[["voom_normalization"]])
dev.off()

write.table(y$E, snakemake@output[["e_values"]], sep=" ")
write.table(y$weights, snakemake@output[["weights"]], sep=" ")
