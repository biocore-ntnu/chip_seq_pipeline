library(edgeR)

df = read.table(gzfile(snakemake@input[["counts"]]), header=T, row.names=1, sep=" ")
## design = read.table(snakemake@input[["design"]], header=T, row.names=1)

print(snakemake@params[["voom_normalization_dge"]])

dge = DGEList(counts=as.matrix(df))
dge <- calcNormFactors(dge, method=snakemake@params[["voom_normalization_dge"]])

pdf(snakemake@output[["plot"]])
y = voom(dge, plot=TRUE) # , normalize.method=snakemake@params[["voom_normalization"]]) #design,
dev.off()

write.table(y$E, snakemake@output[["e_values"]], sep=" ")
write.table(y$weights, snakemake@output[["weights"]], sep=" ")
