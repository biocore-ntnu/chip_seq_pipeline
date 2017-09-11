library(edgeR)

df = read.table(snakemake@input[[1]], header=T, row.names=1, sep=" ", check.names=FALSE)
df <- 2 ^ df
## design = read.table(snakemake@input[["design"]], header=T, row.names=1)

print(snakemake@params[["voom_normalization_dge"]])

dge = DGEList(counts=as.matrix(df))
dge <- calcNormFactors(dge, method=snakemake@params[["voom_normalization_dge"]])

pdf(snakemake@output[["plot"]])
y = voom(dge, lib.size=c(1e6), plot=TRUE) # , normalize.method=snakemake@params[["voom_normalization"]]) #design,
dev.off()

write.table(y$E, snakemake@output[["e_values"]], sep=" ")
write.table(y$weights, snakemake@output[["weights"]], sep=" ")
