design = snakemake@input[["design"]]
evalues = snakemake@input[["e_values"]]
weights = snakemake@input[["weights"]]
contrast = snakemake@wildcards[["contrast"]]
outfile.all = snakemake@output[["all_regions"]]
outfile.cutoff = snakemake@output[["cutoff"]]

library(edgeR)

design = read.table(design, sep=" ", header=1, row.names=1)

df = read.table(evalues, sep=" ", header=1, row.names=1)
weights = read.table(weights, sep=" ", header=1, row.names=1)

lm = lmFit(as.matrix(df), design, weigths=weights)

m = makeContrasts(contrasts=contrast, levels=design)
cf = contrasts.fit(lm, contrasts=m)

fit = eBayes(cf)

tta = topTable(fit, sort.by="P", n=Inf, coef=contrast)
write.table(tta, outfile.all, sep=" ")

ttc = topTable(fit, sort.by="P", n=Inf, coef=contrast, p.value=0.05)
write.table(ttc, outfile.cutoff, sep=" ")
