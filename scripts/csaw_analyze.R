data = readRDS(snakemake@input[["data"]])
design = read.table(snakemake@input[["design"]])
normfacs = unlist(read.table(snakemake@input[["normfacs"]]))
y = readRDS(snakemake@input[["estimate_disp"]])


library(csaw)
library(edgeR)
library(IRanges)

contrast = as.numeric(makeContrasts(snakemake@wildcards[["contrast"]], levels=colnames(design)))
fit <-  glmQLFit(y, design, contrast=contrast, robust=TRUE)
results <- glmQLFTest(fit, contrast=contrast)

# correct for multiple testing

print("correct for multiple testing")

merged <- mergeWindows(rowRanges(data), tol=1000L)
tabcom <- combineTests(merged$id, results$table)

print("writing output file")

df = as.data.frame(merged$region)
row.names = with(df, paste0(seqnames, "_", start, "_", end))

rownames(tabcom) = row.names

write.table(tabcom, snakemake@output[[1]])
