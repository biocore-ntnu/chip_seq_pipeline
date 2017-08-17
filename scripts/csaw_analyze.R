
data = readRDS(snakemake@input[["data"]])
design = read.table(snakemake@input[["design"]])
normfacs = unlist(read.table(snakemake@input[["normfacs"]]))

require(csaw)
require(edgeR)

# filter out uninteresting regions

print("filtering")

keep <- aveLogCPM(asDGEList(data)) >= -1
data <- data[keep,]

# identify db windows

print("identify db windows")
y <- asDGEList(data, norm.factors=normfacs)
y <- estimateDisp(y, design)

fit <-  glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)

# correct for multiple testing


print("correct for multiple testing")

merged <-  mergeWindows(rowranges(data), tol=1000L)
tabcom <-  combineTests(merged$id, results$table)

print("writing output file")

write.table(tabcom, snakemake@output[[1]])
