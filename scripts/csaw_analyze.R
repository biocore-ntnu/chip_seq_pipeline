
require(edgeR)

data = readRDS(snakemake@input[["data"]])
normfacs = readRDS(snakemake@input[["normfacs"]])
design = read.table(snakemake@input[["design"]])

# filter out uninteresting regions

keep <- aveLogCPM(asDGEList(data)) >= -1
data <- data[keep,]

# identify db windows

y <- asDGEList(data, norm.factors=normfacs)
y <- estimateDisp(y, design)

fit <-  glmLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)

# correct for multiple testing

merged <-  mergeWindows(rowranges(data), tol=1000L)
tabcom <-  combineTests(merged$id, results$table)

write.table(tabcom, snakemake@output[[1]])
