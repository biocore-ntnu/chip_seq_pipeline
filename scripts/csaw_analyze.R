# data = readRDS("projects/barbara/data/peaks/csaw/count/count_table.rds")
# design = read.table("projects/barbara/data/peaks/csaw/design.csv")
# normfacs = unlist(read.table("projects/barbara/data/peaks/csaw/normfacs/normfacs.csv"))

data = readRDS(snakemake@input[["data"]])
design = read.table(snakemake@input[["design"]])
normfacs = unlist(read.table(snakemake@input[["normfacs"]]))

library(csaw)
library(edgeR)
library(IRanges)

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

merged <-  mergeWindows(rowRanges(data), tol=1000L)
tabcom <-  combineTests(merged$id, results$table)

print("writing output file")

df = as.data.frame(merged$region)
row.names = with(df, paste0(seqnames, "_", start, "_", end))

rownames(tabcom) = row.names

write.table(tabcom, snakemake@output[[1]])
