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

saveRDS(y, snakemake@output[["disp"]])
saveRDS(data, snakemake@output[["data"]])
