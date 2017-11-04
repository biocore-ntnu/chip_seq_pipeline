data = readRDS(snakemake@input[["data"]]) # the counts from csaw-count
design = read.table(snakemake@input[["design"]])
normfacs = unlist(read.table(snakemake@input[["normfacs"]])) # the normalization factors from csaw normfacs

library(csaw)
library(edgeR)
library(IRanges)

write("Filtering out uninteresting regions", stderr())

keep <- aveLogCPM(asDGEList(data)) >= -1
data <- data[keep,]

write("Identify differentially bound windows", stderr())

y <- asDGEList(data, norm.factors=normfacs)
y <- estimateDisp(y, design)

saveRDS(y, snakemake@output[["disp"]])
saveRDS(data, snakemake@output[["data"]])
