library(csaw)

files = snakemake@input[["chip_input"]]

param <- readParam(minq=snakemake@config[["fastq_quality"]])
binned <- windowCounts(files, bin=TRUE, width=10000, param=param)
normfacs <-  normOffsets(binned)

n = as.data.frame(normfacs)
rownames(n) = files

write.table(n, snakemake@output[[1]])
