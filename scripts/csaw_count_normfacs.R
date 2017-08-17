
library(csaw)

param <- readParam(minq=snakemake@config[["fastq_quality"]])
binned <- windowCounts(snakemake@config[["filenames_as_string"]], bin=TRUE, width=10000, param=param)
normfacs <-  normOffsets(binned)

saveRDS(snakemake@output[[1]])
