
library(csaw)

param = readParam(minq=snakemake@config[["fastq_quality"]])
data = windowCounts(snakemake@config[["filenames_as_string"]], ext=snakemake@config[["fragment_size"]]/2, width=snakemake@config[["window_size"]], spacing=snakemake@config[["window_size"]], param=param, filter=snakemake@config[["csaw_filter"]])

regions = as.data.frame(rowRanges(data))
rownames = with(regions, paste(seqnames, start, end, sep="_"))

counts = as.data.frame(assay(data))

write.table(counts, snakemake@output[["table"]], sep=" ", row.names=rownames, quote=F)

saveRDS(counts, snakemake@output[["rdata"]])
