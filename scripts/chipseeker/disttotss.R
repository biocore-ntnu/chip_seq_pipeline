library(ChIPseeker)

files = snakemake@params[["files"]]
sort_order = snakemake@params[["sort_order"]]
files_sorted = files[sort_order]

annos = lapply(files_sorted, readRDS)

pdf(snakemake@output[[1]])
plotDistToTSS(annos)
dev.off()
