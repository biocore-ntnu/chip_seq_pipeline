library(ChIPseeker)

files = snakemake@params[["files"]]
sort_order = snakemake@params[["sort_order"]]
print(sort_order)
files_sorted = files[sort_order]
print(files_sorted)

annos = lapply(files_sorted, readRDS)
print(names(annos))

pdf(snakemake@output[[1]])
plotAnnoBar(annos)
dev.off()
