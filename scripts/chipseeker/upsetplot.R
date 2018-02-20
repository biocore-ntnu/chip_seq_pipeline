
library(ChIPseeker)

peakAnno = readRDS(snakemake@input[[1]])

pdf(snakemake@output[[1]])
upsetplot(peakAnno)
dev.off()
