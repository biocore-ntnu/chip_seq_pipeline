library(limma)

sample.sheet = read.table(snakemake@input[[1]], col.names=1)
ss = sample.sheet

ssc = ss[ss$ChIP == "ChIP",]
groups = sort(subset(ssc, !duplicated(ssc$Name))$Group)

name = factor(groups, levels=unique(groups))

design = model.matrix(~0 + name)

## makeContrasts

write.table(design, snakemake@output[[1]], sep=" ")
