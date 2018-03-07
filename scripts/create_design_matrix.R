sample.sheet = read.table(snakemake@input[[1]], header=TRUE)
ss = sample.sheet

ssc = ss[ss$ChIP == "ChIP",]
groups = sort(subset(ssc, !duplicated(ssc$Name))$Group)

names = unique(ssc[with(ssc, order(Group)), ]$Name)

group.name = factor(groups, levels=unique(groups))

design = model.matrix(~0 + group.name)

colnames(design) <- gsub("group.name", "", colnames(design))
rownames(design) <- names

write.table(design, snakemake@output[[1]], sep=" ")
