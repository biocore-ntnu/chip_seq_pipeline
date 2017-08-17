groups = snakemake@params[["groups"]]
group = factor(groups)

print(group)
design <- model.matrix(~0+group)

colnames(design) <- gsub("group", "", colnames(design))

write.table(design, snakemake@output[[1]])
