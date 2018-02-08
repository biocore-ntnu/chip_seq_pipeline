library(ggplot2)

df = read.table(snakemake@input[[1]], sep=" ", header=1)

names = snakemake@params[["sort_order"]]

df$Label = factor(df$Label, levels=names)

p = ggplot(df, aes(x=Label, y=Counts, fill=Region)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.5))

pdf(snakemake@output[[1]])
plot(p)
dev.off()
