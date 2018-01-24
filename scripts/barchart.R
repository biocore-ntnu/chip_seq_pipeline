library(ggplot2)

df = read.table(snakemake@input[[1]], sep=" ", header=1)

p = ggplot(df, aes(x=Label, y=Counts, fill=Region)) + geom_bar(stat="identity")

pdf(snakemake@output[[1]])
plot(p)
dev.off()
