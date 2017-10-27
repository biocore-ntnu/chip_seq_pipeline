library(ggplot2)

df = read.table(snakemake@input[[1]], sep=" ")

p = ggplot(df, aes(x=Label, y=Counts, fill=Region)) + geom_bar(stat="identity")

ggsave(snakemake@output[[1]], p)
