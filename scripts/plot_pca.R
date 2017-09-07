library(ggplot2)

f = snakemake@input[[1]]
o = snakemake@output[[1]]

df = read.table(f, header=1)
p = ggplot(df, aes(x=X, y=Y, colour=Group)) + geom_text(aes(label=File))

ggsave(o, p)
