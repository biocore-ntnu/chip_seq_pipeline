library(ggplot2)

coords = snakemake@input[["coords"]]
variance = snakemake@input[["variance_explained"]]
o = snakemake@output[[1]]

df = read.table(coords, header=1)
var = read.table(variance, header=1)
p = ggplot(df, aes(x=X, y=Y, colour=Group)) + geom_text(aes(label=File)) + labs(x = var$X, y = var$Y)

ggsave(o, p)
