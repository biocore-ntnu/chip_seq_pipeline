library(tidyverse)
library(ggplot2)

t = read_delim(snakemake@input[[1]], delim=" ", col_names=TRUE)
p = ggplot(t, aes(x=Sample, y=Ratio)) + geom_violin(scale="count") + facet_wrap(LeftOut ~ Contrast) + geom_boxplot(width=.1, fill="black", outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5)
ggsave(snakemake@output[[1]], p)
