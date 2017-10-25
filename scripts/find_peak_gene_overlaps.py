
from bx.intervals.intersection import IntervalTree

genes = snakemake.input.genes
peaks = snakemake.input.peaks

# distance = snakemake.params

genome = dict()

for chromosome, lines in groupby(open(genes), lambda l: l.split()[0]):
    chromosome_intervaltree = IntervalTree()
    for line in lines:
        start, _, name, quartile = line.split()[1:5]
        start = int(start)
        chromosome_intervaltree.add(start - 5000, start + 5000, (start, name, quartile))

    genome[chromosome] = chromosome_intervaltree
