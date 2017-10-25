from itertools import groupby

from bx.intervals.intersection import IntervalTree

# genes = snakemake.input.genes
# peaks = snakemake.input.peaks

# distance = snakemake.params


def create_intervaltree(genes):

    genome = dict()

    for chromosome, lines in groupby(open(genes), lambda l: l.split()[0]):
        chromosome_intervaltree = IntervalTree()
        for line in lines:
            start, end, name, region_type = line.split()[1:5]
            start, end = int(start), int(end)
            chromosome_intervaltree.add(start, end, (start, name, region_type))

        genome[chromosome] = chromosome_intervaltree

    return genome


def find_gene_peak_overlaps(intervaltrees, peaks):



    pass


        # # check whether read hits inside tss region
        # quartile_counts = defaultdict(int)
        # name_counts = defaultdict(int)
        # for line in open(input.data):

        #     chromosome, start, end, _, _, strand = line.split()
        #     start, end = int(start), int(end)

        #     # extend the read by fragment length
        #     if strand == "-":
        #         read_midpoint = end - 75
        #     else:
        #         read_midpoint = start + 75

        #     for tss, name, quartile in genome[chromosome].find(read_midpoint, read_midpoint + 1):
        #         tss_distance = tss - read_midpoint
        #         tss_distance_bin = tss_distance - (tss_distance % 50)
        #         quartile_counts[quartile, tss_distance_bin] += 1
        #         name_counts[name, tss_distance_bin] += 1
