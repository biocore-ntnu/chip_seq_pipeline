import pickle
from itertools import groupby
from collections import defaultdict
from intervaltree import IntervalTree

def create_intervaltrees(genes):

    genome = dict()

    file_handle = open(genes)
    next(file_handle) # skip header

    for chromosome, lines in groupby(file_handle, lambda l: l.split()[0]):
        chromosome_intervaltree = IntervalTree()
        for line in lines:
            start, end, region_type, _, name = line.split()[1:6]
            start, end = int(start), int(end)
            chromosome_intervaltree[start:end] = (start, name, region_type)

        genome[chromosome] = chromosome_intervaltree

    return genome


if __name__ == "__main__":

    genes = snakemake.input[0]

    its = create_intervaltrees(genes)

    pickle.dump(its, open(snakemake.output[0], "wb"))
