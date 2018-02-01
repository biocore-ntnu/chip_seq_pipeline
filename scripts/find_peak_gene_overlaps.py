from collections import defaultdict
import pickle

import pandas as pd
from numpy import int64

from intervaltree import IntervalTree
from itertools import groupby

# def create_intervaltrees(genes):

#     genome = dict()

#     for chromosome, lines in groupby(open(genes), lambda l: l.split()[0]):
#         chromosome_intervaltree = IntervalTree()
#         for line in lines:
#             print(line)
#             start, end, name, region_type = line.split()[1:5]
#             start, end = int(start), int(end)
#             chromosome_intervaltree[start:end] = (start, name, region_type)
#             raise

#         genome[chromosome] = chromosome_intervaltree

#     return genome


def find_peak_gene_overlaps(intervaltrees, peaks):


    rowdicts = []
    for i, (chromosome, start, end, *_) in peaks.iterrows():
        gene_regions = [i[2] for i in intervaltrees[chromosome][start:end]]
        for _, _, gene_region in gene_regions:
            rowdict = {"Peak": i, "Chromosome": chromosome, "Start": start, "End":
                       end, "Region": gene_region}
            rowdicts.append(rowdict)

    df = pd.DataFrame.from_dict(rowdicts)["Chromosome Start End Peak Region".split()]
    df.index = df.index.astype(int64)

    return df.sort_values("Peak Region".split()).reset_index(drop=True)


def parse_overlap_dataframe(df):

    df.loc[:, "Region"] = df.Region.str.replace("gene", "intron")
    gene_regions_prioritized = "tss tes exon intron".split()
    counts = defaultdict(int)
    for g, gdf in df.groupby("Peak"):
        for gene_region in gene_regions_prioritized:
            if gene_region in list(gdf.Region):
                counts[gene_region] += 1
                break

    outdf = pd.DataFrame.from_dict(counts, orient="index").reset_index()
    outdf.columns = ["Region", "Counts"]

    return outdf.sort_values("Region").reset_index(drop=True)


def create_barchart_data(its, peak_file, label):

    peaks = pd.read_table(peak_file, sep="\s+", header=0, comment="#")

    df = find_peak_gene_overlaps(its, peaks)
    count_df = parse_overlap_dataframe(df)
    count_df.insert(2, "Label", label)

    return count_df.reset_index(drop=True)


if __name__ == "__main__":

    group = snakemake.params.label
    peaks = snakemake.input.peaks
    intervals = snakemake.input.intervals

    interval_trees = pickle.load(open(intervals, "rb"))

    count_df = create_barchart_data(interval_trees, peaks, group)

    count_df.to_csv(snakemake.output[0], sep=" ", index=False, header=True)
