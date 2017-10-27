from itertools import groupby
from collections import defaultdict

import pandas as pd

from bx.intervals.intersection import IntervalTree


def create_intervaltrees(genes):

    genome = dict()

    for chromosome, lines in groupby(open(genes), lambda l: l.split()[0]):
        chromosome_intervaltree = IntervalTree()
        for line in lines:
            start, end, name, region_type = line.split()[1:5]
            start, end = int(start), int(end)
            chromosome_intervaltree.add(start, end, (start, name, region_type))

        genome[chromosome] = chromosome_intervaltree

    return genome


def find_peak_gene_overlaps(intervaltrees, peaks):

    rowdicts = []
    for i, (chromosome, start, end) in peaks.iterrows():
        gene_regions = intervaltrees[chromosome].find(start, end)
        for _, _, gene_region in gene_regions:
            rowdict = {"Peak": i, "Chromosome": chromosome, "Start": start, "End":
                    end, "Region": gene_region}
            rowdicts.append(rowdict)

    return pd.DataFrame.from_dict(rowdicts)["Chromosome Start End Peak Region".split()]


def parse_overlap_dataframe(df):

    df.loc[:, "Region"] = df.Region.str.replace("gene", "intergenic")
    gene_regions_prioritized = "tss tes exon intergenic".split()
    counts = defaultdict(int)
    for g, gdf in df.groupby("Peak"):
        for gene_region in gene_regions_prioritized:
            if gene_region in list(gdf.Region):
                counts[gene_region] += 1
                break

    outdf = pd.DataFrame.from_dict(counts, orient="index").reset_index()
    outdf.columns = ["Region", "Counts"]

    return outdf.sort_values("Region")


def create_barchart_data(genes, peak_file, label):

    its = create_intervaltrees(genes)

    peaks = pd.read_table(peak_file, sep="\s+")

    df = find_peak_gene_overlaps(its, peaks)
    count_df = parse_overlap_dataframe(df)
    count_df.insert(2, "Label", label)

    return count_df


if __name__ == "__main__":

    group = snakemake.params.label
    genes = snakemake.input.genes

    count_df = create_barchart_data(genes, peak_files_dict, group)

    count_df.to_csv(snakemake.output[0], sep=" ", index=False, header=True)
