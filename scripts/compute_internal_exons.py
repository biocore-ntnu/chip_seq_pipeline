import sys
import pandas as pd


def compute_internal_exons(exons):
    """Avoid first and last exon for each transcript."""

    colnames = "Chromosome Start End Name Score Strand Gene".split() #  Transcript

    # print(exons)

    info = exons.Name.str.split(":|\.", expand=True).iloc[:,[1, 2, 3]]
    info.columns = "name transcript exon".split()
    info.exon = info.exon.astype(int)

    exons = pd.concat([exons, info], axis=1).sort_values(["transcript", "exon"])

    internal_exons_per_gene = []
    for gene, gdf in exons.groupby(["Gene", "transcript"]):
        r = gdf.loc[~gdf.exon.isin([gdf.exon.min(), gdf.exon.max()])]
        internal_exons_per_gene.append(r)

    exons = pd.concat(internal_exons_per_gene)

    return exons[colnames].reset_index(drop=True)



if __name__ == "__main__":
    exons = pd.read_table(sys.argv[1], header=None,
                          names="Chromosome Start End Name Score Strand Gene".split(), index_col=False)

    outdf = compute_internal_exons(exons)

    outdf.to_csv(sys.stdout, sep="\t", header=False, index=False)
