import sys
import pandas as pd


def compute_internal_exons(df):
    """Avoid first and last exon for each transcript."""

    internal_exons = []
    exons = df.loc[df.Type == "exon"]
    for transcript, tdf in exons.groupby("Transcript"):

        r = tdf.sort_values("ExonNumber").iloc[1:-1]

        internal_exons.append(r)

    exons = pd.concat(internal_exons)

    return exons.reset_index(drop=True)



if __name__ == "__main__":

    exons = pd.read_table(snakemake.input[0], header=0, index_col=False)

    outdf = compute_internal_exons(exons)

    outdf.to_csv(snakemake.output[0], sep="\t", header=0, index=False)
