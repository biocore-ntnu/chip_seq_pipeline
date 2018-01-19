import sys
import pandas as pd


def compute_internal_exons(exons):
    """Avoid first and last exon for each transcript."""

    colnames = "Chromosome Start End Name Score Strand Gene Transcript".split()

    info = exons.Name.str.split(":|\.", expand=True).iloc[:,[1, 2, 3]]
    info.columns = "name transcript exon".split()
    info.exon = info.exon.astype(int)

    exons = pd.concat([exons, info], axis=1).sort_values(["transcript", "exon"])

    exons = exons.groupby(["name", "transcript"]).apply(lambda r: r.iloc[1:-1])

    exons = exons.rename(index=str, columns={"name": "Gene", "transcript": "Transcript"})

    exons.loc[:, "Transcript"] = exons.Transcript.astype(int)

    return exons[colnames].reset_index(drop=True)



if __name__ == "__main__":
    exons = pd.read_table(sys.argv[1], header=None,
                          names="Chromosome Start End Name Score Strand Gene".split(), index_col=False)

    outdf = compute_internal_exons(exons)

    outdf.to_csv(sys.stdout, sep="\t", header=False, index=False)
