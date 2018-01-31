
import pandas as pd


def remove_all_but_one_transcript_per_gene(df):

    transcripts_with_most_exons = []

    for gene, gdf in df.groupby("Name"):

        max_exon = gdf.ExonNumber.max()
        max_transcript = gdf.loc[gdf.ExonNumber == max_exon].Transcript.iloc[0]

        max_rows = gdf.loc[gdf.Transcript == max_transcript]

        transcripts_with_most_exons.append(max_rows)

    return pd.concat(transcripts_with_most_exons).reset_index(drop=True)


if __name__ == "__main__":

    df = pd.read_table(snakemake.input[0], sep="\s+")

    result = remove_all_but_one_transcript_per_gene(df)

    result.to_csv(snakemake.output[0], sep="\t", index=False, header=True)
