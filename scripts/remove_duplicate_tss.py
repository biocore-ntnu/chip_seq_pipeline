
import pandas as pd

def remove_duplicate_tss(df):

    tdf = df.loc[df.Type == "gene"]

    indexes_to_keep = []
    for _, cdf in tdf.groupby("Chromosome"):
        forward = tdf.loc[tdf.Strand == "+"]
        reverse = tdf.loc[tdf.Strand == "-"]

        new_tdf = pd.concat([forward.Start, reverse.End])
        df_ix = new_tdf.drop_duplicates()
        indexes_to_keep.append(df_ix)

    ix = pd.concat(indexes_to_keep).index
    tdf = tdf.ix[ix]

    df_no_overlapping_tss = df.loc[df.Transcript.isin(tdf.Transcript)]

    return df_no_overlapping_tss.reset_index(drop=True)




if __name__ == "__main__":

    df = pd.read_table(snakemake.input[0], header=0)

    df.to_csv(snakemake.output[0], sep="\t", index=False, header=False)
