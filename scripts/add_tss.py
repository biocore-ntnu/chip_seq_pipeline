import pandas as pd

def add_tss(df, barchart_tss_length):

    pd.options.mode.chained_assignment = None
    tss_pos = df.loc[(df.Strand == "+") & (df.Type == "gene")]

    tss_neg = df.loc[(df.Strand == "-") & (df.Type == "gene")]
    tss_neg.loc[:, "Start"] = tss_neg.End

    tss = pd.concat([tss_pos, tss_neg])
    tss.loc[:, "End"] = tss.Start
    tss.loc[:, "Type"] = "tss"

    tes_pos = df.loc[(df.Strand == "+") & (df.Type == "gene")]
    tes_pos.loc[:, "Start"] = tes_pos.End

    tes_neg = df.loc[(df.Strand == "-") & (df.Type == "gene")]

    tes = pd.concat([tes_pos, tes_neg])
    tes.loc[:,"Type"] = "tes"
    tes.loc[:, "End"] = tes.Start

    tes.loc[:, "End"] = tes.Start + barchart_tss_length
    tes.loc[:, "Start"] = tes.Start - barchart_tss_length

    tss.loc[:, "End"] = tss.Start + barchart_tss_length
    tss.loc[:, "Start"] = tss.Start - barchart_tss_length

    outdf = pd.concat([df, tss, tes])
    outdf = outdf.sort_values("Chromosome")

    pd.options.mode.chained_assignment = "warn"

    return outdf.reset_index(drop=True)


if __name__ == "__main__":

    df = pd.read_table(input[0], sep="\t", header=None)

    barchart_tss_length = snakemake.config["barchart_tss_length"]

    outdf = add_tss(df, barchart_tss_length)

    outdf.to_csv(output[0], sep="\t", index=False)
