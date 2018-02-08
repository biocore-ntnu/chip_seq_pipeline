import pandas as pd

def extract_correct_gene_type(df, gene_type):

    df = df.loc[df.Type == gene_type]

    return df


if __name__ == "__main__":

    df = pd.read_table(snakemake.input[0])

    gene_type = snakemake.wildcards["gene_type"]

    df = extract_correct_gene_type(df, gene_type)

    df.to_csv(snakemake.output[0], sep="\t", index=False)
