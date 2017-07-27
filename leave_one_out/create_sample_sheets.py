import pandas as pd

def create_sample_sheet(sample_sheet):
    ss = sample_sheet

    dfs = []
    for group, gdf in ss.groupby("Group"):

        samples = list(gdf.loc[gdf.ChIP != "Input"].Name.drop_duplicates())

        for sample in samples:
            ndf = gdf.loc[gdf.Name != sample].reset_index(drop=True)
            ndf.Group = sample + "_lo"
            dfs.append(ndf)

    return pd.concat(dfs).reset_index(drop=True)


def create_config_files():
    pass
