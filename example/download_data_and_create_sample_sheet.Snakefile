import pandas as pd


from os.path import abspath


ss = pd.read_table("example_data_links.txt", sep=" ")
print(ss.ChIP.str.lower())
print(ss)


# rule download_example_data:
#     output:
#         "data/{chip}_{sample}_{type}.bed.gz"
#     params:
#         url = lambda w: ss.loc[(ss.CellType == w.type) & (ss.Sample == int(w.sample)) & (ss.ChIP == w.chip)].URL.iloc[0]
#     shell:
#         "curl {params.url} > {output[0]}"


infiles = {(chip, sample, celltype): "data/{chip}_{sample}_{type}.bed.gz".format(sample=sample, type=celltype, chip=chip) for chip, sample, celltype in zip(ss.ChIP, ss.Sample, ss.CellType)}
print(infiles)


rule create_sample_sheet:
    input:
        infiles.values()
    output:
        "sample_sheet.txt"
    run:
        rowdicts = []

        for (chip, sample, celltype), f in infiles.items():
            rowdicts.append({"Name": celltype + "_" + str(sample) + "_" + chip, "File": abspath(f), "Group": celltype, "ChIP": chip, "Mate": 1})

        df = pd.DataFrame.from_dict(rowdicts)["File Name Group ChIP Mate".split()]

        df.to_csv(output[0], sep=" ", index=False)
