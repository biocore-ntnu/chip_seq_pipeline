import os
import re
import pandas as pd

from snakemake.io import expand


def read_sample_sheet(sample_sheet):

    return pd.read_table(sample_sheet, sep="\s+")


def sample_vs_group(ss, chip):

    ss = ss.loc[ss.ChIP == chip]

    samples_unique = ss.drop_duplicates("Name")["Name Group".split()]
    d = samples_unique.set_index("Name").to_dict()["Group"]

    groups = set(ss.Group)

    rowdicts = []
    for sample, group in d.items():
        groups_copy = set(list(groups))
        groups_copy.remove(group)
        for g in groups_copy:
            rowdict = {"OtherGroup": g, "Sample": sample}
            rowdicts.append(rowdict)

    df = pd.DataFrame.from_dict(rowdicts)
    return df["Sample OtherGroup".split()].sort_values("Sample OtherGroup".split()).reset_index(drop=True)


def read_rna_seq_sample_sheet(sample_sheet):

    if os.path.exists:
        return read_sample_sheet(sample_sheet)


def _get_samples(sample_sheet, chip, group=None):

    chip_idx = sample_sheet.ChIP.str.contains(chip, flags=re.IGNORECASE)
    if not group:
        samples = sample_sheet.loc[chip_idx].Name.drop_duplicates()
    else:
        group_idx = sample_sheet.Group == group
        samples = sample_sheet.loc[group_idx & chip_idx].Name.drop_duplicates()

    return list(samples)


def find_filetype(ss):

    filetypes = ss.File.str.replace(".gz$", "").str.split(".", expand=True).iloc[:, -1]
    filetypes = list(filetypes.drop_duplicates())
    assert len(filetypes) == 1, "More than one filetype in sample sheet: " + ", ".join(filetypes)

    filetype = filetypes[0]

    assert filetype in ["bed", "bam", "fastq", "fq"]

    return filetype


def correct_cs_files(sample_sheet, prefix, chip, extension, config, group):
    "Gets the chip seq files. By default gets all files."

    samples = _get_samples(sample_sheet, chip, group)

    filetype = find_filetype(sample_sheet)

    if filetype == extension:
        return list(sample_sheet.loc[sample_sheet.Name.isin(samples)].File)

    if extension == "bam" :
        fs = "{prefix}/data/bam/{sample}.sorted.bam"
    elif config["paired_end"]:
        fs = "{prefix}/data/bam/{sample}.bedpe"
    else:
        fs = "{prefix}/data/bam/{sample}.bed"

    files = expand(fs, sample=samples, extension=extension, prefix=prefix)

    return files
