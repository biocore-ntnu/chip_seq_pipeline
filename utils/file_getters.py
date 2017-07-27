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

def _get_samples(sample_sheet, group, chip):
    group_idx = sample_sheet.Group == group
    chip_idx = sample_sheet.ChIP.str.contains(chip, flags=re.IGNORECASE)
    samples = sample_sheet.loc[group_idx & chip_idx].Name.drop_duplicates()

    return list(samples)


def correct_cs_files(sample_sheet, group, prefix, chip, extension, config):

    samples = _get_samples(sample_sheet, group, chip)

    if extension == "bam" and not config["bam"]:
        fs = "{prefix}/data/align/{sample}.bam"
    elif extension == "bam" and config["bam"]:
        ss = sample_sheet
        files = list(ss.loc[ss.Name.isin(samples)].File)
        fs = files
    elif config["paired_end"]:
        fs = "{prefix}/data/align/{sample}.bedpe"
    else:
        fs = "{prefix}/data/align/{sample}.bed"

    return expand(fs, sample=samples, extension=extension, prefix=prefix)



def get_bam_or_bed(sample, no_multi_aligners, file_format):

    f = "{prefix}/data/align/{sample}.{file_format}"

    return f.format(**vars())
