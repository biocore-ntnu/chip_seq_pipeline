import re

from os import environ
from os.path import basename

import pandas as pd

def merge_colnames_sample_sheet(columns, ss):

    names = pd.Series(columns, name="FullName")
    basenames = pd.Series([basename(f).split(".")[0] for f in columns], name="BaseName")

    names_basenames = pd.concat([names, basenames], 1)
    m = names_basenames.merge(ss, left_on="BaseName", right_on="Name").drop("BaseName", axis=1)

    return m

def error_if_not_using_tmux():

    if not environ.get("TMUX", ""):
        raise Exception("Not using TMUX!")

def get_if_not_empty(config, key, default):
    "like dict.get(key, default) only that it checks that the entry is not empty '' [] {} too."

    if key in config and config[key]:
        return config[key]
    else:
        return default

# not automatically discovered since they give errors when sample sheet only has one group
multi_group_targets = ["log2_ratio_group_vs_group_bigwig_chip_only",
                       "bigwig_log2ratio_sample_vs_group",
                       "log2_ratio_input_normalized_group_vs_input_normalized_group"]


def fetch_main_targets(snakefile="Snakefile"):
    "Get all targets in main Snakefile"

    r = re.compile("^rule (.*):")

    targets = []
    for line in open(snakefile):
        match = r.search(line)
        if match:
            targets.append(match.group(1))

    return targets

if __name__ == "__main__":

    targets = fetch_main_targets()

    print(targets)
