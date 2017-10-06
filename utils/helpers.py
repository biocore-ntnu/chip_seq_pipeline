import re

from os import environ
from os.path import basename

import pandas as pd


from snakemake.io import expand

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
multi_group_targets = ["log2_ratio_sample_vs_group_bigwig",
                       "log2_ratio_group_vs_group_bigwig",
                       "log2_ratio_group_vs_group_heatmap"]


def fetch_main_targets(snakefile="Snakefile"):
    "Get all targets in main Snakefile"

    r = re.compile("^rule (.*):")

    targets = []
    for line in open(snakefile):
        match = r.search(line)
        if match:
            targets.append(match.group(1))

    return targets


def expand_zip(template, zip_dict, regular_dict):

    regular_vars = regular_dict.keys()

    new_template = template

    for k, v in regular_dict.items():
        if isinstance(v, str):
            new_template = new_template.replace("{" + k + "}", v)

    new_templates = []
    for k, v_list in regular_dict.items():

        if not isinstance(v_list, str):
            for v in v_list:
                r = "{" + k + "}"
                new_templates.append(new_template.replace(r, v))

    final_paths = []
    for t in new_templates:
        final_paths.extend(expand(t, zip, **zip_dict))

    return final_paths


if __name__ == "__main__":

    targets = fetch_main_targets()

    print(targets)
