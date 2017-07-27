import re

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
