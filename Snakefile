import os

from leave_one_out.create_sample_sheets import create_sample_sheet
from utils.file_getters import read_sample_sheet, sample_vs_group
from utils.helpers import expand_zip
from snakemake.shell import shell

shell.executable("bash")

from itertools import product
import pandas as pd


from itertools import combinations


def make_contrasts(groups):

    cs = {"minus".join([a, b]): "-".join([a, b]) for a, b in combinations(groups, 2)}

    return cs



# if the config dict is empty, no config file was given on the command line
# perhaps this should give an error instead?
if not config:
    configfile: "config.yaml"


if config["tmux"]:
    from utils.helpers import error_if_not_using_tmux
    error_if_not_using_tmux()


prefix = config["prefix"]
sample_sheet = read_sample_sheet(config["sample_sheet"])
ss = sample_sheet


def find_filetype(ss):

    filetypes = ss.File.str.replace(".gz$", "").str.split(".", expand=True).iloc[:, -1]
    filetypes = list(filetypes.drop_duplicates())
    assert len(filetypes) == 1, "More than one filetype in sample sheet: " + ", ".join(filetypes)

    return filetypes[0]


def find_merge_lanes(ss, filetype):

    if filetype not in ["fastq", "fq"]:
        return False
    elif len(ss.Name) == len(ss.Name.drop_duplicates()):
        return False
    else:
        return True


# def find_paired_end(ss, filetype):

#     if filetype in ["bedpe"]:
#         return True
#     elif set(ss.Mate) == {1, 2}:
#         return True
#     elif filetype == "bam":
#         for f in ss.File:
#             result = subprocess.check_output("samtools view -c -f 1 {}".format(f))
#             print(result)

#     raise

filetype = find_filetype(ss)
merge_lanes = find_merge_lanes(ss, filetype)
# paired_end = find_paired_end(ss, filetype)

loo_ss = create_sample_sheet(ss)
loo_groups = list(loo_ss.Group.drop_duplicates())


if config.get("external_control_sample_sheet", "") and config["external_control_sample_sheet"]:
    ec_ss = pd.read_table(config["external_control_sample_sheet"], sep="\s+", header=0)
    ec_groups = list(ec_ss.Group.drop_duplicates())
    ec_samples = list(ec_ss.Name.drop_duplicates())
else:
    ec_ss, ec_samples, ec_groups = pd.DataFrame(), [], []


to_include = ["download/annotation", "download/chromsizes",
              "deeptools/bamcompare","deeptools/bigwig_compare",
              "deeptools/heatmap", "deeptools/profileplot",
              "deeptools/computematrix", "deeptools/bamcoverage",
              "deeptools/bigwig", "deeptools/multi_bigwig_summary",
              "deeptools/plot_coverage", "pca/pca",
              "deeptools/plot_fingerprint", "merge_lanes/merge_lanes",
              "compute_tss/compute_tss", "trim/atropos", "align/hisat2",
              "sort_index_bam/sort_index_bam", "bamtobed/bamtobed",
              "chip_seq/epic", "chip_seq/macs2", "chip_seq/csaw",
              "epic/epic_merge", "epic/epic_blacklist", "epic/epic_cluster",
              "epic/epic_count", "leave_one_out/compute_chip_over_input",
              "leave_one_out/violin_plots", "normalize/average_input",
              "normalize/divide_chip_input", "voom/voom", "limma/limma"] #, "voom/voom"]


path_prefix = config["prefix"]
chip_samples = list(ss[ss.ChIP == "ChIP"].Name.drop_duplicates())
input_samples = list(ss[ss.ChIP == "Input"].Name.drop_duplicates())
groups = list(sample_sheet.Group.drop_duplicates())
first_group = groups[0]
all_but_first_group = groups[1:]
more_than_one_group = None
if all_but_first_group:
    more_than_one_group = groups[1]

contrasts = make_contrasts(groups).values()

regions = ["CDS", "exon", "five_prime_UTR", "gene", "start_codon",
           "stop_codon", "stop_codon_redefined_as_selenocysteine", "three_prime_UTR",
           "transcript", "internal_exon"]

regular_regions = config["regions"] if config["regions"] else []

custom_regions = list(config["custom_regions"]) if config["custom_regions"] else []
all_regions = regular_regions + custom_regions


wildcard_constraints:
    sample = "({})".format("|".join(chip_samples + input_samples + ec_samples)),
    group = "({})".format("|".join(groups + loo_groups + ec_groups)),
    group1 = "({})".format("|".join(groups + loo_groups + ec_groups)),
    group2 = "({})".format("|".join(groups + loo_groups + ec_groups)),
    chip = "(chip|input|log2ratio|ChIP|Input)",
    region_type = "({})".format("|".join(regions + custom_regions)),
    caller = "({})".format("|".join(config["peak_callers"])),
    contrast = "({})".format("|".join(contrasts))


for rule in to_include:
    include: "rules/{rule}.rules".format(rule=rule)


# rule all:
#     input:
#         expand("{prefix}/data/peaks/csaw/{contrast}.raw", prefix=prefix, contrast=contrasts)


rule log2_ratio_heatmaps:
    input:
        expand("{prefix}/data/heatmap/{region_type}/{chip}/scale_regions/{chip}_{group}_{region_type}.png",
                group=groups, region_type=all_regions, chip="log2ratio", prefix=prefix),


rule input_heatmaps:
    input:
        expand("{prefix}/data/heatmap/{region_type}/{chip}/scale_regions/{chip}_{group}_{region_type}.png",
                group=groups, region_type=all_regions, chip="input", prefix=prefix)


rule chip_heatmaps:
    input:
        expand("{prefix}/data/heatmap/{region_type}/{chip}/scale_regions/{chip}_{group}_{region_type}.png",
               group=groups, region_type=all_regions, chip="input", prefix=prefix)


rule peaks:
    input:
        expand("{prefix}/data/peaks/{cs_caller}/{group}.csv", group=list(set(sample_sheet.Group)),
               cs_caller=config["peak_callers"], prefix=prefix)

rule input_profileplots:
    input:
        expand("{prefix}/data/profileplot/{region_type}_{chip}_scale_regions_{group}_profile_plot.png",
               group=groups, region_type=all_regions, chip="input", prefix=prefix)


rule log2_ratio_profileplots:
    input:
        expand("{prefix}/data/profileplot/{region_type}_{chip}_scale_regions_{group}_profile_plot.png",
                group=groups, region_type=all_regions, chip="log2ratio", prefix=prefix)


rule chip_profileplots:
    input:
        expand("{prefix}/data/profileplot/{region_type}_{chip}_{scaled}_{group}_reference_plot.png", scaled="tss tes".split(),
                group=groups, region_type=all_regions, chip="chip", prefix=prefix)


rule group_merged_chip_vs_merged_input_bigwig:
    input:
        expand("{prefix}/data/bigwigcompare/group_{group}_chip_vs_input.bigwig",
               group=groups, prefix=prefix)


rule chip_sample_vs_merged_input_bigwigs:
    input:
        expand("{prefix}/data/bigwigcompare/sample_{sample}_vs_merged_input.bw", sample=sample_sheet.loc[sample_sheet.ChIP == "ChIP"].Name, prefix=prefix)


rule chip_bigwigs:
    input:
        expand("{prefix}/data/bigwig/{sample}.bigwig", sample=sample_sheet.loc[sample_sheet.ChIP == "ChIP"].Name, prefix=prefix)


rule input_bigwigs:
    input:
        expand("{prefix}/data/bigwig/{sample}.bigwig", sample=sample_sheet.loc[sample_sheet.ChIP == "Input"].Name, prefix=prefix)


rule merged_input_bigwigs:
    input:
        expand("{prefix}/data/merged_bigwig/{group}_Input.bigwig", group=groups, prefix=prefix)


rule merged_chip_bigwigs:
    input:
        expand("{prefix}/data/merged_bigwig/{group}_ChIP.bigwig", group=groups, prefix=prefix)


rule limma_:
    input:
        expand("{prefix}/data/limma/{caller}_{contrast}_cutoff.toptable",
               caller=config["peak_callers"], prefix=prefix,
               contrast=contrasts)


rule pca_chip_vs_merged_input:
    """Create a PCA of all the ChIP samples log2-divided by the merged input from
       the samples' respective groups."""
    input:
        expand("{prefix}/data/plot_pca/pca_{multibigwig}.pdf",
               prefix=prefix, multibigwig="chip_vs_merged_input")


rule pca_individual:
    """Create a PCA of all the ChIP-samples and all Input-samples in the same plot.
    The data are RPKM-normalized."""
    input:
        expand("{prefix}/data/plot_pca/pca_{multibigwig}.pdf",
                prefix=prefix, multibigwig="individual")

rule pca_limma:
    """Create a PCA of the data that is input to limma for differential expression.
       This data has gone through several rounds of normalization."""
    input:
        expand("{prefix}/data/plot_pca/{caller}.pdf",
               prefix=prefix, caller=config["peak_callers"])


rule fingerprint_plot:
    input:
        expand("{prefix}/data/plot_fingerprint/{group}_fingerprint_deeptools.pdf",
               prefix=prefix, group=groups)

rule coverage_plot:
    input:
        expand("{prefix}/data/plot_coverage/{group}_coverage.pdf",
               prefix=prefix, group=groups)


loo_file = "{prefix}/data/loo/chip_over_input/{logfc}_{group}_{caller}_{contrast}_lo_info.ratios"

rule leave_one_out:
    input:
        expand(loo_file,
               prefix=prefix, group=loo_ss.Group.drop_duplicates(),
               caller=config["peak_callers"], contrast=contrasts,
               logfc="above below".split())


if not len(ss.Group.drop_duplicates()) == 1:

    svsg = sample_vs_group(ss, "ChIP")

    rule log2_ratio_sample_vs_group_bigwig:
        input:
            expand("{prefix}/data/bigwigcompare/sample_{sample}_vs_group_{group}.bigwig", zip, sample=svsg.Sample, group=svsg.OtherGroup, prefix=prefix)

    # cartesian product without duplicates
    group1, group2 = [], []
    for g1, g2 in product(groups, groups):
        if g1 != g2:
            group1.append(g1)
            group2.append(g2)


    rule log2_ratio_group_vs_group_bigwig:
        input:
            expand("{prefix}/data/bigwigcompare/group_{group1}_vs_group_{group2}.bigwig", zip,
                    group1=group1, group2=group2, prefix=prefix)

    rule log2_ratio_group_vs_group_heatmap:
        input:
               expand_zip("{prefix}/data/heatmap/{region_type}/{chip}/scale_regions/group_vs_group/{group1}_vs_{group2}_{region_type}.png",
                          zip_dict={"group1": group1, "group2": group2},
                          regular_dict={"region_type": all_regions, "chip": "log2ratio",
                                        "prefix": prefix})



# bad: should not hinge on having design_matrix available;
# we should make one for them
# if config["design_matrix"]:
#     design_matrix = pd.read_table(config["design_matrix"], sep="\s+")

#     rule differential_analysis:
#         input:
#             expand("{prefix}/data/limma/{group}_{coefficient}.toptable", group=groups,
#                 coefficient=design_matrix)

# else:
#     print("No design matrix given so target differential_analysis not available.")
