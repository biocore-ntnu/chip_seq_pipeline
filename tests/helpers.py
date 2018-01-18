from subprocess import run

import yaml

# from tempfile import TemporaryDirectory
# from os.path import abspath, dirname, join

import pandas as pd


def create_sample_sheet_files(sample_sheet):

    df = pd.read_table(sample_sheet, sep=" ")
    for f in df.File:
        cmd = "touch {f}".format(f=f)
        print(cmd)
        run(cmd)


def run_dag(targets, configfile, sample_sheet, snakefile="Snakefile", extras="", dryrun=True):

    # abs_configfile = abspath(configfile)
    # abs_sample_sheet = abspath(sample_sheet)

    # print(abs_configfile)
    # print(abs_sample_sheet)

    # with TemporaryDirectory() as tempdir:
    if dryrun:
        dryrun = "n"
    else:
        dryrun = ""

    cmd = ("snakemake -s {snakefile} --configfile {configfile}"
            " -{dryrun}p {targets} -F {extras} --config sample_sheet={sample_sheet}").format(**locals())
    print(cmd)
    return run(cmd, shell=True).returncode
