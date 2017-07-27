from subprocess import run

import pandas as pd


def create_sample_sheet_files(sample_sheet):

    df = pd.read_table(sample_sheet, sep=" ")
    for f in df.File:
        cmd = "touch {f}".format(f=f)
        print(cmd)
        run(cmd)

def run_dag(targets, config_file, snakefile="Snakefile", extras=""):

    cmd = ("snakemake -s {snakefile} --configfile {config_file}"
           " -np {targets} -F {extras}").format(**locals())
    print(cmd)
    return run(cmd, shell=True).returncode
