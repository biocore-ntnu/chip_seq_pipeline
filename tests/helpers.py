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


def run_dag(targets, configfile, sample_sheet, snakefile="Snakefile", extras="", dryrun=True, ncores=1, forceall=True, configs=""):

    force = "-F" if forceall else ""

    if (dryrun == False) and ncores != 1:
        cores = "-j {ncores}".format(ncores=ncores)
    else:
        cores = ""

    # with TemporaryDirectory() as tempdir:
    if dryrun:
        dry = "n"
    else:
        dry = ""

    cmd = ("snakemake -s {snakefile} {cores} --configfile {configfile}"
            " -{dry}p {targets} {force} {extras} --config sample_sheet={sample_sheet} {configs}").format(**locals())
    print(cmd)
    return run(cmd, shell=True).returncode
