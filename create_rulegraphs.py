from subprocess import run

from utils.helpers import fetch_main_targets

targets = fetch_main_targets()

# TODO: add outfolder option

for target in targets:
    cmd = ("snakemake {target} --rulegraph --configfile config.yaml"
     " -Fnp | dot -T svg > {target}_rulegraph.svg ").format(target=target)
    print(cmd)
    run(cmd, shell=True)
