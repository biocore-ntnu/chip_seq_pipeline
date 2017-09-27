from snakemake.shell import shell
import re

targets = re.findall("^rule (.*?):", open("Snakefile").read(), flags=re.MULTILINE)

for target in targets:
    cmd = "snakemake -n --configfile example/config.yaml --rulegraph {target} | dot -Tpng > docs/img/rulegraphs/{target}.png".format(**locals())
    print(cmd)
    shell(cmd)
