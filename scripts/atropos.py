from snakemake.shell import shell

shell.executable("bash")


adapters = ""
if snakemake.config["adapters"]:
    adapters = " -b ".join(snakemake.config["adapters"])

min_length = snakemake.config["min_read_length"]
quality = snakemake.config["fastq_quality"]
log = snakemake.log_fmt_shell()

if snakemake.config["paired_end"]:
    assert len(snakemake.input) == 2, "Paired end mode needs two input files, got: " + ", ".join(snakemake.input)

    command = "atropos --threads {snakemake.threads} -q {quality} -m {min_length} -b {adapters} -o {snakemake.output[0]} -p {snakemake.output[1]} {snakemake.input} {log}"

    shell(command)

else:

    assert len(snakemake.input) == 1, "Paired end mode needs one input file, got: " + ", ".join(snakemake.input)

    command = "atropos --threads {snakemake.threads} -q {quality} -m {min_length} -b {adapters} -o {snakemake.output[0]} {snakemake.input} {log}"

    shell(command)
