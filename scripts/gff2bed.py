from snakemake.shell import shell

if snakemake.wildcards.region_type != "internal_exon":
    cmd = r"""gff2bed < {snakemake.input[0]} | grep -E '\b{snakemake.wildcards.region_type}\b' |
            perl -a -ne'/gene_id=(.*?);/ && print "$F[0] $F[1] $F[2] $F[3] $F[4] $F[5] $1\n"' | tr " " "\t" > {snakemake.output[0]}"""
    print(cmd.format(**vars()))
    shell(cmd)

else:
    cmd = r"""gff2bed < {snakemake.input[0]} | grep -E '\bexon\b' |
            perl -a -ne'/gene_id=(.*?);/ && print "$F[0] $F[1] $F[2] $F[3] $F[4] $F[5] $1\n"' | tr " " "\t" > {snakemake.output[0]}.tmp"""
    print(cmd.format(**vars()))
    shell(cmd)
    shell("python scripts/compute_internal_exons.py {snakemake.output[0]}.tmp > {snakemake.output[0]}")
    shell("rm {snakemake.output[0]}.tmp")
