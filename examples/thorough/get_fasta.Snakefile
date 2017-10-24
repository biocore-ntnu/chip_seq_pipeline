from glob import glob
from snakemake.shell import shell


shell.executable("bash")


fasta = "/mnt/cargo/genomes/hg19.fa"
prefix = "/mnt/scratch/endrebak/chip_seq_fasta"

samples = [f.replace(".bed.gz", "") for f in glob("*.bed.gz")]



rule all:
    input:
        expand("{prefix}/{sample}.fq.gz", sample=samples, prefix=prefix)


rule get_fasta:
    input:
        "{sample}.bed.gz"
    output:
        "{prefix}/{sample}.fa"
    shell:
        "bedtools getfasta -bed {input[0]} -fi {fasta} > {output[0]}"


rule remove_ns:
    input:
        "{prefix}/{sample}.fa"
    output:
        "{prefix}/{sample}_no_n.fa"
    run:
        import sys
        from Bio import SeqIO
        handle = open(input[0], "rU")
        filtered = (record for record in SeqIO.parse(handle, "fasta") if record.seq.count('N') == 0)
        output_handle = open(output[0], "w")
        SeqIO.write(filtered, output_handle, "fasta")
        output_handle.close()
        handle.close()


rule fasta_to_fastq_gz:
    input:
        "{prefix}/{sample}_no_n.fa"
    output:
        "{prefix}/{sample}.fq.gz"
    shell:
        "perl fasta_to_fastq.pl {input[0]} | gzip -9 > {output[0]}"
