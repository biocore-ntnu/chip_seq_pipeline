# ChIP-Seq pipeline

[![Build Status](https://travis-ci.org/biocore-ntnu/chip_seq_pipeline.svg?branch=master)](https://travis-ci.org/biocore-ntnu/chip_seq_pipeline)

This is a collection of pipelines to analyse ChIP-Seq data.

Currently the pipeline has the following targets:

`peaks`, `bigwigs`. `profileplots`, `heatmaps` and `tss_tes_plots`.

* peaks uses the epic or macs2 peak callers to find peaks in your ChIP-Seq data.
* bigwigs creates unstranded bigwig files of your ChIP-Seq data
* profileplots creates a scaled plot of of genomic regions, such as TSSes or exons
* heatmaps creates heatmaps of your data
* tss\_tes\_plots creates unscaled plots of the regions before and after your tss_regions

## Usage

`snakemake --configfile <configfile_path> --use-conda -j <number_of_threads> <target_name>`

## Configuration

#### config.yaml

The pipeline has a config.yaml with many switches to customize the pipeline,
depending on whether the ChIP-Seq peaks are broad or narrow and whether the
reads are paired-end or not.

#### sample sheet

The config file takes a sample sheet which describes your experiment. It has the
columns File, Name, Group, ChIP, Mate.

* File is the path to the fastq file
* Name is the name of the sample the file belongs to
* Group is the name of the group the file belongs to. Files belonging to the
  same group are analyzed together.
* ChIP should be either ChIP or Input
* Mate is 1 or 2, which means that the mate is either mate1 or mate2. For single
  end, Mate should be 1.

## Install

```
git clone https://github.com/endrebak/chip_seq_pipeline
```
