Configuration files
===================

bincs requires two files to run, namely a configuration file and a
sample sheet.

A config file describes the pipeline settings, such as genome build, annotation
files and ChIP-Seq callers. The sample sheet describes the layout of your
experiment and the files which contain your biological data.

Sample sheet
------------

Let's look at the sample sheet first.

This file contains five columns: File, Name, Group, ChIP and Mate.


.. code-block:: csv

   File Name Group ChIP Mate
   download/ChIP_1_fibroblast.bed.gz fibroblast_1_ChIP fibroblast ChIP 1
   download/ChIP_2_fibroblast.bed.gz fibroblast_2_ChIP fibroblast ChIP 1
   download/ChIP_3_fibroblast.bed.gz fibroblast_3_ChIP fibroblast ChIP 1
   download/ChIP_1_keratinocyte.bed.gz keratinocyte_1_ChIP keratinocyte ChIP 1
   download/ChIP_2_keratinocyte.bed.gz keratinocyte_2_ChIP keratinocyte ChIP 1
   download/ChIP_3_keratinocyte.bed.gz keratinocyte_3_ChIP keratinocyte ChIP 1
   download/ChIP_1_melanocyte.bed.gz melanocyte_1_ChIP melanocyte ChIP 1
   download/ChIP_2_melanocyte.bed.gz melanocyte_2_ChIP melanocyte ChIP 1
   download/ChIP_3_melanocyte.bed.gz melanocyte_3_ChIP melanocyte ChIP 1
   download/Input_1_fibroblast.bed.gz fibroblast_1_Input fibroblast Input 1
   download/Input_2_fibroblast.bed.gz fibroblast_2_Input fibroblast Input 1
   download/Input_3_fibroblast.bed.gz fibroblast_3_Input fibroblast Input 1
   download/Input_1_keratinocyte.bed.gz keratinocyte_1_Input keratinocyte Input 1
   download/Input_2_keratinocyte.bed.gz keratinocyte_2_Input keratinocyte Input 1
   download/Input_3_keratinocyte.bed.gz keratinocyte_3_Input keratinocyte Input 1
   download/Input_1_melanocyte.bed.gz melanocyte_1_Input melanocyte Input 1
   download/Input_2_melanocyte.bed.gz melanocyte_2_Input melanocyte Input 1
   download/Input_3_melanocyte.bed.gz melanocyte_3_Input melanocyte Input 1


* **File**

   This is the path to the ChIP-Seq file. Valid filetypes are fastq(.gz), bed/bedpe(.gz)
   and bam.

* **Name**

   This is the name of the sample the file contains.

* **Group**

   This is the group the sample belongs to. The group variable is used to denote
   which files should be analyzed together and often corresponds to an
   experimental condition, timepoint or tissue type. For example, when calling
   peaks, you do not want to call peaks for all the data pooled, but rather find
   peaks for each group, fibroblast, keratinocyte and melanocyte so that you
   know which peaks are specific to which condition.

* **ChIP**

   This variable tells whether the file is a ChIP file or an Input file. Since
   ChIP-Seq data typically contain ~90-95% noise, having a background file for
   statistical comparisons is extremely important.

   Note that there is no requirement about the number of Input files being the
   same as the number of ChIP files. Also, ChIP and Input files are not matched,
   i.e., fibroblast_1_ChIP is not matched against fibroblast_1_Input above. In
   all analyses, the background signal is pooled, since it is supposed to be
   just background noise.

* **Mate**

   Mate is used to denote whether a file contains reads from the first mate or
   second in paired end data. It is only applicable when the file format used is
   fastq.


Configuration file
------------------

The configuration file contains a wealth of variables which indicate
user-configurable settings.

General 
 ~~~~~~~

Settings common to all the targets.

* **prefix**

.. code-block:: yaml

   # (Required)
   prefix: my_project

Path to store the intermediate and final results. Will be created if it does not exist.


* **tmux**

.. code-block:: yaml

   # (Not required)
   tmux: False

If set, ensures that Snakemake won't start unless you use the tmux session manager.

Sample Sheets 
 ~~~~~~~~~~~~~

The sample sheet describes the sequence files to use, their type (ChIP/Input) and which experimental condition they belong to.

* **sample_sheet**

.. code-block:: yaml

   # (Required)
   sample_sheet: sample_sheet.txt

Main sample sheet.


* **external_control_sample_sheet**

.. code-block:: yaml

   # (Not required)
   external_control_sample_sheet: 

A sample sheet of external controls. For example the same ChIP used in a different species.
Can be used to find non-specific binding by aligning these files to the genome used and seeing which bins get a lot of reads. These are candidates for blacklisting and tested according to a Poisson model.

Genomes And Annotations 
 ~~~~~~~~~~~~~~~~~~~~~~~

Which genomes and annotations to use

* **genome**

.. code-block:: yaml

   # (Required)
   genome: hg38

Which genome to use. Must be the UCSC name. Mostly used to automatically fetch genome and chromosome sizes.

Alignment 
 ~~~~~~~~~

These flags and options are only of interest if your filetype is fastq.

* **hisat2_index_prefix**

.. code-block:: yaml

   # (Not required)
   hisat2_index_prefix: 

The prefix of the hisat2 genome index. Premade indexes can be found at ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data These indexes are used to align the fastqs to the genome. If the input is not fastq, these files are not used.


* **hisat2_extra_flags**

.. code-block:: yaml

   # (Not required)
   hisat2_extra_flags: --no-spliced-alignment -k 1 --no-discordant --no-mixed

Additional flags to give to hisat2.


* **keep_multi_aligning_reads**

.. code-block:: yaml

   # (Not required)
   keep_multi_aligning_reads: True

If this flag is set, multi aligning reads from the alignment are removed (i.e. those with a NH:i field of 1.)


* **adapters**

.. code-block:: yaml

   # (Not required)
   adapters: 

A list of adapters to remove from your fastqs.


* **min_read_length**

.. code-block:: yaml

   # (Not required)
   min_read_length: 14

The minimum length of your reads after adapter removal.

Sequencing 
 ~~~~~~~~~~

heyo

* **paired_end**

.. code-block:: yaml

   # (Required)
   paired_end: False

Whether or not the reads are paired end.


* **fragment_length**

.. code-block:: yaml

   # (Required)
   fragment_length: 150

Estimated size of the sequenced fragments.


* **fastq_quality**

.. code-block:: yaml

   # (Not required)
   fastq_quality: 20

Minimum required fastq quality.

Heatmaps And Profileplots 
 ~~~~~~~~~~~~~~~~~~~~~~~~~

Settings for heatmaps and profileplots.

* **annotation_gff3**

.. code-block:: yaml

   # (Not required)
   annotation_gff3: 

The gzipped gencode annotation file to use. It is used to find the regions to be plotted in heatmaps and profileplots. If the path starts with www/http/ftp it will be downloaded, otherwise the path will be interpreted to be local.
This kind of high-quality annotation only exists for humans and mice. If your genome is any other, the appropriate regions will be downloaded from UCSC refGene for your genome.


* **tss_distance_gene**

.. code-block:: yaml

   # (Not required)
   tss_distance_gene: 3000

The distance to be shown upstream and downstream of any region with 'gene' in the name in the heatmaps and profileplots.


* **tss_distance_other**

.. code-block:: yaml

   # (Not required)
   tss_distance_other: 500

The distance to be shown upstream and downstream of any region without 'gene' in the name in the heatmaps and profileplots.


* **regions**

.. code-block:: yaml

   # (Not required)
   regions: ['genes', 'exons']

`regions` is a list of the regions to be plotted in the heatmaps and profileplots.
The following regions are valid: ```'CDS', 'exon', # 'five_prime_UTR', 'gene', 'start_codon', 'stop_codon', 'stop_codon_redefined_as_selenocysteine', 'three_prime_UTR', 'transcript', 'internal_exon'```
Using it requires that you either have the configuration variable `remote_annotation_gff3` or `local_annotation_gff3` defined.


* **custom_regions**

.. code-block:: yaml

   # (Not required)
   custom_regions: 

`custom_regions` is a map between region names and bed-files used to define the custom regions.
Custom region files can also be used to draw different lines and colors for different regions in the same plot. This can be done by adding a seventh column called deepTools_group to the file. If used, the file must be sorted on the column deepTools_group.

Differential Enrichment 
 ~~~~~~~~~~~~~~~~~~~~~~~

limma and stuff

* **contrasts**

.. code-block:: yaml

   # (Not required)
   contrasts: 

Contrasts to test for differential enrichment between groups. The names used must be the same as the column names in the design_matrix. By default each group in the design matrix is tested against the others in pairs. By default a design matrix with an intercept where each group in the sample sheet is a column in the design matrix. Command used is model.matrix(~0 + groups).


* **design_matrix**

.. code-block:: yaml

   # (Not required)
   design_matrix: 

The file used to describe the experimental setup.

Peak Calling 
 ~~~~~~~~~~~~

peak calling settings

* **peak_callers**

.. code-block:: yaml

   # (Not required)
   peak_callers: ['macs2', 'epic']

The peak callers used to find peaks. Valid values are epic and macs2.


* **epic_gap**

.. code-block:: yaml

   # (Not required)
   epic_gap: 3

The epic gaps-allowed parameter.


* **epic_window_size**

.. code-block:: yaml

   # (Not required)
   epic_window_size: 200

The epic window-size parameter.


* **macs2_broad**

.. code-block:: yaml

   # (Not required)
   macs2_broad: 

The macs2 --broad parameter.


* **macs2_extra_flags**

.. code-block:: yaml

   # (Not required)
   macs2_extra_flags: 

Any other flags you wish to pass to macs2


* **csaw_window_size**

.. code-block:: yaml

   # (Not required)
   csaw_window_size: 200

The window size csaw should use.


* **csaw_filter**

.. code-block:: yaml

   # (Not required)
   csaw_filter: 1

The minimum required count sum across libraries for a bin to be included in the analysis.

