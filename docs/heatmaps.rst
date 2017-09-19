Heatmaps
========

bincs can be used to create heatmaps of predefined regions in your data.

There are three targets for creating heatmaps: log2\_ratio\_heatmaps,
chip_heatmaps and input_heatmaps.

The region\_types list in the config file can be used to select which regions
should be graphed. Using this list relies on having a gff3 annotation file set
in your config, either under local_annotation_gff3 or remote_annotation_gff3.

.. code-block:: yaml

   # custom region types can be any of the following: "CDS", "exon",
   # "five_prime_UTR", "gene", "start_codon", "stop_codon",
   # "stop_codon_redefined_as_selenocysteine", "three_prime_UTR", "transcript",
   # "internal_exon"
   region_types:
     - exon
     - gene
     - internal_exon

If a you do not have a gff3 annotation available for your genome, you can use
the region\_files option in the config. This is a list of the region type name
and then the path to a bed file denoting the regions.

.. code-block:: yaml

   region_files:
     quantiles_internal_exon: test_data/WT.internal_exon.bed
     quantiles_exon: test_data/WT.exon.bed
     quantiles_gene: test_data/WT.gene.bed

The bed files should have a header and at least six columns. You can give an
optional seventh column which must be called deepTools_group. It is used by
deeptools to show each group separately in the profileplots.

.. code-block:: bash

   #chrom  start   end     name    score   strand  deepTools_group
   chr2    241252956       241253035       exon:ENST00000391975.5:11       .       -       0
   chr9    133103746       133103833       exon:ENST00000424572.1:5        .       -       0
   chr17   32208106        32208309        exon:ENST00000584692.1:3        .       +       0
   chr17   32207511        32207563        exon:ENST00000584692.1:2        .       +       0
   chr17   82236728        82236872        exon:ENST00000584689.5:3        .       +       0
   chr17   82235982        82236231        exon:ENST00000584689.5:2        .       +       0
   chr9    133104262       133104331       exon:ENST00000424572.1:4        .       -       0
   chr9    133105931       133106016       exon:ENST00000424572.1:3        .       -       0
   chr9    133106644       133106748       exon:ENST00000424572.1:2        .       -       0
   ...
   chr5    122391088       122391191       exon:ENST00000509154.6:3        .       +       75-100
   chr5    122336787       122336904       exon:ENST00000509154.6:2        .       +       75-100
   chr1    160282038       160282200       exon:ENST00000392220.2:5        .       -       75-100
   chr1    160282416       160282502       exon:ENST00000392220.2:4        .       -       75-100
   chr1    160282943       160283109       exon:ENST00000392220.2:3        .       -       75-100
   chr1    160283529       160283639       exon:ENST00000392220.2:2        .       -       75-100
   chr12   98832028        98832136        exon:ENST00000552748.5:2        .       -       75-100
   chr12   98829173        98829353        exon:ENST00000552748.5:3        .       -       75-100
   chr4    59429   59556   exon:ENST00000509152.3:2        .       +       75-100
   chr1    36307769        36307825        exon:ENST00000505871.6:3        .       +       75-100
