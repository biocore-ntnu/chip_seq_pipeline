Peak calling
============

bincs can be used to find enriched regions in ChIP-Seq data using multiple
ChIP-Seq callers. It currently supports the peak callers macs2 and epic.

Example Output
--------------

The output will vary a little bit from peak caller to peak caller, but the three
first columns should always be chromosome, start and end. The output displays
the regions that were considered enriched for ChIP-Signal according to the peak
caller used.

.. code-block:: csv

   Chromosome Start End ChIP Input Score Log2FC P FDR
   chr1 87000 91199 648 70 162.33972644375527 2.5380417048293062 2.4622088068635633e-265 1.6146152112914843e-264
   chr1 234200 235999 261 72 60.63375010999239 1.18545571401722 4.51020459224775e-32 8.383370665861949e-32
   chr1 603600 604999 174 28 49.184656715786296 1.9630632926807723 2.3741640028312945e-49 5.6529520186485316e-49
   chr1 794200 843799 11250 2139 2609.3036128480867 1.7223913327163656 0.0 0.0
   chr1 845000 869399 10179 793 1632.4189178303059 3.0096058786178026 0.0 0.0
   chr1 870200 879599 3034 235 593.3292390015133 3.017963142571682 0.0 0.0
   chr1 893200 899599 616 158 167.01962876419088 1.2904805114074465 7.967564262501827e-84 2.5885901149853e-83
   chr1 901400 911999 2270 330 568.2288212275157 2.1096290868221095 0.0 0.0

The output-files are stored in `{prefix}/data/peaks/{cs_caller}/{group}.csv`.

Target
------

peaks
"""""

This target calls peaks with the options set in the config file.

Options
-------

You can configure which peak callers to use in the config file under the peak
calling header. There you will also see several configurable flags for each peak
caller. These are specific to each peak caller and you can read more about them
in the software's respective documentation.
