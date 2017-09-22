PCA
===

bincs can create three types of PCA for you.

Targets
-------

pca_individual
""""""""""""""

This target creates a PCA of all the ChIP-samples and all Input-samples in the same plot.
The data are RPKM-normalized.

The output path is `{prefix}/data/plot_pca/pca_{multibigwig}.pdf`

pca_chip_vs_merged_input
""""""""""""""""""""""""

This target creates a PCA of all the ChIP samples log2-divided by the merged
input. The normalization is the same as for the bigwig target
:ref:`chip_sample_vs_merged_input_bigwig`.

pca_limma
"""""""""

This target creates a PCA of the data that is input to limma for differential
expression. This data has gone through several rounds of normalization before
being input to limma.

The output target is `{prefix}/data/plot_pca/{caller}.pdf`
