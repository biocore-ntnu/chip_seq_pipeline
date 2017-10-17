PCA
===

bincs can do principal component analysis (PCA) of three types of data.

Example Output
--------------

.. figure:: img/pca/fibroblast.png

Targets
-------

Bincs can create two types of PCA for your raw data, and one for the normalized
data that is input to limma.

pca_chip_vs_merged_input
""""""""""""""""""""""""

Create a PCA of all the ChIP samples log2-divided by the merged input from the
samples' respective groups.

The output of this target is `{prefix}/data/plot_pca/pca_chip_vs_merged_input.pdf`.

pca_individual
""""""""""""""

Create a PCA of all the ChIP-samples and all Input-samples in the same plot. The
data are RPKM-normalized.

The output of this target is `{prefix}/data/plot_pca/pca_individual.pdf`.

pca_limma
"""""""""

Create a PCA of the data that is input to limma for differential expression.
This data has gone through several rounds of normalization.

The output of this target is `{prefix}/data/plot_pca/pca_{caller}.pdf`.
