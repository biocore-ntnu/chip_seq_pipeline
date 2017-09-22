Bigwigs
=======

bincs can create many types of bigwigs for you. A bigwig is a binary file that
can be used to visualize your data in genome browsers.

Targets
-------

There are several targets to produce bigwigs that display your data in different ways.

In the equations below, C is the total number of ChIP-files while I is the total
number of Input-files.

group_merged_chip_vs_merged_input_bigwig
""""""""""""""""""""""""""""""""""""""""

This target creates a bigwig for each group, where the pooled ChIP is log2-divided by the pooled Input.

.. math::

   log_{2} \frac{RPKM(\sum_{i=1}^{C} {ChIP_i}}{RPKM(\sum_{j=1}^{I} {Input_j})}

.. _chip_sample_vs_merged_input_bigwig:

chip_sample_vs_merged_input_bigwigs
"""""""""""""""""""""""""""""""""""

This target creates a bigwig of each ChIP-sample log2-divided against the
merged input.

.. math::

   log_{2} \frac{RPKM_{ChIP_i}}{RPKM(\sum_{j=1}^{I} {Input_j})}

The output of this target is `{prefix}/data/bigwigcompare/sample_{sample}_vs_merged_input.bw`

chip_bigwigs, input_bigwigs
"""""""""""""""""""""""""""

These targets create an RPKM-normalized bigwig of each ChIP (or Input) file
in your sample sheet.

The output of both these targets are stored in `{prefix}/data/bigwig/{sample}.bigwig`

merged_chip_bigwigs, merged_input_bigwigs
"""""""""""""""""""""""""""""""""""""""""

These targets create an RPKM-normalized bigwig of the merged ChIP (or Input)
files in your sample sheet.

.. math::

   RPKM (\sum_{i=1}^{C} {ChIP_i})

or

.. math::

   RPKM (\sum_{j=1}^{I} {Input_j})

The output of these targets is stored in
`{prefix}/data/merged_bigwig/{group}_ChIP.bigwig` and
`{prefix}/data/merged_bigwig/{group}_Input.bigwig`, respectively.

log2_ratio_group_vs_group_bigwig
""""""""""""""""""""""""""""""""

This target creates a log2-ratio bigwig for each combination input-normalized
group against each other. It requires more than one group in your sample sheet.

C1 and C2 is the total number of ChIP-files in group 1 and 2, respectively. I1
and I2 is the same for input.


.. math::

   log_{2} \frac{log_{2} \frac{RPKM(\sum_{i=1}^{C1} {ChIP Group1_i})}{RPKM(\sum_{j=1}^{I1} {Input Group1_j})}}
                {log_{2} \frac{RPKM(\sum_{k=1}^{C2} {ChIP Group2_k})}{RPKM(\sum_{l=1}^{I2} {Input Group2_l})}}


The output of this target is stored in `{prefix}/data/bigwigcompare/group_{group1}_vs_group_{group2}.bigwig`
