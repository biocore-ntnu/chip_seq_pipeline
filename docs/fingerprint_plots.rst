Fingerprint plots
=================

Fingerprint plots are used to determine whether the ChIP-signal can be separated
from the Input-signal.

It uses the deeptools tool plotFingerprint_.

.. _plotFingerprint: http://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html

Example Output
--------------

Targets
-------

fingerprint_plot
""""""""""""""""

This target is used to create a fingerprint for each group. The output file is
`{prefix}/data/plot_fingerprint/{group}_fingerprint_deeptools.pdf`
