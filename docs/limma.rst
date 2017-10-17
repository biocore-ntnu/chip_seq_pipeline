Differential analysis
=====================

bincs can perform differential analysis of your ChIP-Seq data. The regions
tested for differential enrichment are the de-novo regions found by the ChIP-Seq
callers specified in the config file. Therefore, the limma target will output
one list of differentially enriched regions per ChIP-Seq caller.

Example Output
--------------

Targets
-------

limma
"""""

Use limma to find de-novo regions that are differentially enriched according to
the sample sheet given.

Options
-------

The settings that affect the differential enrichment target are the following:

ChIP-Seq callers
~~~~~~~~~~~~~~~~

The different ChIP-Seq callers used will affect which de-novo regions are found.

Contrasts
~~~~~~~~~


The contrasts to test for differential enrichment between groups. The names used
must be the same as the column names in the design_matrix. By default the
contrasts are that each group in the design matrix is tested against the others
in pairs. So if you have the three groups WT, KO1 and KO2 in your sample sheet,
the contrasts WT-KO1, WT-KO2 and KO1-KO2 are tested by default.

.. code-block:: yaml

   contrasts:
     - 0h-3h
     - 0h-6h
     - 3h-6h


Design matrix
~~~~~~~~~~~~~
