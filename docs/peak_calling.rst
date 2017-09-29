Peak calling
============

bincs can be used to find enriched regions in ChIP-Seq data using multiple
ChIP-Seq callers. It currently supports the peak callers macs2 and epic.

To call peaks use the target `peaks`:

.. code-block:: bash

   snakemake peaks

Choose peak callers by adding them to the list cs\_callers in your config.yaml:

.. code-block:: yaml

   peak_callers:
     - epic
     - macs2

The output of the peak target is stored in the folder
`{prefix}/data/peaks/{peak_caller}/{group}.csv`



.. bincs also supports CSAW, but since it finds regions that differ between
    conditions it is documented under differential enrichment and can not be used
    with the peaks target.
