Quick-start
================================

bincs is written as a Snakemake_ pipeline.

.. _Snakemake: http://snakemake.readthedocs.io/en/stable/

bincs requires a config file and a sample sheet to run.

To obtain the example files you should go to the folder example and run the command

.. code-block:: bash

   snakemake -s download_data_and_create_sample_sheet.Snakefile

This will download some data from gencode and create an appropriate sample sheet.
