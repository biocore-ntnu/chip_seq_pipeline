Quick-start
================================

bincs is written as a Snakemake_ pipeline.

.. _Snakemake: http://snakemake.readthedocs.io/en/stable/

See the `installation page`_ for instructions on how to install Snakemake.

.. _`installation page`: installation.html

The data for this quick-start example can be found at zenodo_.

.. _zenodo: https://zenodo.org/record/1008923#.Wd3cTxOCxlc

Go to the subfolder examples/quick_start under bincs and download the dataset.

.. code-block:: bash

   cd bincs/examples/quick_start
   wget https://zenodo.org/record/1008923/files/subsample_example_data.tar
   tar -xvf subsample_example_data.tar

bincs requires a config file and a sample sheet to run.
