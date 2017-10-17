Quick-start
===========

bincs is written as a Snakemake_ pipeline.

.. _Snakemake: http://snakemake.readthedocs.io/en/stable/

See the `installation page`_ for instructions on how to install Snakemake.

.. _`installation page`: installation.html

The data for this quick-start example can be found at zenodo_.

.. _zenodo: https://zenodo.org/record/1008923#.Wd3cTxOCxlc

Go to the subfolder examples/quick_start under bincs and download the dataset.

.. code-block:: bash

   cd bincs/examples/quick_start/files
   wget https://zenodo.org/record/1008923/files/subsample_example_data.tar
   tar -xvf subsample_example_data.tar

bincs requires a config file to know what settings to use, and a sample sheet to
know what files to use.

They are both in the :file:`examples/quick_start/` folder - they are called
:file:`config.yaml` and :file:`sample_sheet.txt`.

Now, let us invoke Snakemake with the config file given. Go to the main bincs
folder and then run Snakemake like so:

.. code-block:: bash

   cd ../../..
   snakemake -pr -j 48 --use-conda --configfile examples/quick_start/config.yaml peaks

This tells Snakemake to run the peaks target with 48 cores. If you have less
cores available, just change this parameter.

Running it might take a little while but in the end you should see the output:

.. code-block:: bash

  localrule peaks:
      input: projects/quick_start/data/peaks/epic/melanocyte.csv,
             projects/quick_start/data/peaks/macs2/melanocyte.csv,
             projects/quick_start/data/peaks/epic/keratinocyte.csv,
             projects/quick_start/data/peaks/macs2/keratinocyte.csv,
             projects/quick_start/data/peaks/epic/fibroblast.csv,
             projects/quick_start/data/peaks/macs2/fibroblast.csv

These are the files created by our invocation of Snakemake.

Let us try one more:

.. code-block:: bash

   snakemake -pr -j 48 --use-conda --configfile examples/quick_start/config.yaml log2_ratio_heatmaps
