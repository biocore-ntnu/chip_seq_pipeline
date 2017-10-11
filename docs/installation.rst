Installation
================================

To use bincs, you need snakemake and the conda package manager.

First we will install the conda pacakge manager. Go to
https://www.anaconda.com/download/ and select the Python 3+ version.

Then follow the installation instructions.

Next, we will need to install snakemake. We will do this using anaconda.

.. code-block:: bash

   # install snakemake from the channel bioconda
   conda install -c bioconda snakemake

Finally you can install bincs with the command

.. code-block:: bash

   git clone git@github.com:biocore-ntnu/chip_seq_pipeline.git bincs
