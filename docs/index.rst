.. bincs

bincs: the biocore-NTNU ChIP-Seq pipeline
===========================================

bincs is a highly-configurable pipeline for ChIP-Seq analysis that allows you
to perform a wealth of different quality controls, create many different kinds
of graphs and perform peak calling with multiple peak callers.

bincs does not require any programming skill to use, all you need to do is
fill out a sample sheet and config file. It does not have any external
dependencies except conda and snakemake; all the required software is downloaded
automatically.

Novel features
--------------

* bincs allows you to do leave one out analyses to find out which peak caller
  finds the most reproducible peaks in your data
* it allows you to easily perform statistical analyses of different experimental
  conditions



-----


.. toctree::
   :caption: Getting Started
   :hidden:
   :maxdepth: 2

   installation
   quick_start
   basic_intro

.. toctree::
   :caption: Configuration Files
   :hidden:
   :maxdepth: 2

   configuration_files

.. toctree::
   :caption: Quality Control
   :hidden:
   :maxdepth: 2



.. toctree::
   :caption: Peak calling
   :hidden:
   :maxdepth: 2

   peak_calling


.. toctree::
   :caption: Graphs
   :hidden:
   :maxdepth: 2

   heatmaps




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
