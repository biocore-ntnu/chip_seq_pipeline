.. epic documentation master file, created by
   sphinx-quickstart on Wed Jul 19 11:24:44 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

manmoth: the biocore-NTNU ChIP-Seq pipeline
===========================================

manmoth is a highly-configurable pipeline for ChIP-Seq analysis that allows you
to perform a wealth of different quality controls, create many different kinds
of graphs and perform peak calling with multiple peak callers.

manmoth does not require any programming skill to use, all you need to do is
fill out a sample sheet and config file. It does not have any external
dependencies except conda and snakemake; all the required software is downloaded
automatically.

All the analyses

Novel features
--------------

* manmoth allows you to do leave one out analyses to find out which peak caller
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

..
   .. toctree::
      :caption: epic
      :hidden:
      :maxdepth: 2

      basic_intro
      options
      output_files

   .. toctree::
      :caption: epic-tools
      :hidden:
      :maxdepth: 2

      epic_merge
      epic_cluster
      epic_count
      epic_blacklist


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
