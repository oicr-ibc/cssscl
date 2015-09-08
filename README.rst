CSSSCL: a taxonomic classifier for DNA sequences
======

**Project Leader** Vincent Ferretti

**Author** Ivan Borozan 

About:
======

**CSSSCL** is a python package that uses Combined Sequence Similarity Scores for accurate taxonomic CLassification of long and short reads.

Downloading and using cssscl is free, if you use cssscl or its code in your work 
please acknowledge cssscl by citing Borozan I, Ferretti V. "CSSSCL: a python package that uses Combined Sequence Similarity Scores for accurate taxonomic CLassification of long and short sequence reads.", Bioinformatics 

This is important for us since obtaining grants is one significant way to fund planning 
and implementation for our project. Also if you find cssscl useful in your research feel 
free to let us know.  


Getting started: 
================


.. contents::
    :local:
    :depth: 1
    :backlinks: none


====================
Tested environments 
====================

Distributor ID: Debian / Ubuntu__
Description: Debian GNU/Linux 8.1 (jessie) / Ubuntu 12.04.3 LTS__ 
Release: 8.1 64-bit / 12.04 64-bit__
Codename: jessie / precise__


=================================
Dependencies on Debian and Ubuntu
=================================

- Python: Python 2.7.3+ supported. No support for Python 3 at the moment.

In order to compile cssscl on Debian GNU/Linux 8.1 and Ubuntu 12.04 LTS the following packages need to be installed:

.. code-block:: bash

   $ sudo apt-get update
   $ sudo apt-get install build-essential python2.7 python2.7-dev g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev

Note: If you are testing cssscl using a VM please make sure that you have at least 1024 MB of RAM.


============
Installation
============

If any of the following packages: **jellyfish**, **BLAST** or **plzip** are already INSTALLED on your system make sure that they are in your 
executable search path (i.e. PATH variable) (as shown in the examples below):

- BLAST
.. code-block:: bash
  
   # e.g. PATH_TO_YOUR_BLAST=/home/user_x/blast/ncbi-blast-2.2.30+/bin
   $ export PATH=$PATH:PATH_TO_YOUR_BLAST 

- jellyfish
.. code-block:: bash

   # e.g. PATH_TO_YOUR_jellyfish=/home/user_x/jellyfish-1.1.10/bin
   $ export PATH=$PATH:PATH_TO_YOUR_jellyfish 
 
- plzip
.. code-block:: bash

   # e.g. PATH_TO_YOUR_plzip=/home/user_x/plzip-1.1/plzip
   $ export PATH=$PATH:PATH_TO_YOUR_plzip

To install the cssscl package you now have two options:

1. Install the cssscl package using the Python's Virtual Environment tool to keep the dependencies required by the cssscl package in a separate directory and to keep your global python dist- or site-packages directory clean and manageable.
 * Download the cssscl package
 .. code-block:: bash 
     $ wget --no-check-certificate https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/code/cssscl-1.0.tar.gz
     $ tar -zxvf cssscl-1.0.tar.gz
     $ cd cssscl-1.0
 * CHECK THAT ALL PACKAGES NECESSARY TO RUN THE cssscl ARE INSTALLED AND ARE AVAILABLE
 .. code-block:: bash 
     $ ./cssscl_check_pre_installation.sh
     Note: Run the 'cssscl_check_pre_installation.sh' script to check if all third party software is installed (namely pip, plzip,
     BLAST, jellyfish and mongoDB), the script will also install them if necessary. The script will also check if: python (and  
     python-dev), libxml2-dev, libxslt-dev, gfortran, libopenblas-dev and liblapack-dev are installed. 
     All the third party executables such as blastn, plzip and jellyfish will be installed in the cssscl-1.0/src/bin/ directory.  	     


=====
Usage
=====

=====================
License and Copyright
=====================
Licensed under the GNU General Public License, Version 3.0. See LICENSE for more details.

Copyright 2015 The Ontario Institute for Cancer Research.

===============
Acknowledgement
===============

This project is supported by the Ontario Institute for Cancer Research
(OICR) through funding provided by the government of Ontario, Canada.

