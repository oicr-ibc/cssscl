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

.. contents::
    :local:
    :depth: 1
    :backlinks: none


====================
Tested environments 
====================

Distributor ID: Debian / Ubuntu
Description:    Debian GNU/Linux 8.1 (jessie) / Ubuntu 12.04.3 LTS
Release:        8.1 64-bit / 12.04 64-bit
Codename:       jessie / precise


=================================
Dependencies on Debian and Ubuntu
=================================

- Python_ - Python 2.7.3+ supported. No support for Python 3 at the moment.

In order to compile cssscl on Debian GNU/Linux 8.1 and Ubuntu 12.04 LTS the following packages need to be installed:
--------------------------------------------------------------------------------------------------------------------
.. code-block:: bash

   $ sudo apt-get update
   $ sudo apt-get install build-essential python2.7 python2.7-dev g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev

Note: If you are testing cssscl using a VM please make sure that you have at least 1024 MB of RAM.


============
Installation
============

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

