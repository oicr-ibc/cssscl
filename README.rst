CSSSCL: a taxonomic classifier for DNA sequences.
=================================================

**Project Leader** Vincent Ferretti

**Author** Ivan Borozan 


About:
======

**CSSSCL** is a python package that uses Combined Sequence Similarity Scores for accurate taxonomic CLassification of long and short reads.

Downloading and using cssscl is free, if you use cssscl or its code in your work please acknowledge cssscl by citing Borozan I, Ferretti V. "*CSSSCL: a python package that uses Combined Sequence Similarity Scores for accurate taxonomic CLassification of long and short sequence reads. Bioinformatics 2015.*"

This is important for us since obtaining grants is one significant way to fund planning and implementation for our project. Also if you find cssscl useful in your research feel free to let us know.  


Getting started: 
================


.. contents::
    :local:
    :depth: 1
    :backlinks: none


====================
Tested environments 
====================


| Distributor ID: Debian/Ubuntu
| Description: Debian GNU/Linux 8.1 (jessie) / Ubuntu 12.04.3 LTS 
| Release: 8.1 64-bit / 12.04 64-bit 
| Codename: jessie / precise


=================================
Dependencies on Debian and Ubuntu
=================================

Python: Python 2.7.3+ supported. No support for Python 3 at the moment.

In order to compile cssscl on Debian GNU/Linux 8.1 and Ubuntu 12.04 LTS the following packages need to be installed:

.. code-block:: bash

     $ sudo apt-get update
     $ sudo apt-get install build-essential python2.7 python2.7-dev g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev

Note: If you are testing cssscl using a VM please make sure that you have at least 1024 MB of RAM.


============
Installation
============

**Note:** if any of the following packages: **jellyfish**, **BLAST** or **plzip** are already INSTALLED on your system make sure that they are in your executable search path (i.e. PATH variable) (as shown in the examples below):

*BLAST*

.. code-block:: bash

     # e.g. PATH_TO_YOUR_BLAST=/home/user_x/blast/ncbi-blast-2.2.30+/bin
     $ export PATH=$PATH:PATH_TO_YOUR_BLAST 

*jellyfish*

.. code-block:: bash

     # e.g. PATH_TO_YOUR_jellyfish=/home/user_x/jellyfish-1.1.10/bin
     $ export PATH=$PATH:PATH_TO_YOUR_jellyfish 
 
*plzip*

.. code-block:: bash

     # e.g. PATH_TO_YOUR_plzip=/home/user_x/plzip-1.1/plzip
     $ export PATH=$PATH:PATH_TO_YOUR_plzip


To install the cssscl package you have two options:
-------------------------------------------------------

**Option A**

Install the cssscl package using the **Python's Virtual Environment** tool to keep the dependencies required by the cssscl package in a separate directory and to keep your global python dist- or site-packages directory clean and manageable.

1. Download the cssscl package

  .. code-block:: bash 

     $ tar -zxvf cssscl-1.0.tar.gz
     $ cd cssscl-1.0

2. Check that all packages necessary to run the cssscl are installed and are avaialble 

  .. code-block:: bash 

     $ ./cssscl_check_pre_installation.sh

**Note:** Run the 'cssscl_check_pre_installation.sh' script to check if all third party software is installed (namely pip, plzip, BLAST, jellyfish and mongoDB), the script will also install them if necessary. The script will also check if: python (and python-dev), libxml2-dev, libxslt-dev, gfortran, libopenblas-dev and liblapack-dev are installed. All the third party executables such as blastn, plzip and jellyfish will be installed in the cssscl-1.0/src/bin/ directory.  	     

3. Create a virtual environment for the cssscl program (e.g. name it 'csssclvenv')

  .. code-block:: bash 
 
     $ virtualenv csssclvenv

4. To begin using the virtual environment, it first needs to be activated:

  .. code-block:: bash 

     $ source csssclvenv/bin/activate

5. INSTALL cssscl as root 

  .. code-block:: bash 

     $ sudo pip install .
    
Note: this will install all the python modules necessary for running the cssscl package in the 'cssscl-1.0/csssclvenv/' directory. 


6. Configure mongodb

 .. code-block:: bash 

     $ cssscl configure 
    
Accept all the values prompted by default by pressing [ENTER]  

7. If you are done working in the virtual environment for the moment, you can deactivate it:

  .. code-block:: bash 

     $ deactivate


**Option B**
    
Install the cssscl package directly to your python global dist- or site-packages directory (CAUTION: some of the python packages on  your system might be updated if required by the cssscl package) 
            
1. Download the cssscl package 
   
   .. code-block:: bash 

     $ tar -zxvf cssscl-1.0.tar.gz
     $ cd cssscl-1.0

2. Check that all packages necessary to run the cssscl are installed and are avaialble 
	      
   .. code-block:: bash 

     $ ./cssscl_check_pre_installation.sh

3. INSTALL cssscl   

   .. code-block:: 
   
     $ sudo pip install .        


4. Configure mongodb

 .. code-block:: bash 

     $ cssscl configure 

Accept all the values prompted by default by pressing [ENTER]  


=================
Uninstall cssscl 
=================

**Note:** this will only work if you installed cssscl with the cmd 'sudo pip install .' as shown above.
          
 .. code-block:: bash 

     $ cd cssscl-1.0/
     $ ./cssscl_uninstall.sh 


=====
Usage
=====

**To test the classifier we have provided taxon and test data for you to download, as shown from the links provided below:**

Download taxon data:

 .. code-block:: bash 

     $ wget --no-check-certificate https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/taxon.tar.gz
     $ tar -zxvf taxon.tar.gz
    

Download test/train data:

 .. code-block:: bash 

     $ wget --no-check-certificate https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/test_data.tar.gz
     $ tar -zxvf test_data.tar.gz


**To run the cssscl classifier**

1. First build the necessary databases from the training set.

 .. code-block:: bash 
    
     $ cssscl build_dbs -btax -c -blast -nt 2 PATH_TO/test_data/TRAIN.fa PATH_TO/taxon/

By default all databases will be outputted to the DIR where the TRAIN.fa resides (note that all paths provided in the examples above are using absolute/full paths to the files/directories). The above command will build three databases (blast, compression and kmer dbs) for sequences in the training set.

The cssscl's ``build_dbs`` module requires two positional arguments to be provided: 


     | 1. a file in the fasta format (e.g. TRAIN.fa as in the example above) that specifies the collection of reference genomes composing the training set.
     |
     | 2. a directory (taxon/ in the example above) that specifies the location where the taxon data is stored (more specifically the directory should contain the following files: gi_taxid_nucl.dmp, names.dmp and nodes.dmp, these files can be downloaded from the NCBI taxonomy database at ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/).


The additional optional arguments used above have the following meaning:
------------------------------------------------------------------------

     | -btax, --build_taxonomy_data
     |            taxa


For more information please consult the cssscl's build_dbs help page by typing:

 .. code-block:: bash 

     $ cssscl build_dbs --help



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

