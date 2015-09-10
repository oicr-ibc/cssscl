CSSSCL: a taxonomic classifier for DNA sequences.
=================================================

**Project Leader** Vincent Ferretti

**Author** Ivan Borozan 


About:
======

**CSSSCL** is a python package that uses Combined Sequence Similarity Scores for accurate taxonomic CLassification of long and short reads.

.. Downloading and using ``cssscl`` is free. If you use ``cssscl`` or its code in your work, please acknowledge ``cssscl`` by citing Borozan I and Ferretti V. "*CSSSCL: a python package that uses Combined Sequence Similarity Scores for accurate taxonomic CLassification of long and short sequence reads. Bioinformatics 2015.*"

.. This is important for us since obtaining grants is one significant way to fund planning and implementation for our project. Also if you find ``cssscl`` useful in your research feel free to let us know.  


Tested environments 
====================


   | Distributor ID: Debian/Ubuntu
   | Description: Debian GNU/Linux 8.1 (jessie) / Ubuntu 12.04.3 LTS 
   | Release: 8.1 64-bit / 12.04 64-bit 
   | Codename: jessie / precise


**We have setup tree ways for installing and testing the** ``cssscl`` **package:**


1. `Quick deployment using Docker (small file) <https://github.com/oicr-ibc/cssscl/wiki/Quick-deployment-and-testing-using-Docker>`_.

2. `Quick deployment using a VM (large file ~ 5GB) <https://github.com/oicr-ibc/cssscl/wiki/Quick-deployment-and-testing-using-a-VM>`_.

3. System wide installation using the source code (see installation instructions below).



Getting started
===============

.. contents::
    :local:
    :depth: 2
    :backlinks: none


============
Installation
============


Dependencies on Debian and Ubuntu
---------------------------------

**If you are using the Python's Virtual Environment**

In order to compile ``cssscl`` on Debian GNU/Linux 8.1 and Ubuntu 12.04 LTS the following packages need to be installed:

.. code-block:: bash

     $ sudo apt-get update
     $ sudo apt-get install build-essential g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev

**see** `Option A <https://github.com/oicr-ibc/cssscl#install-cssscl-option-a>`_ **below**

**If you are not using the Python's Virtual Environment**

Python: Python 2.7.3+ supported. No support for Python 3 at the moment.

In order to compile ``cssscl`` on Debian GNU/Linux 8.1 and Ubuntu 12.04 LTS the following packages need to be installed:

.. code-block:: bash

     $ sudo apt-get update
     $ sudo apt-get install build-essential python2.7 python2.7-dev g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev

**see** `Option B <https://github.com/oicr-ibc/cssscl#install-cssscl-option-b>`_ **below**


Install ``cssscl`` Option A
---------------------------

Install the ``cssscl`` package using the **Python's Virtual Environment** tool to keep the dependencies required by the ``cssscl`` package in a separate directory and to keep your global python dist- or site-packages directory clean and manageable as shown below:

**Note:** if any of the following packages: **jellyfish**, **BLAST** or **plzip** are **already installed** on your system make sure that they are in your executable search path (i.e. PATH variable) (as shown in the examples below):


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


1. Download the ``cssscl`` package

  .. code-block:: bash 
   
     # use wget 
     $ wget --no-check-certificate https://github.com/oicr-ibc/cssscl/archive/master.tar.gz
     $ tar -zxvf master.tar.gz; mv cssscl-master cssscl 
     # or use git clone  
     $ git clone git@github.com:oicr-ibc/cssscl.git


2. Check that all packages necessary to run the ``cssscl`` are installed and are available by running the ``cssscl_check_pre_installation.sh`` script (only for Ubuntu/Debian distributions). 

  .. code-block:: bash 
    
     $ cd cssscl
     $ ./cssscl_check_pre_installation.sh

**Note:** for more information regarding the ``cssscl_check_pre_installation.sh`` script see `here <https://github.com/oicr-ibc/cssscl/wiki/cssscl_check_pre_installation>`_.

3. In the ``cssscl``  ``directory`` create a virtual environment (e.g. name it ``csssclvenv``)

  .. code-block:: bash 
 
     $ virtualenv csssclvenv


4. To begin using the virtual environment, it first needs to be **activated**:

  .. code-block:: bash 

     $ source csssclvenv/bin/activate


5. Install ``cssscl`` as root 

  .. code-block:: bash 

     $ sudo pip install .
    
**Note:** this will install all the python modules necessary for running the ``cssscl`` package in the ``cssscl/csssclvenv/`` directory. 


6. Configure ``cssscl``

 .. code-block:: bash 

     $ cssscl configure 
    

Accept all the values prompted by default by pressing [ENTER]  


**Note:** If you are done working in the virtual environment, you can deactivate it as shown below. 

  .. code-block:: bash 

     $ deactivate

If you would like to run the ``cssscl`` program again (and you have deactivated the python virtual environment) you will need to **activate** it again as shown above. 


Install ``cssscl`` Option B
---------------------------

Install the ``cssscl`` package directly to your python global dist- or site-packages directory as shown below (**CAUTION: some of the python packages on your system might be updated if required by the** ``cssscl`` **package**):

**Note:** if any of the following packages: **jellyfish**, **BLAST** or **plzip** are **already installed** on your system make sure that they are in your executable search path (i.e. PATH variable) (as shown in the examples below):


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

       
1. Download the ``cssscl`` package 
   
   .. code-block:: bash 

     # use wget 
     $ wget --no-check-certificate https://github.com/oicr-ibc/cssscl/archive/master.tar.gz
     $ tar -zxvf master.tar.gz; mv cssscl-master cssscl 
     # or use git clone  
     $ git clone git@github.com:oicr-ibc/cssscl.git

2. Check that all packages necessary to run the ``cssscl`` are installed and are avaialble by running the ``cssscl_check_pre_installation.sh`` script (only for Ubuntu/Debian distributions). 
	      
   .. code-block:: bash 

     $ cd cssscl
     $ ./cssscl_check_pre_installation.sh

**Note:** for more information regarding the ``cssscl_check_pre_installation.sh`` script please see `here <https://github.com/oicr-ibc/cssscl/wiki/cssscl_check_pre_installation>`_.


3. Install ``cssscl`` as root  

   .. code-block:: 
   
     $ sudo pip install .        


4. Configure ``cssscl`` 

   .. code-block:: bash 

     $ cssscl configure 

Accept all the values prompted by default by pressing [ENTER]  


Additional instructions for non-automated installation of third party software necessary for running the ``cssscl`` package
---------------------------------------------------------------------------------------------------------------------------

In case the **cssscl_check_pre_installation.sh** script (see the installation subsections above) fails please read the info below for the manual installation of individual third party software:

Necessary Python modules: 

- BioPython_ - Tools for biological computation.
- PyMongo_ - Python module needed for working with MongoDB (PyMongo = 2.8)
- Sklearn_ - Machine Learning in Python
- Numpy_ - NumPy is the fundamental package for scientific computing with Python
- Cython_ - Cython is an optimising static compiler for both the Python programming language and the extended Cython programming language (based on Pyrex)
- SciPy_ - SciPy is a Python-based ecosystem of open-source software for mathematics, science, and engineering. In particular, these are some of the core packages:

.. _Python: http://www.python.org
.. _BioPython: http://biopython.org/wiki/Main_Page
.. _PyMongo: http://api.mongodb.org/python/2.8/
.. _Sklearn: http://scikit-learn.org/stable/
.. _Numpy: http://www.numpy.org/
.. _Cython: http://cython.org/
.. _SciPy: http://www.scipy.org/


**Installing python modules using pip manually:**

 .. code-block:: bash 

     $ pip install cython
     $ pip install numpy
     $ pip install pymongo==2.8
     $ pip install biopython
     $ pip install scikit-learn
     $ pip install scipy    

**Third party software:**

**BLAST (version 2.2.30+ and higher)**
Basic Local Alignment Search Tool.
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

**JELLYFISH (version 1.1.+ but not 2.0.+)**
JELLYFISH is a tool for fast, memory-efficient counting of k-mers in DNA.
http://www.cbcb.umd.edu/software/jellyfish/

**PLZIP (version 1.1+)**
Plzip is a massively parallel (multi-threaded) lossless data compressor based on the lzlib compression library, with a user interface similar to the one of lzip, bzip2 or gzip. 
http://download.savannah.gnu.org/releases/lzip/plzip/

**Note:** that the classification results in the paper were obtained using: Plzip 1.1 using Lzlib 1.5

**To compile Plzip 1.1 and Lzlib 1.5:**

1. Donwload lzlib-1.5.tar.gz 

.. code-block:: bash 

     $ wget --no-check-certificate http://download.savannah.gnu.org/releases/lzip/lzlib/lzlib-1.5.tar.gz 

2. Install lzlib-1.5:

.. code-block:: bash 

     $ gunzip lzlib-1.5.tar.gz
     $ tar -xvf lzlib-1.5.tar
     $ cd lzlib-1.5
     $ ./configure
     $ make
     $ make install


3. Donwload Plzip 1.1 

.. code-block:: bash 

     $ wget --no-check-certificate  http://download.savannah.gnu.org/releases/lzip/plzip/plzip-1.1.tar.gz

4. Install Plzip

.. code-block:: bash 

     $ gunzip plzip-1.1.tar.gz
     $ tar -xvf plzip-1.1.tar 
     $ cd plzip-1.1 
     $ ./configure
     $ make
     $ make install

For more information about plzip consult:
http://www.nongnu.org/lzip/manual/plzip_manual.html

and for memory required to compress and decompress: 
http://www.nongnu.org/lzip/manual/plzip_manual.html#Memory-requirements


**Make sure that JELLYFISH, BLAST and Plzip are in your executable search path (see the examples below):**

.. code-block:: bash 

     # for example 
     $ export PATH=$PATH:PATH_TO_BLAST/blast/ncbi-blast-2.2.30+/bin
     $ export PATH=$PATH:PATH_TO_jellyfish/jellyfish-1.1.10/bin
     $ export PATH=$PATH:PATH_TO_plzip/plzip-1.1/plzip
   

**Install MongoDB**

*Ubuntu*

You will first need to install Mongodb (ignore mongodb installation if mongodb is already installed jump to 2. Set up cssscl):

MongoDB should be installed using the following set of instructions (see also mongodb installation):

First add the 10gen GPG key, the public gpg key used for signing these packages. It should be possible to import the key into apt's public keyring with a command like this:

.. code-block:: bash 

     $ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv 7F0CEB10

Add this line verbatim to your /etc/apt/sources.list:

.. code-block:: bash 

     $ deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen

In order to complete the installation of the packages, you need to update the sources and then install the desired package

.. code-block:: bash 

     $ sudo apt-get update 
     $ sudo apt-get install mongodb-10gen=2.4.14


*Debian*

.. code-block:: bash 

     $ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv 7F0CEB10
     $ echo 'deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen' | tee -a /etc/apt/sources.list
     $ apt-get update 
     $ apt-get install mongodb-10gen=2.4.14


Uninstall ``cssscl`` 
---------------------

**Note:** this will only work if you installed cssscl with the cmd ``sudo pip install .`` as shown in the Installation section above. 
          
 .. code-block:: bash 

     $ cd cssscl/
     $ ./cssscl_uninstall.sh 


==========
User Guide
==========

Download taxon and test data
----------------------------

Download taxon data:

 .. code-block:: bash 

     $ wget --no-check-certificate https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/taxon.tar.gz
     $ tar -zxvf taxon.tar.gz
    

Download test/train data:

 .. code-block:: bash 

     $ wget --no-check-certificate https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/test_data.tar.gz
     $ tar -zxvf test_data.tar.gz


Example 1 - run the ``cssscl`` classifier without the optimization using the taxon data and the test set provided
-----------------------------------------------------------------------------------------------------------------

1. Build the necessary databases from the training set

 .. code-block:: bash 
     
     $ cssscl build_dbs -btax -c -blast -nt 2 PATH_TO/test_data/TRAIN.fa PATH_TO/taxon/

(the whole process should take ~ 37 min using 2 CPUs)

By default all databases will be outputted to the directory where the TRAIN.fa resides (note that all paths provided in the examples above are using absolute/full paths to the files/directories). The above command will build three databases (blast, compression and the kmer database) for sequences in the training set.

The ``cssscl's`` ``build_dbs`` module requires two positional arguments to be provided: 

      | 1. a **file** in the fasta format (e.g. TRAIN.fa as in the example above) that specifies the collection of reference genomes composing the training set.
      |
      | 2. a **directory** (taxon/ in the example above) that specifies the location where the taxon data is stored (more specifically the directory should contain the following files: gi_taxid_nucl.dmp, names.dmp and nodes.dmp, these files can be downloaded from the NCBI taxonomy database at ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/).

The information about the additional optional arguments used in the command line above is provided `here <https://github.com/oicr-ibc/cssscl/wiki/build_dbs>`_.

For more information please consult the ``cssscl's``  ``build_dbs`` help page by typing:

 .. code-block:: bash 

      $ cssscl build_dbs --help


2. Perform the classification the test test set

 .. code-block:: bash 

      # use cssscl to classify sequences in TEST.fa 
      $ cssscl classify -c -blast blastn -tax genus -nt 2 PATH_TO/test_data/test/TEST.fa PATH_TO/test_data/
 
(the whole process should take ~ 29 min using 2 CPUs)

Note that in the above example the output file ``cssscl_results_genus.txt`` with classification results will be located in the directory where the TEST.fa resides. 

**Note**: For the `test set data <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/test_data.tar.gz>`_ provided above the values of the parameters used in the model have already been optimized and are included as part of the test set (see the ``optimum_kmer`` directory in the ``test_set/`` directory provided). Thus for the test dataset the optimization is not required to be performed prior to running the classifier. On how to run the classifier by performing the optimization stage first please see the step 3 below. 

The ``cssscl's``  ``classify`` module requires two positional arguments to be provided: 

      | 1. a **file** with test data with sequences in the FASTA format for classification (e.g. TEST.fa as in the example above)
      |
      | 2. a **directory** where the databases (built using the training set) reside


**Note**: This will run the classifier with all the similarity measures (including the compression and the blast measure) as described in:  Borozan I et al. "*Integrating alignment-based and alignment-free sequence similarity measures for biological sequence classification.*"  Bioinformatics. 2015 Jan 7. pii: btv006.


The information about the additional optional arguments used in the command line above is provided `here <https://github.com/oicr-ibc/cssscl/wiki/classify>`_.


For more information please consult the ``cssscl's``  ``classify`` help page by typing 

 .. code-block:: bash 

      $ cssscl classify --help 


Example 2 - perform the classification by optimizing the ``cssscl's`` parameter values first
--------------------------------------------------------------------------------------------

1. Build the necessary databases from the training set

**Note**: Only do this is you did not already built the database in Example 1 above.

 .. code-block:: bash 
     
     $ cssscl build_dbs -btax -c -blast -nt 2 PATH_TO/test_data/TRAIN.fa PATH_TO/taxon/

(the whole process should take ~ 37 min using 2 CPUs)


2. Perform the classification by optimizing the ``cssscl's`` parameter values first

 .. code-block:: bash 

      $ cssscl classify -c -blast blastn -opt -tax genus -nt 8 PATH_TO/test_data/test/TEST.fa PATH_TO/test_data/


More information about the optimization can be found `here <https://github.com/oicr-ibc/cssscl/wiki/optimization>`_. 

Note that the optimization phase will take considerably longer when ``-c`` (compression) argument is used as mentioned in the section **Note regarding the compression measure** below.

The information about the additional optional arguments used in the command line above is provided `here <https://github.com/oicr-ibc/cssscl/wiki/classify_opt>`_.


Note regarding the compression measure
--------------------------------------

The use of the compression measure will slow down considerably the optimization and the classification parts because of the running 
time complexity ~ O(n*n) (for the optimization phase) and  ~ O(n*m) for the classification phase, where n and m are respectively 
the number of sequences in the training and test sets. Thus the compression measure should only be used with smaller genome 
databases (e.g. viruses) and/or with smaller datasets (i.e. smaller number of reads/contigs to classify).



==================
Supplementary data
==================

1. **Accompanying supplementary file** to the Bioinformatics 2015 paper "*CSSSCL: a python package that uses Combined Sequence Similarity Scores for accurate taxonomic CLassification of long and short sequence reads. Bioinformatics 2015*" can be found `here <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/supplementary_data.pdf>`_.

2. **Test data:**

      Genome sequences: `test data <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/test_data.tar.gz>`_

      Taxon Data: `taxon <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/taxon.tar.gz>`_


3. **Links to the three full datasets used to generate the results presented in Table 1 on pg.2 of the manuscript are shown below:**

      `Viral <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/viral/train_test_viral_full_data.tar.gz>`_ - Viral sequences (full dataset) used in the paper.

      `Bacterial I <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/bacterial1/bacterial1.tar.gz>`_ - dataset I Bacterial sequences (full dataset) used in the paper.

      `Bacterial II <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/data/bacterial2/bacterial2.tar.gz>`_ - dataset II Bacterial sequences (full dataset) used in the paper. 



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

