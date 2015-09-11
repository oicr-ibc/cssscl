============
Installation Guide
============


Option A - Installing  ``cssscl`` using the Python's Virtual Environment
---------------------------

We recommend to install the ``cssscl`` package using the **Python's Virtual Environment** tool to keep the dependencies required by the ``cssscl`` package in a separate directory and to keep your global python dist- or site-packages directory clean and manageable as shown below:

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

*Step 1*. Install dependencies on Debian and Ubuntu

In order to compile ``cssscl`` on Debian GNU/Linux 8.1 and Ubuntu 12.04 LTS the following packages need to be installed:

.. code-block:: bash

     $ sudo apt-get update
     $ sudo apt-get install build-essential g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev

*Step 2*. Download the ``cssscl`` package

  .. code-block:: bash 
   
     # use wget 
     $ wget --no-check-certificate https://github.com/oicr-ibc/cssscl/archive/master.tar.gz
     $ tar -zxvf master.tar.gz; mv cssscl-master cssscl 
     
  or use git clone, note that ``sudo apt-get install git`` is required for git access
     
  ..  code-block:: bash  
     
     # use git clone 
     $ git clone https://github.com/oicr-ibc/cssscl.git


*Step 3*. Check that all packages necessary to run the ``cssscl`` are installed and are available by running the ``cssscl_check_pre_installation.sh`` script (only for Ubuntu/Debian distributions). 

  .. code-block:: bash 
    
     $ cd cssscl
     $ ./cssscl_check_pre_installation.sh

**Note:** when prompted follow instructions to export when ``source cssscl/scripts/export.sh`` shows on the screen.

**Note:** for more information regarding the ``cssscl_check_pre_installation.sh`` script see `here <https://github.com/oicr-ibc/cssscl/wiki/cssscl_check_pre_installation>`_.

*Step 4*. In the ``cssscl``  ``directory`` create a virtual environment (e.g. name it ``csssclvenv``)

  .. code-block:: bash 
 
     $ virtualenv csssclvenv


*Step 5*. To begin using the virtual environment, it first needs to be **activated** as shown below:

  .. code-block:: bash 

     $ source csssclvenv/bin/activate


*Step 6*. Install ``cssscl`` as root 

  .. code-block:: bash 

     $ sudo pip install .
    
**Note:** this will install all the python modules necessary for running the ``cssscl`` package in the ``cssscl/csssclvenv/`` directory. 


*Step 7*. Configure ``cssscl``

 .. code-block:: bash 

     $ cssscl configure 
    

Accept all the values prompted by default by pressing [ENTER]  


**Note:** If you are done working in the virtual environment, you can deactivate it as shown below. 

  .. code-block:: bash 

     $ deactivate

If you would like to run the ``cssscl`` program again (and you have deactivated the python virtual environment) you will need to **activate** it again as shown above. 


Option B - Install ``cssscl`` without using the Python's Virtual Environment 
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
     
*Step 1*. Install dependencies on Debian and Ubuntu

Python: Only Python 2.7.3+ is supported. No support for Python 3 at the moment.

In order to compile ``cssscl`` on Debian GNU/Linux 8.1 and Ubuntu 12.04 LTS the following packages need to be installed:

.. code-block:: bash

     $ sudo apt-get update
     $ sudo apt-get install build-essential python2.7 python2.7-dev g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev
       
*Step 2*. Download the ``cssscl`` package 
   
   .. code-block:: bash 

     # use wget 
     $ wget --no-check-certificate https://github.com/oicr-ibc/cssscl/archive/master.tar.gz
     $ tar -zxvf master.tar.gz; mv cssscl-master cssscl 

  or use git clone, note that ``sudo apt-get install git`` is required for git access
     
  ..  code-block:: bash  
     
     # use git clone 
     $ git clone https://github.com/oicr-ibc/cssscl.git



*Step 3*. Check that all packages necessary to run the ``cssscl`` are installed and are avaialble by running the ``cssscl_check_pre_installation.sh`` script (only for Ubuntu/Debian distributions). 
	      
   .. code-block:: bash 

     $ cd cssscl
     $ ./cssscl_check_pre_installation.sh


**Note:** when prompted follow instructions to export when ``source cssscl/scripts/export.sh`` shows on the screen.

**Note:** for more information regarding the ``cssscl_check_pre_installation.sh`` script please see `here <https://github.com/oicr-ibc/cssscl/wiki/cssscl_check_pre_installation>`_.


*Step 4*. Install ``cssscl`` as root  

   .. code-block:: 
   
     $ sudo pip install .        


*Step 5*. Configure ``cssscl`` 

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

*Step 1*. Donwload lzlib-1.5.tar.gz 

.. code-block:: bash 

     $ wget --no-check-certificate http://download.savannah.gnu.org/releases/lzip/lzlib/lzlib-1.5.tar.gz 

*Step 2*. Install lzlib-1.5:

.. code-block:: bash 

     $ gunzip lzlib-1.5.tar.gz
     $ tar -xvf lzlib-1.5.tar
     $ cd lzlib-1.5
     $ ./configure
     $ make
     $ make install


*Step 3*. Donwload Plzip 1.1 

.. code-block:: bash 

     $ wget --no-check-certificate  http://download.savannah.gnu.org/releases/lzip/plzip/plzip-1.1.tar.gz

*Step 4*. Install Plzip

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

MongoDB should be installed using the following set of instructions:

*Ubuntu*

First add the 10gen GPG key, the public gpg key used for signing these packages. It should be possible to import the key into apt's public keyring with a command like this:

.. code-block:: bash 

     $ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv 7F0CEB10

Add this line verbatim to your ``/etc/apt/sources.list``:

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

