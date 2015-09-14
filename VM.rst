**Note:** If you are testing the ``cssscl`` package using this Oracle VirtualBox VM, please make sure that you have at least 1024 MB of RAM.

We have setup an OVF-formatted virtual machine (VM) running Ubuntu provisioned with Python 2.7.3 (including all the python modules), MongoDB, BLAST, plzip and jellyfish for the quick testing of the ``cssscl`` program. The VM also includes taxon and test data.

**Procedure for setting up the VM**

   | 1. Download the .ova file from `here <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/cssscl_opt.ova>`_
   | 2. In Workstation, open VirtualBox software, select File > Import appliance.
   | 3. Browse to the .ova file and click Open.
   | 4. Follow the instructions, and click Import. Workstation performs OVF specification conformance and virtual hardware compliance checks. A status bar indicates the progress of the import process.
   | 5. If the import fails, click Retry to try again, or click Cancel to cancel the import.
   | If you retry the import, Workstation relaxes the OVF specification conformance and virtual hardware compliance checks and you might not be able to use the virtual machine in Workstation.
   | After Workstation successfully imports the OVF virtual machine, the virtual machine appears in the virtual machine library.
   | 6. Then to install and run ``cssscl`` follow the instructions below.

**Run cssscl on a VM**

First login to the VM

   | precise64 login:vagrant
   | Password:vagrant


1. Download the ``cssscl`` package and install the program 

  .. code-block:: bash 
    
    # download the cssscl package
    $ wget --no-check-certificate https://github.com/oicr-ibc/cssscl/archive/master.tar.gz
    $ tar -zxvf master.tar.gz; mv cssscl-master cssscl 
    $ cd cssscl/
    # check the installation 
    $ ./cssscl_check_pre_installation.sh
    # install
    $ sudo pip install .

2. Configure ``cssscl`` :

  .. code-block:: bash 

    $ cssscl configure 

Accept all the values prompted by default by pressing [ENTER]  
 

3. Build the necessary databases from the training set:

  .. code-block:: bash 

    $ cssscl build_dbs -btax -c -blast -nt 2 /home/vagrant/test_data/TRAIN.fa /home/vagrant/taxon/

(the whole process should take ~ 37 min using 2 CPUs)

By default all databases will be outputted to the ``directory`` where the train.fa resides.

  .. code-block:: bash 

    $ cssscl build_dbs --help

For more information about cssscl please read the REAME.rst or the INSTALL.rst provided with this package.


4. Perform the classification using ``cssscl``:

  .. code-block:: bash 

    $ cssscl classify -c -blast blastn -tax genus -nt 2 /home/vagrant/test_data/test/TEST.fa /home/vagrant/test_data/


(the whole process should take ~ 29 min using 2 CPUs)

Note that in the above example the output file ``cssscl_results_genus.txt`` with classification results will be located in the directory where the TEST.fa resides. 

Note that for the test set data the parameters of the model have already been optimized and are included as part of the test set data, thus optimization is not required to be performed prior to running the classifier.

This will run the classifier with all the similarity measures (including the compression and the blast measure) described in Borozan et al. "Integrating alignment-based and alignment-free sequence similarity measures for biological sequence classification."  Bioinformatics. 2015 Jan 7. pii: btv006. 

For more information about the cssscl ``classify`` please consult its help page by typing: 

  .. code-block:: bash 

    $ cssscl classify --help 

For more information about ''cssscl'' please read the README.rst and INSTALL.rst provided with this package.


