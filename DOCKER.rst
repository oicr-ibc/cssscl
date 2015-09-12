We have setup a Dockerfile to create an image and a container that runs ubuntu 12.04 provisioned with Python 2.7.3, MongoDB, BLAST, plzip, jellyfish and the ``cssscl`` program for the quick deployment and testing.

**Procedure:**

  | a. Make sure that you have Docker installed on your system.
  |    On how to install Docker on your system please consult the docker installation guide `here <https://docs.docker.com/installation/>`_.
  | b. Download the `ubuntucsss.tar.gz <https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/code_2xx/ubuntucsss.tar.gz>`_
     file that contains the Dockerfile.
  | c. tar -zxvf ubuntucsss.tar.gz
  | d. To run the ``cssscl`` program follow the instructions below:


**Use docker to run cssscl**

1. If necessary add a docker user

    .. code-block:: bash 
      
      $ sudo adduser -a [username] docker
      # then Logout

2.  Login and ``cd`` to the ``ubuntucsss`` ``directory``

    .. code-block:: bash 

      $ cd ./ubuntucsss   

3. Build the docker image using the Dockerfile located in the ubuntucss directory as shown below

    .. code-block:: bash 

      $ docker build -t cssscl/ubuntucsss .

   or if you did not add the docker user run it with the sudo command 

    .. code-block:: bash 

      $ sudo docker build -t cssscl/ubuntucsss .


4. Now use the docker image cssscl/ubuntucsss to create and run a container:

    .. code-block:: bash 
   
      $ docker run -ti cssscl/ubuntucsss /bin/bash       

   or if you did not add the docker user run it with the sudo command 

    .. code-block:: bash 

      $ sudo docker run -ti cssscl/ubuntucsss /bin/bash 

   Note: you could setup any number of cpus for the container to use as shown below:

    .. code-block:: bash 
    
      $ docker run -ti --cpuset-cpus="0-7" cssscl/ubuntucsss /bin/bash       

   or if you did not add the docker user run it with the sudo command 

    .. code-block:: bash 

      $ sudo docker run -ti --cpuset-cpus="0-7" cssscl/ubuntucsss /bin/bash

   This will run the container with 8 cpus

5. Inside the running container start the mongo database as shown below:

    .. code-block:: bash 
    
      $ mongod --fork --logpath /data/db/log


6. Configure cssscl :

    .. code-block:: bash 

      $ cssscl configure 


Accept all the values prompted by default by pressing [ENTER]  
 

**Run the** ``cssscl`` **classifier**


7. Build the necessary databases from the training set:

    .. code-block:: bash 

      $ cssscl build_dbs -btax -c -blast -nt 2 /home/test_data/TRAIN.fa /home/taxon/

(the whole process should take ~ 37 min using 2 CPUs)

By default all databases will be outputted to the ``directory`` where the train.fa resides (note that all paths provided need to be absolute/full paths to the files/directories).

For more information about the ``cssscl`` ``build_dbs`` please consult its help page by typing:

    .. code-block:: bash

      $ cssscl build_dbs --help

8. Perform the classification using ``cssscl`` :

    .. code-block:: bash

      $ cssscl classify -c -blast blastn -tax genus -nt 2 /home/test_data/test/TEST.fa /home/test_data/

(the whole process should take ~ 29 min using 2 CPUs)

**Note**: in the above example the output file ``cssscl_results_genus.txt`` with classification results will be located in the directory where the TEST.fa resides. 

Note that for the test set data the parameters of the model have already been optimized and are included as part of the test set data, thus optimization is not required to be performed prior to running the classifier.

This will run the classifier with all the similarity measures (including the compression and the blast measure) described in Borozan et al. *"Integrating alignment-based and alignment-free sequence similarity measures for biological sequence classification."*  Bioinformatics. 2015 Jan 7. pii: btv006. 

For more information about the ``cssscl`` ``classify`` please consult its help page by typing: 
 
    .. code-block:: bash 

      $ cssscl classify --help 
