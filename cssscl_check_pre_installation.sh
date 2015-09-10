#!/bin/bash

# Copyright 2014-2015, Ivan Borozan <ivan.borozan@gmail.com>
#
# This file is part of the CSSSCL taxonomic sequence classification system.
#
# CSSSCL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CSSSCL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CSSSCL.  If not, see <http://www.gnu.org/licenses/>.


set -u  # Protect against uninitialized vars.
#set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands


vercomp () {
    if [[ $1 == $2 ]]
    then
	local test=0
	echo "$test"
	#return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
	    local test=1
	    echo "$test"
            #return 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
	    local test=2
	    echo "$test"
            #return 2
        fi
    done
    return 0
}


ARCH=$(uname -m | sed 's/x86_//;s/i[3-6]86/32/')

if [ "$ARCH" != "64" ]; then
    echo "You need a 64-bit computer architecture in order to run cssscl. Exiting"
    exit 1
fi


# UBUNTU/DEBIAN 

if [ -f /etc/debian_version ]; then
    OS_VER=($(uname -v))
    if [[ $OS_VER =~ 'Ubuntu' ]]; then
	DISTRO='Ubuntu'
    else
	DISTRO='Debian'
    fi
else
    DISTRO='Other'
fi


if [ $DISTRO == "Other" ]; then
    echo "The automated installation script might not work since your system is neither Debian nor Ubuntu."
    echo "You will need to install BLAST, jellyfish and plzip manually please see the INSTALL.rst and README.rd file for more info."
    exit 1
elif [ $DISTRO == "Debian" ]; then
    echo "Found $DISTRO distribution on the system."
    echo "Please make sure that this user has sudo privilages, if not:"
    echo "su; adduser 'USER_NAME'; then log-out and login again"
    echo "Does the user have sudo privilages? [Y/n]"
    read sudo_user
    if [ "$sudo_user" == "" ] || [ "$sudo_user" == 'Y' -o "$sudo_user" == 'y' ]; then
	echo "OK"
    else
	"Please make sure that the user has sudo privilages. Exiting"
	exit 1	
    fi
else    
    echo "Found $DISTRO distribution on the system."
    echo "Does the user have sudo privilages? [Y/n]"
    read sudo_user
    if [ "$sudo_user" == "" ] || [ "$sudo_user" == 'Y' -o "$sudo_user" == 'y' ]; then
	echo "OK"
	printf "\n"	
    else
	"Please make sure that the user has sudo privilages. Exiting"
	exit 1	
    fi    
fi

# print this warning 

echo  "---------------------------------------------------------------------------------------------------------------------------"
echo  "| If you have either one of the following programs installed : blast, plzip, jellyfish or mongoDB                         |" 
echo  "| that are not in your PATH you need to add them to your executable search PATH variable before you run this script       |"
echo  '| for example : export PATH=$PATH:PATH_TO_YOUR_program                                                                    |'
echo  "| for more details read the instructions in the INSTALL.rst file to continue press [ENTER] to exit press [n]              |"
echo  "---------------------------------------------------------------------------------------------------------------------------" 
read warning
if [ "$warning" == "" ]; then
    echo "OK"
    printf "\n"
elif [ "$warning" == "n" ]; then
    exit 1
fi

pwd=$PWD

# create a log file

if [ ! -f $pwd/log_install.txt ]; then 
    touch $pwd/log_install.txt 
else
    echo "|The log file $pwd/log_install.txt already exists. Checking for the installed cssscl programs.|"
    printf "\n"
    if [ -n "$(type -p blastn)" ] && [ -n "$(type -p jellyfish)" ] && [ -n "$(type -p plzip)" ] && [ -n "$(type -p mongo)" ]; then
	echo "----------------------------------------------------------------------------------------------------"
	echo "| Checked that all the programs necessary for running the cssscl are installed and are in the PATH.|"
	echo "----------------------------------------------------------------------------------------------------"
	exit 0
    else
	status=0
	installed_progs=()
	check_blast=$(grep blastn $pwd/log_install.txt)
	check_jellyfish=$(grep jellyfish $pwd/log_install.txt)
	check_plzip=$(grep plzip $pwd/log_install.txt)
	check_mongo=$(grep mongodb $pwd/log_install.txt)
        # first check which programs were installed 
	if [ "$check_jellyfish" ] && [ -z "$(type -p jellyfish)" ]; then
	   installed_progs+=('jellyfish')
        fi
	if [ "$check_blast" ] && [ -z "$(type -p blastn)" ]; then
	   installed_progs+=('blast') 
	fi
	if [ "$check_plzip" ] && [ -z "$(type -p plzip)" ]; then
	   installed_progs+=('plzip') 
	fi
	installed=${#installed_progs[@]}
        if [ $installed -gt 0 ]; then
	    echo "-----------------------------------------------------"
	    echo "| CHECKED installed programs                        |"
	    echo "-----------------------------------------------------"
	    echo "| The installation log file $pwd/log_install.txt indicates that :"
	    echo "| ${installed_progs[*]} have been installed by this program in $pwd/src/bin but are not currently in the PATH."
	    echo "| NOTE :"
	    echo "| YOU NEED TO ADD $pwd/src/bin to YOUR PATH VARIABLE BY SOURCING THE SCRIPT AS SHOWN BELOW:"
            echo "---------------------------------"
            echo "  source $pwd/scripts/export.sh  "
	    echo "---------------------------------"
        else
	    # re-install 
	    status=1
	fi
	if [ "$check_mongo" ] && [ -z "$(type -p mongo)" ]; then
	    echo "-----------------------------------------------------"
	    echo "| CHECKING mongodb.                                   |"
	    echo "-----------------------------------------------------"
	    echo "| The installation log file $pwd/log_install.txt indicates that :"
	    echo "| mongoDB has been installed but is not currently in the PATH."
	    echo "| Check that mongoDB has been properly installed or/and is in the PATH." 
	fi
	if [ "$status" -eq 0 ]; then 
	    exit 0
	else
	    echo "-----------------------------------"
	    echo "|Will need to re-install programs.|"
	    echo "-----------------------------------"
	fi	
    fi
fi


#echo $pwd

# check that python exists
if type -p python; then
   python_VERSION=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
   if python -c 'import sys; sys.exit(1 if (hex(sys.hexversion)>=0x20703f0 and hex(sys.hexversion)<0x03000000) else 0)'; then   
       echo "Found python version $python_VERSION"
   else
       echo "CSSSCL needs python version 2.7.3+"
       echo "Exiting."
       exit 1
   fi
else
   echo "Python not installed. CSSSCL needs python version 2.7.3+" 
fi


# install the necessary packages for debian/Ubuntu  
python_VERSION2=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')

check_distro_pack=$(grep distro_pac $pwd/log_install.txt)

if [ "$python_VERSION2" == "2.7" ] && [ "$check_distro_pack" == "" ]; then
    sudo apt-get update
    sudo apt-get install build-essential python2.7-dev g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev
    echo "distro_pac" >> $pwd/log_install.txt
elif [ "$python_VERSION2" == "2.8" ] && [ "$check_distro_pack" == "" ]; then
    sudo apt-get update
    sudo apt-get install build-essential python2.8-dev g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev
    echo "distro_pac" >> $pwd/log_install.txt
elif [ "$python_VERSION2" == "2.9" ] && [ "$check_distro_pack" == "" ]; then
    sudo apt-get update
    sudo apt-get install build-essential python2.9-dev g++ libxml2-dev libxslt-dev gfortran libopenblas-dev liblapack-dev
    echo "distro_pac" >> $pwd/log_install.txt
fi


# check that pip is installed 
if type -p pip; then
    pip_VERSION=$(pip --version | awk 'NR==1{print $2}')
    pip_VERSION_EXPECTED="6.0"
    res=($(vercomp $pip_VERSION_EXPECTED $pip_VERSION))
    if [ "$res" == 0 -o "$res" == 2 ]; then
	echo "Found correct version of pip v$pip_VERSION installed"
	# install virtualenv
	echo "Installing virtualenv using pip"
	sudo -H pip install virtualenv
    else 
	echo "Found pip v$pip_VERSION cssscl needs pip v6.+"
	echo "Install pip using sudo? sudo $pwd/scripts/install_pip.sh [Y/n]"
	read pip
	if [ "$pip" == "" ] || [ "$pip" == 'Y' -o "$pip" == 'y' ]; then  
	    $pwd/scripts/install_pip.sh
	else
	    echo "Skipping pip installation and exiting sice pip v6.+ is required"
	    exit 1 
       fi	
    fi
else
    echo "Could not find pip: will need to install pip"
    echo "Install pip using sudo? sudo $pwd/scripts/install_pip.sh [Y/n]"
    read pip
    if [ "$pip" == "" ] || [ "$pip" == 'Y' -o "$pip" == 'y' ]; then 
	$pwd/scripts/install_pip.sh
    else
	echo "Skipping pip installation and exiting sice pip v6.+ is required"
	exit 1 
    fi
fi




# check that blastn is installed 
# first check the log_file
#check_blast=$(grep blastn $pwd/log_install.txt)
#if [ -z "$check_blast" ]; then

if type -p blastn; then
    BLASTN_VERSION=$(blastn -version | awk 'NR==1{print $2}')
    BLASTN_VERSION_EXPECTED="2.2.30"
    if [[ $BLASTN_VERSION =~ '+'$ ]]; then
	#echo "Found blastn v$BLASTN_VERSION in the PATH"
	BLASTN_VERSION=${BLASTN_VERSION%+}
	res=($(vercomp $BLASTN_VERSION_EXPECTED $BLASTN_VERSION))
	if [ "$res" == 1 ]; then
	   echo "Found blastn v$BLASTN_VERSION in the PATH" 
	   echo "Will install v$BLASTN_VERSION_EXPECTED in $pwd/src/bin [Y/n]" 
	   read blastn
	   if [ "$blastn" == "" ] || [ "$blastn" == 'Y' -o "$blastn" == 'y' ]; then  
	       $pwd/scripts/install_blast.sh
	       echo "blastn" >> $pwd/log_install.txt
	   else
	       echo "Skipping blastn installation."
	       #exit 1 
	   fi
	else
	    echo "Found correct version of blastn v$BLASTN_VERSION installed"
	fi
    else
	echo "Found legacy blastn version" 
	echo "cssscl requires  blastn v2.2.30+" 
	echo "Will install v$BLASTN_VERSION v2.2.30+ $pwd/src/bin [Y/n]"
	read blastn
	if [ "$blastn" == "" ] || [ "$blastn" == 'Y' -o "$blastn" == 'y' ]; then  
	    $pwd/scripts/install_blast.sh 
	    echo "blastn" >> $pwd/log_install.txt
	else
	   echo "Skipping blastn installation and exiting sice at at least blast+/2.2.30 is required"
	   exit 1 
	fi    
    fi
else
   echo "cssscl requires at least blast+/2.2.30" 
   echo "blastn not detected will install blastn v2.2.30+ in $pwd/src/bin [Y/n]"
   read blastn
   if [ "$blastn" == "" ] || [ "$blastn" == 'Y' -o "$blastn" == 'y' ]; then  
       $pwd/scripts/install_blast.sh 
       echo "blastn" >> $pwd/log_install.txt
   else
	   echo "Skipping blastn installation and exiting since at least blast+/2.2.30 is required"
	   exit 1 
   fi       
fi



# check that jellyfish is installed 
# first check the log_file
#check_jellyfish=$(grep jellyfish $pwd/log_install.txt)
#if [ -z "$check_jellyfish" ]; then

if type -p jellyfish; then
    JELLYFISH_VERSION=$(jellyfish --version | awk '{print $2}')
    JELLYFISH_VERSION_EXPECTED="1.1"
    res=($(vercomp $JELLYFISH_VERSION_EXPECTED $JELLYFISH_VERSION))
    if [[ $JELLYFISH_VERSION =~ ^1\. ]]
    then
        if [ "$res" == 1 ]; then
	    echo "Found jellyfish v$JELLYFISH_VERSION"
	    echo "Will install v$JELLYFISH_VERSION_EXPECTED in $pwd/src/bin [Y/n]" 
	    read jellyfish
	    if [ "$jellyfish" == "" ] || [ "$jellyfish" == 'Y' -o "$jellyfish" == 'y' ]; then  
		if type -p g++; then
		    echo "Found the g++ compiler"
		else
		    echo "The g++ compiler not found"
		    echo "Installing the g++ compiler using sudo? [Y/n]"
		    read g2p
		    if [ "$g2p" == "" ] || [ "$g2p" == 'Y' -o "$g2p" == 'y' ]; then
			sudo apt-get install g++
		    fi
		fi
		$pwd/scripts/install_jellyfish.sh
		echo "jellyfish" >> $pwd/log_install.txt
	    else
		echo "Skipping jellyfish installation and exiting sice cssscl requires jellyfish 1.1+ to be installed"
		exit 1 
	    fi
	else
	    echo "Found correct version of jellyfish v$JELLYFISH_VERSION installed"
	fi
    else
	echo "Found jellyfish v$JELLYFISH_VERSION in PATH"
	echo "cssscl requires jellyfish version 1.1+"
	echo "Will need to install jellyfish version 1.1 in $pwd/src/bin [Y/n]"
	read jellyfish
	if [ "$jellyfish" == "" ] || [ "$jellyfish" == 'Y' -o "$jellyfish" == 'y' ]; then
	    if type -p g++; then
		echo "Found the g++ compiler"
	    else
		echo "The g++ compiler not found"
		echo "Installing the g++ compiler using sudo? [Y/n]"
		read g2p
		if [ "$g2p" == "" ] || [ "$g2p" == 'Y' -o "$g2p" == 'y' ]; then
		    sudo apt-get install g++
		fi
	    fi
	    $pwd/scripts/install_jellyfish.sh 
	    echo "jellyfish" >> $pwd/log_install.txt
	else
	   echo "Skipping jellyfish installation and exiting sice cssscl requires jellyfish 1.1+ to be installed"
	   exit 1 
	fi       	    
    fi
else
    echo "cssscl requires jellyfish version 1.1+ to be installed"
    echo "Will need to install jellyfish version 1.1 in $pwd/src/bin [Y/n]"
    read jellyfish
    if [ "$jellyfish" == "" ] || [ $jellyfish == 'Y' -o $jellyfish == 'y' ]; then
	if type -p g++; then
	    echo "Found the g++ compiler"
	else
	    echo "The g++ compiler not found"
	    echo "Installing the g++ compiler using sudo? [Y/n]"
	    read g2p
	    if [ "$g2p" == "" ] || [ "$g2p" == 'Y' -o "$g2p" == 'y' ]; then
		sudo apt-get install g++
	    fi
	fi
	$pwd/scripts/install_jellyfish.sh 
	echo "jellyfish" >> $pwd/log_install.txt
    else
       echo "Skipping jellyfish installation and exiting sice cssscl requires jellyfish 1.1+ to be installed"
       exit 1 
    fi
fi


# check that plzip is installed 
#check_plzip=$(grep plzip $pwd/log_install.txt)
#if [ -z "$check_plzip" ]; then

if type -p plzip; then
    PLZIP_VERSION=$(plzip --version | awk 'NR==1{print $2}')
    if [[ $PLZIP_VERSION =~ 1.1$ ]];then
	echo "Found correct version of plzip v$PLZIP_VERSION installed"
    else
	echo "Found plzip v$PLZIP_VERSION"
	echo "cssscl requires plzip version 1.1"
	if type -p g++; then
	    echo "Found the g++ compiler"
	else
	    echo "The g++ compiler not found"
	    echo "Installing the g++ compiler using sudo? [Y/n]"
	    read g2p
	    if [ "$g2p" == "" ] || [ "$g2p" == 'Y' -o "$g2p" == 'y' ]; then
		sudo apt-get install g++
	    fi
	fi
	echo "Installing plzip version 1.1"
	$pwd/scripts/install_plzip.sh
	echo "plzip" >> $pwd/log_install.txt
    fi
else
    echo "CSSSCL requires plzip version 1.1 to be installed"
    echo "Installing plzip version 1.1 in $pwd/src/bin"
    if type -p g++; then
	echo "Found the g++ compiler"
    else
	echo "The g++ compiler not found"
	echo "Installing the g++ compiler using sudo? [Y/n]"
	read g2p
	if [ "$g2p" == "" ] || [ "$g2p" == 'Y' -o "$g2p" == 'y' ]; then
	    sudo apt-get install g++
	fi
    fi
    $pwd/scripts/install_plzip.sh 
    echo "plzip" >> $pwd/log_install.txt
fi


# check that mongodb is installed 
if type -p mongo; then
    mongo_VERSION=$(mongo --version | awk '{print $4}')
    mongo_VERSION_EXPECTED="2.4"
    res=($(vercomp $mongo_VERSION_EXPECTED $mongo_VERSION))
    if [[ $mongo_VERSION =~ ^2\. ]];then
	    echo "Found mongodb v$mongo_VERSION"
    else
	echo "Found mongo v$mono_VERSION"
	echo "cssscl requires mongo v2.+ to be installed"
	echo "Please install a new version of mongodb v.2+"
    fi
else
    echo "No mongodb found"
    echo "cssscl requires mongo v2.+ to be installed"
    if [ $DISTRO == 'Other' ]; then
	echo "For now only DEBIAN/UBUNTU automated installation supported. Please install mongodb v2.4+ (and < v3.0) manually on your system"
	echo "Exiting"
	exit 1 
    fi
    echo "Will try to install mongodb on your system"
    echo "Install mongodb using sudo? [Y/n]"
    read mongodb
    if [ "$mongodb" == "Y" -o "$mongodb" == "y" ] || [ "$mongodb" == "" ] && [ "$DISTRO" == 'Ubuntu' ];then
	   sudo $pwd/scripts/install_mongodb_ubuntu.sh
	   echo "mongodb" >> $pwd/log_install.txt
    elif [ "$mongodb" == "Y" -o "$mongodb" == "y" ] || [ "$mongodb" == "" ] && [ "$DISTRO" == 'Debian' ];then
	sudo $pwd/scripts/install_mongodb_debian.sh
	echo "mongodb" >> $pwd/log_install.txt
    else
	echo "Skipping mongodb installation and exiting sice cssscl requires mongodb to be installed" 	
	exit 1
    fi
fi


check_path=$(grep path $pwd/log_install.txt)

# check if any of these were installed
check_blast=$(grep blastn $pwd/log_install.txt)
check_jellyfish=$(grep jellyfish $pwd/log_install.txt)
check_plzip=$(grep plzip $pwd/log_install.txt)
check_mongo=$(grep mongodb $pwd/log_install.txt)

# add "path" the log file if either of the above (blastn...etc) installed 
if [ -z "$check_path" ]; then
    if [ "$check_jellyfish" ] || [ "$check_blast" ] || [ "$check_plzip" ]; then 
	echo "path" >> $pwd/log_install.txt
	if [ -f $HOME/.bashrc ]; then
	    echo "#adding to the PATH in $HOME/.bashrc for the cssscl programs:" | tee -a $HOME/.bashrc
	    echo "export PATH='$pwd/src/bin:$PATH'" | tee -a $HOME/.bashrc
	fi
	if [ -f $HOME/.profile ]; then
	    echo "#adding to the PATH in $HOME.profile for the cssscl programs:" | tee -a $HOME/.profile
	    echo "export PATH='$pwd/src/bin:$PATH'" | tee -a $HOME/.profile
	fi
	eval '$pwd/scripts/export.sh'	
    fi
fi
 

if [ -n "$(type -p blastn)" ] && [ -n "$(type -p jellyfish)" ] && [ -n "$(type -p plzip)" ] && [ -n "$(type -p mongo)" ]; then
    echo "----------------------------------------------------------------------------------------------------"
    echo "| Checked that all the programs necessary for running the cssscl are installed and are in the PATH.|"
    echo "----------------------------------------------------------------------------------------------------"
    exit 0
else
    # first check the log file for programs that were installed 
    if [ -e $pwd/log_install.txt ]; then
	installed_progs=()
	if [ "$check_jellyfish" ] && [ -z "$(type -p jellyfish)" ]; then
	   installed_progs+=('jellyfish')
        fi
	if [ "$check_blast" ] && [ -z "$(type -p blastn)" ]; then
	   installed_progs+=('blast') 
	fi
	if [ "$check_plzip" ] && [ -z "$(type -p plzip)" ]; then
	    installed_progs+=('plzip') 
	fi
	installed=${#installed_progs[@]}
        if [ $installed -gt 0 ]; then
	    echo "-----------------------------------------------------"
	    echo "| CHECKING installed programs                       |"
	    echo "-----------------------------------------------------"
	    echo "| The installation log file $pwd/log_install.txt indicates that :"
	    echo "| ${installed_progs[*]} have been installed by this program in $pwd/src/bin but are not currently in the PATH."
	    echo "| YOU NEED TO ADD $pwd/src/bin to YOUR PATH VARIABLE BY SOURCING THE SCRIPT AS SHOWN BELOW:"
            echo "---------------------------------"
            echo "| source $pwd/scripts/export.sh  "
	    echo "---------------------------------"
	fi
	if [ "$check_mongo" ] && [ -z "$(type -p mongo)" ]; then
	    echo "-----------------------------------------------------"
	    echo "| CHECKING mongodb.                                   |"
	    echo "-----------------------------------------------------"
	    echo "| The installation log file $pwd/log_install.txt indicates that :"
	    echo "| mongo has been installed but is not currently in the PATH."
	    echo "| Check that mongoDB has been properly installed or/and is in the PATH." 
	fi     
    else
	echo "---------------------------------------------------------"
	echo "|Can not find the log installation file in $pwd. Exiting.|"
	echo "---------------------------------------------------------"
	exit 1	
    fi
fi