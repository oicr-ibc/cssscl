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
# RUN THIS SCRIPT with sudo.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

pwd=$PWD
#echo $pwd

echo "Donwloading lzlib-1.5.tar from http://download.savannah.gnu.org/releases/lzip/lzlib/lzlib-1.5.tar.gz"
echo "Will need to install lzlib-1.5 using sudo: [Y/n]"
read lzlib
if [ "$lzlib" == "" ] || [ "$lzlib" == 'Y' -o "$lzlib" == 'y' ]; then
	    suffix='scripts'
	    install_dir=${pwd%$suffix}
            mkdir -p $install_dir/src
	    #wget --no-check-certificate -P /tmp/ http://download.savannah.gnu.org/releases/lzip/lzlib/lzlib-1.5.tar.gz
	    wget --no-check-certificate -P /tmp/ https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/source/lzlib-1.5.tar.gz
	    tar -xzf /tmp/lzlib-1.5.tar.gz -C $install_dir/src
	    cd $install_dir/src/lzlib-1.5
	    ./configure 
	    make
	    sudo make install
	    cd $pwd
	    #echo "Downloading plzip1.1 from http://download.savannah.gnu.org/releases/lzip/plzip/plzip-1.1.tar.gz"
	    echo "Downloading plzip1.1 from https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/source/plzip-1.1.tar.gz"
	    #wget --no-check-certificate -P /tmp/ http://download.savannah.gnu.org/releases/lzip/plzip/plzip-1.1.tar.gz
	    wget --no-check-certificate -P /tmp/ https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/source/plzip-1.1.tar.gz
	    tar -xzf /tmp/plzip-1.1.tar.gz -C $install_dir/src
	    cd $install_dir/src/plzip-1.1
	    #./configure --prefix=$HOME --datarootdir=$install_dir/src/plzip-1.1  --infodir=$install_dir/src/plzip-1.1 --mandir=$install_dir/src/plzip-1.1
	    ./configure --prefix=$install_dir/src --datarootdir=$install_dir/src/plzip-1.1  --infodir=$install_dir/src/plzip-1.1 --mandir=$install_dir/src/plzip-1.1
	    make
	    make install
	else
	    exit 1
fi