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

echo "Downloading BLAST from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz"
#echo "Downloading BLAST from https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/source/ncbi-blast-2.2.30+-x64-linux.tar.gz"
suffix='scripts'
install_dir=${pwd%$suffix}
#echo $install_dir/bin/
mkdir -p $install_dir/src/bin
wget --no-check-certificate -P /tmp/ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz
#wget --no-check-certificate -P /tmp/ https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/source/ncbi-blast-2.2.30+-x64-linux.tar.gz
tar -xzf /tmp/ncbi-blast-2.2.30+-x64-linux.tar.gz -C $install_dir/src
#mkdir -p $HOME/bin/
#cp $install_dir/src/ncbi-blast-2.2.30+/bin/* $HOME/bin/
cp $install_dir/src/ncbi-blast-2.2.30+/bin/* $install_dir/src/bin/
