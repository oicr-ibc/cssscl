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
#echo "Donwloading jellyfish from https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/source/jellyfish-1.1.11.tar.gz"
echo "Donwloading jellyfish from https://github.com/gmarcais/Jellyfish/releases/download/v1.1.12/jellyfish-1.1.12.tar.gz"
suffix='scripts'
install_dir=${pwd%$suffix}
mkdir -p $install_dir/src
wget --no-check-certificate -P /tmp/ https://github.com/gmarcais/Jellyfish/releases/download/v1.1.12/jellyfish-1.1.12.tar.gz
#wget --no-check-certificate -P /tmp/ https://collaborators.oicr.on.ca/vferretti/borozan_cssscl/source/jellyfish-1.1.11.tar.gz
tar -xzf /tmp/jellyfish-1.1.12.tar.gz -C $install_dir/src
cd $install_dir/src/jellyfish-1.1.12
./configure --prefix=$install_dir/src
make
make install
#mkdir -p $HOME/bin/
#ln -s $install_dir/src/bin/jellyfish $HOME/bin/jellyfish
#ln -s $install_dir/src/bin/jellyfish $pwd/bin/jellyfish
