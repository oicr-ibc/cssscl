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

pwd=$PWD

echo "Do you really want to uninstall the cssscl package [Y/n]?"
read uninstall
if [ "$uninstall" == "" ] || [ "$uninstall" == 'Y' -o "$uninstall" == 'y' ]; then
    echo "Removing the cssscl package using pip." 
    sudo -H pip uninstall cssscl 
    echo "Removing the $pwd directory." 
    rm -rf $pwd
    echo "Please remove the changes to the PATH variable associated with the cssscl package in the $HOME/.bashrc and/or $HOME/.profile files."
else
    exit 1
fi