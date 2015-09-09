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


#if type -p pip; then
#    echo found pip executable in PATH
#else
#echo "no pip: will need to install pip using sudo"
wget --no-check-certificate -P /tmp/ https://bootstrap.pypa.io/get-pip.py
echo "If setuptools (or distribute) is not already installed, get-pip.py will install setuptools for you."
sudo -H python /tmp/get-pip.py
# install virtualenv
sudo -H pip install virtualenv
#fi

