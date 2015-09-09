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

apt-key adv --keyserver keyserver.ubuntu.com --recv 7F0CEB10
#echo 'deb http://downloads-distro.mongodb.org/repo/debian-sysvinit dist 10gen' | tee /etc/apt/sources.list.d/mongodb.list
echo 'deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen' | tee -a /etc/apt/sources.list
apt-get update 
apt-get install mongodb-10gen=2.4.14