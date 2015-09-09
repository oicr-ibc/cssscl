#!/usr/bin/env python


# Copyright 2015(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.i


from setuptools import setup
import sys, os


if sys.version_info < (2, 7, 3):
    exec('raise Error, "Python 2.7.3 or later is required"')


def read(*path):
        return open(os.path.join(os.path.abspath(os.path.dirname(__file__)), *path)).read()


VERSION = '1.0'
README = read('README.rst')
NEWS = read('NEWS.rst')
install_requires = ['cython', 'numpy', 'pymongo>=2.8,<3.0', 'biopython', 'scikit-learn', 'scipy']

if sys.version_info < (2, 7, 3):
    install_requires.append('argparse')

# Get classifiers from http://pypi.python.org/pypi?%3Aaction=list_classifiers
classifiers = """
Development Status :: 5 - Production/Stable
License :: OSI Approved :: GNU General Public License (GPL)
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
Programming Language :: Python :: 2.7
Topic :: Scientific/Engineering :: Bio-Informatics
Operating System :: Unix
"""

config = {
    'name': 'cssscl',
    'version': VERSION,
    'description': 'Combining sequence similarity scores for biological sequence classification',
    'long_description': README + '\n\n' + NEWS,
    'license': 'GNU General Public License, Version 3.0',
    'author': 'Ivan Borozan',
    'author_email': 'ivan.borozan@gmail.com',
    #'url': 'https://github.com/cssscl/cssscl',
    #'download_url': 'https://github.com/cssscl/cssscl',
    'classifiers': filter(None, classifiers.split("\n")),
    'scripts': ['bin/cssscl'],
    'packages': ['cssscl'],
    'zip_safe': True,
   'install_requires': install_requires
}


def setup_package():
    """Setup Package"""

    setup(**config)


if __name__ == '__main__':
    setup_package()
