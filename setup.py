# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import re

from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.

__version__ = '%s'
"""

long_description = """
skgGenomeTracks aims to produce high-quality genome browser tracks that
are highly customizable. Currently, it is possible to plot:

 * bigwig
 * bed (many options)
 * bedgraph
 * links (represented as arcs)
 * Hi-C matrices (if [HiCExplorer](http://hicexplorer.readthedocs.io) is installed)
 * DNA methylation

skgGenomeTracks can make plots with or without Hi-C data. The following is an example output of skggenometracks from [Ramírez et al. 2017](https://www.nature.com/articles/s41467-017-02525-w)
"""


def update_version_py():
    if not os.path.isdir(".git"):
        print("This does not appear to be a Git repository.")
        return
    try:
        p = subprocess.Popen(["git", "describe",
                              "--tags", "--always"],
                             stdout=subprocess.PIPE)
    except EnvironmentError:
        print("unable to run git, leaving skggenometracks/_version.py alone")
        return
    stdout = p.communicate()[0]
    if p.returncode != 0:
        print("unable to run git, leaving skggenometracks/_version.py alone")
        return
    ver = stdout.strip()
    f = open(os.path.join("skggenometracks", "_version.py"), "w")
    f.write(VERSION_PY % ver)
    f.close()
    print("set skggenometracks/_version.py to '%s'" % ver)


def get_version():
    try:
        f = open(os.path.join("skggenometracks", "_version.py"))
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


class sdist(_sdist):

    def run(self):
        # update_version_py()
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)

# Install class to check for external dependencies from OS environment


class install(_install):

    def run(self):
        # update_version_py()
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return

    def checkProgramIsInstalled(self, program, args, where_to_download,
                                affected_tools):
        try:
            subprocess.Popen([program, args],
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE)
            return True
        except EnvironmentError:
            # handle file not found error.
            # the config file is installed in:
            msg = "\n**{0} not found. This " \
                  "program is needed for the following "\
                  "tools to work properly:\n"\
                  " {1}\n"\
                  "{0} can be downloaded from here:\n " \
                  " {2}\n".format(program, affected_tools,
                                  where_to_download)
            sys.stderr.write(msg)

        except Exception as e:
            sys.stderr.write("Error: {}".format(e))


install_requires_py = ["numpy >=1.16",
                       "matplotlib ==3.1.1",
                       "intervaltree >=2.1.0",
                       "pyBigWig >=0.3.4",
                       "future >=0.17.0",
                       "hicmatrix >=9",
                       "pysam >=0.14",
                       "pytest",
                       "gffutils >=0.9"
                       ]

if sys.version_info[0] == 2 or (sys.version_info[0] == 3 and sys.version_info[1] == 4):
    install_requires_py.append("configparser >= 3.5.0")

setup(
    name='skgGenomeTracks',
    version=get_version(),
    author='Fidel Ramírez, Vivek Bhardwaj, Joachim Wolf, Björn Grüning',
    author_email='deeptools@googlegroups.com',
    packages=find_packages(exclude=['tests']),
    scripts=['bin/skg_make_tracks_file', 'bin/skgGenomeTracks',
             'bin/sgt', 'bin/gene2bed', 'bin/bismark2mr'],
    include_package_data=True,
    package_dir={'skggenometracks': 'skggenometracks'},
    url='http://skggenometracks.readthedocs.io',
    license='LICENSE.txt',
    description='Command-line tool to make beautiful and reproducible genome browser snapshots',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=install_requires_py,
    zip_safe=False,
    python_requires='>=3.6.*, <4',
    cmdclass={'sdist': sdist, 'install': install}
)
