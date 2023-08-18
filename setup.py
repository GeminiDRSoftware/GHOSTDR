#!/usr/bin/env python3

"""
Setup script for ghostdr

Based on the setup.py from the gemini_python distribution.

In this package:
    ghostdr           : Recipes and Primitives for GHOST
    ghost_instruments : AstroData subclass for GHOST
    

Usage:
  python setup.py install
"""

import os.path
import os
import re
import glob
import sys

from setuptools import setup, find_packages

ADLIB_PACKAGES = ['ghostdr']
gitdir = re.compile('.git')
fitsfile = re.compile('.fits$')
txtfile = re.compile('.txt$')
dotpy = re.compile('.py$')
dotsh = re.compile('.sh$')

# PACKAGE_DATA
PACKAGE_DATA = {}

for p in ADLIB_PACKAGES:
    PACKAGE_DATA[p] = []
    fitslist = []
    txtlist = []
    for root, dirs, files in os.walk(os.path.join(p, 'ghost', 'lookups')):
        if not gitdir.search(root) and len(files) > 0:
            fitsfiles = [f for f in files if fitsfile.search(f)]
            dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
            PACKAGE_DATA[p].extend( map((lambda f: os.path.join(dest, f)), fitsfiles) )
            fitslist.extend( map((lambda f: os.path.join(dest, f)), fitsfiles) )

            txtfiles = [f for f in files if txtfile.search(f)]
            dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
            PACKAGE_DATA[p].extend( map((lambda f: os.path.join(dest, f)), txtfiles) )
            txtlist.extend( map((lambda f: os.path.join(dest, f)), txtfiles) )


# SCRIPTS
SCRIPTS = [os.path.join('simulator', 'testsim.py')]
for root, dirs, files in os.walk('utils'):
    if not gitdir.search(root) and len(files) > 0:
        files = [f for f in files if dotpy.search(f) or dotsh.search(f)]
        dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
        SCRIPTS.extend( map((lambda f: os.path.join('utils', dest, f)), files) )

# PACKAGES
PACKAGES = find_packages(where='.', exclude=['ghostdr.ghost.test', 'ghostdr.ghost.recipes.test'])
PACKAGES.append('ghostdr.ghost.lookups.BPM')
PACKAGE_DIR = {"": "."}
SIMPACKAGES = find_packages(where=os.path.join('.', 'simulator'))
PACKAGES.extend(SIMPACKAGES)
for p in SIMPACKAGES:
    PACKAGE_DIR[p] = os.path.join('.', 'simulator', p)
    PACKAGE_DATA[p] = []
    for root, dirs, files in os.walk(os.path.join('simulator', p, 'data')):
        if not gitdir.search(root) and len(files) > 0:
            dest = root.split('/',2)[2] if len(root.split('/',2)) > 2 else ""
            PACKAGE_DATA[p].extend( map((lambda f: os.path.join(dest, f)), files) )

setup(
    name='ghostdr',
    version='1.0.0',
    description='GHOST Data Reduction',
    author="ANU & Gemini Observatory",
    url="https://ghost-drtutorial.readthedocs.io/en/release-3.0.x/index.html",
    #author='RSAA Software Engineering',
    #author_email='ghostdr@mso.anu.edu.au',
    #url='http://rsaa.anu.edu.au',
    maintainer='SUSD',
    classifiers=[
        'Development Status :: Alpha',
        'Intended Audience :: Alpha Testers',
        'Operating System :: Linux :: RHEL',
        'Programming Language :: Python',
        'Topic :: Gemini',
        'Topic :: Data Reduction',
        'Topic :: Astronomy',
    ],
    packages=PACKAGES,
    package_dir=PACKAGE_DIR,
    package_data=PACKAGE_DATA,
    scripts=SCRIPTS,
)
