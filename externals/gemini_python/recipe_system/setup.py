#!/usr/bin/env python

"""
Setup script for recipe_system.

The tools and modules in this package are developed by the Gemini Data Processing Software Group.

In this module:
    The Recipe System

Usage:
  python setup.py install --prefix=/astro/iraf/i686/gempylocal
  python setup.py sdist
"""

import os.path
import re
import sys

from distutils.core import setup

svndir = re.compile('.svn')

MODULENAME = 'recipe_system'

# PACKAGES and PACKAGE_DIRS
SUBMODULES = ['adcc',
              'adcc.servers',
              'cal_service',
              'reduction']
PACKAGES = [MODULENAME]
for m in SUBMODULES:
    PACKAGES.append('.'.join([MODULENAME,m]))
PACKAGE_DIRS = {}
PACKAGE_DIRS[MODULENAME] = '.'

# PACKAGE_DATA
PACKAGE_DATA = {}
PACKAGE_DATA[MODULENAME] = []
# reactivate this once ReleaseNote, README, INSTALL have content.
#for s in ['.']+SUBMODULES:
#    PACKAGE_DATA[MODULENAME].extend([os.path.join(s,'ReleaseNote'),
#                                     os.path.join(s,'README'),
#                                     os.path.join(s,'INSTALL'),
#                                     ])

for root, dirs, files in os.walk(os.path.join('adcc', 'client', 'adcc_faceplate')):
    if not svndir.search(root) and len(files) > 0:
        PACKAGE_DATA[MODULENAME].extend( map((lambda f: os.path.join(root, f)), files) )
 
# We don't have unit tests for the RS yet.       
#for root, dirs, files in os.walk('tests'):
#    if not svndir.search(root) and len(files) > 0:
#        PACKAGE_DATA[MODULENAME].extend( map((lambda f: os.path.join(root, f)), files) )


# DATA_DIRS and DATA_FILES
# no doc yet
DATA_FILES = []
#DOC_DIR = os.path.join('share','recipe_system')
#for root, dirs, files in os.walk('doc'):
#    if not svndir.search(root) and len(files) > 0:
#        dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
#        DOC_FILES = map((lambda f: os.path.join(root,f)), files)      
#        DATA_FILES.append( (os.path.join(DOC_DIR,dest), DOC_FILES) )

# SCRIPTS
RS_SCRIPTS = [ os.path.join('apps','adcc'),
               os.path.join('apps','listprimitives'),
               os.path.join('apps','reduce'),
             ]
if "sdist" in sys.argv:
    #ASTRODATA_SCRIPTS contains the name of the links which might not be dereferenced during sdist
    #Therefore, here we package the .py those links point to.  During "install" the links are
    #dereferenced, always, as far as I can tell, so there's no need for the .py then.
    PYFILES = []
    dotpy = re.compile(".py$")
    for script in RS_SCRIPTS:
        if not dotpy.match(script):
            PYFILES.append(''.join([script,'.py']))
    RS_SCRIPTS.extend(PYFILES)

#SOMEOTHER_SCRIPTS.extend(['another', 'script'])
SCRIPTS = []
SCRIPTS.extend(RS_SCRIPTS)

EXTENSIONS = None

setup ( name='recipe_system',
        version='0.2.0',
        description='AstroData',
        author='Gemini Data Processing Software Group',
        url='http://www.gemini.edu',
        maintainer='Gemini Data Processing Software Group',
        packages=PACKAGES,
        package_dir=PACKAGE_DIRS,
        package_data=PACKAGE_DATA,
        data_files=DATA_FILES,
        scripts=SCRIPTS,
        ext_modules=EXTENSIONS,
        classifiers=[
            'Development Status :: Alpha',
            'Intended Audience :: Alpha Testers',
            'Operating System :: Linux :: RHEL',
            'Programming Language :: Python',
            'Topic :: Gemini',
            'Topic :: Data Reduction',
            'Topic :: Astronomy',
            ],
       )