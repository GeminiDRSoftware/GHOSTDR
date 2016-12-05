#!/usr/bin/env python

"""
Setup script for ghostdr

Based on the setup.py from the gemini_python distribution.

In this package:
    astrodata_GHOST  : Recipes, Primitives and astrodata configuration for GHOST
    pyghost          : GHOST simulator
    polyfit          : GHOST find apertures and polynomial model fitting
    

Usage:
  python setup.py install
"""

import os.path
import os
import re
import glob
import sys

from setuptools import setup

svndir = re.compile('.hg')
fitsfile = re.compile('.fits$')
dotpy = re.compile('.py$')

PACKAGENAME = 'ghostdr'

ADLIB_PACKAGES = ['GHOST']   # This is be used to form 'astrodata_GHOST'
RECIPE_MODULES=[]
PIF_MODULES=[]
ADCONFIG_MODULES=[]
slash = re.compile('/')
for p in ADLIB_PACKAGES:
    if os.path.isdir(os.path.join('astrodata_'+p,'RECIPES_'+p)): 
        RECIPE_MODULES.append('astrodata_'+p+'.RECIPES_'+p)
        RECIPE_MODULES.append('astrodata_'+p+'.RECIPES_'+p+'.primitives')
    if os.path.isdir(os.path.join('astrodata_'+p,'PIF_'+p)):
        PIF_MODULES.append('astrodata_'+p+'.PIF_'+p)
        PIF_MODULES.append('astrodata_'+p+'.PIF_'+p+'.pif'+p.lower())
        PIFROOT = os.path.join('astrodata_'+p,'PIF_'+p,'pif'+p.lower())
        for root, dirs, files in os.walk(PIFROOT):
            if not svndir.search(root) and len(files) > 0:
                pifmodules = map((lambda d: slash.sub('.','/'.join([root,d]))),\
                                 filter((lambda d: not svndir.search(d)), dirs))
                PIF_MODULES.extend( pifmodules )
    if os.path.isdir(os.path.join('astrodata_'+p, 'ADCONFIG_'+p)):
        ADCONFIG_MODULES.append('astrodata_'+p+'.ADCONFIG_'+p)
        ADCONFIG_MODULES.append('astrodata_'+p+'.ADCONFIG_'+p+'.descriptors')
        if os.path.isdir(os.path.join('astrodata_'+p, 'ADCONFIG_'+p, 'lookups')):
            ADCONFIG_MODULES.append('astrodata_'+p+'.ADCONFIG_'+p+'.lookups')
            LUTROOT = os.path.join('astrodata_'+p,'ADCONFIG_'+p,'lookups')
            print LUTROOT
            for root, dirs, files in os.walk(LUTROOT):
                print root, dirs, files
                if not svndir.search(root) and len(files) > 0:
                    lutmodules = map((lambda d: slash.sub('.','/'.join([root,d]))),\
                                     filter((lambda d: not svndir.search(d)), dirs))
                    ADCONFIG_MODULES.extend( lutmodules )

SUBMODULES = []
SUBMODULES.extend(ADCONFIG_MODULES)

PACKAGES = []
PACKAGES.extend(SUBMODULES)
for p in ADLIB_PACKAGES:
    PACKAGES.append('astrodata_'+p)

PACKAGE_DIRS = {}
#PACKAGE_DIRS[''] = '.'

# PACKAGE_DATA
PACKAGE_DATA = {}

for s in SUBMODULES:
    PACKAGE_DATA[s] = []
#    PACKAGE_DATA[s] = ['ReleaseNote',
#                       'README',
#                       'INSTALL'
#                       ]
for p in ADLIB_PACKAGES:
    #PACKAGE_DATA['astrodata_'+p] = \
    #                  ['ReleaseNote',
    #                   'README',
    #                   'INSTALL',
    #                   os.path.join('ADCONFIG_'+p,'structures','*.py'),
    #                   ]
    PACKAGE_DATA['astrodata_'+p] = []
    for root, dirs, files in os.walk(os.path.join('astrodata_'+p,'ADCONFIG_'+p,'lookups')):
        if not svndir.search(root) and len(files) > 0:
            files = [f for f in files if fitsfile.search(f)]
            dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
            PACKAGE_DATA['astrodata_'+p].extend( map((lambda f: os.path.join(dest, f)), files) )
    for root, dirs, files in os.walk(os.path.join('astrodata_'+p,'ADCONFIG_'+p,'lookups','source_detection')):
        if not svndir.search(root) and len(files) > 0:
            files = [f for f in files if not dotpy.search(f)]
            dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
            PACKAGE_DATA['astrodata_'+p].extend( map((lambda f: os.path.join(dest, f)), files) )
    for root, dirs, files in os.walk(os.path.join('astrodata_'+p,'ADCONFIG_'+p,'descriptors')):
        if not svndir.search(root) and len(files) > 0:
            dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
            PACKAGE_DATA['astrodata_'+p].extend( map((lambda f: os.path.join(dest, f)), files) )
    for root, dirs, files in os.walk(os.path.join('astrodata_'+p,'ADCONFIG_'+p,'classifications')):
        if not svndir.search(root) and len(files) > 0:
            dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
            PACKAGE_DATA['astrodata_'+p].extend( map((lambda f: os.path.join(dest, f)), files) )

    PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p] = []
    if os.path.isdir(os.path.join('astrodata_'+p,'RECIPES_'+p)):
        PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p].append('recipe.*')
        PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p].append(os.path.join('subrecipes','recipe.*'))
        PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p].append(os.path.join('demos','recipe.*'))
        PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p].append(os.path.join('tests','recipe.*'))

# DATA_DIRS and DATA_FILES
DATA_FILES = []

SCRIPTS = []
ENTRY_POINTS = {}

# Add the simulator
PACKAGES.append('pyghost')
PACKAGE_DIRS['pyghost'] = 'simulator/pyghost/pyghost'
ENTRY_POINTS['console_scripts'] = ['ghost-sim=pyghost.tests:run']
PACKAGE_DATA['pyghost'] = ['data/*']


EXTENSIONS = None

setup ( name='ghostdr',
        version='0.1.0',
        description='GHOST Data Reduction',
        author='RSAA Software Engineering',
        author_email='ghostdr@mso.anu.edu.au',
        url='http://rsaa.anu.edu.au',
        maintainer='Gemini Data Processing Software Group',
        packages=PACKAGES,
        package_dir=PACKAGE_DIRS,
        package_data=PACKAGE_DATA,
        data_files=DATA_FILES,
        scripts=SCRIPTS,
	entry_points=ENTRY_POINTS,
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
