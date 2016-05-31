#!/usr/bin/env python

"""
Setup script for gemini_python.

The tools and modules in this package are developed by the Gemini Data Processing Software Group.

In this package:
    astrodata : AstroData
    astrodata_Gemini : Recipes, Primitives and astrodata configurations for Gemini
    gempy : Gemini toolkits
    package_sample : An example, template for a new astrodata_<something> package.
    recipe_system : The Recipe System
    #iqtool : IQ assessment.  Used internally by the SOS-DAs
    

Usage:
  python setup.py install --prefix=/astro/iraf/rhux-x86_64-glibc2.5/gempylocal
  python setup.py sdist
"""

import os.path
import os
import re
import glob
import sys

from distutils.core import setup

svndir = re.compile('.svn')
fitsfile = re.compile('.fits$')
dotpy = re.compile('.py$')

PACKAGENAME = 'gemini_python'

# PACKAGES and PACKAGE_DIRS
#    Note: KL not sure what to do about astrodata.adutils.reduceutils.pyjamaface
ASTRODATA_MODULES = ['astrodata',
                     'astrodata.eti',
                     'astrodata.interface',
                     'astrodata.utils',
                     'astrodata.utils.future']
#FITSSTORE_MODULES = ['fitsstore']
GEMPY_MODULES = ['gempy',
                 'gempy.adlibrary',
                 'gempy.gemini',
                 'gempy.gemini.eti',
                 'gempy.library']

RS_MODULES = ['recipe_system',
              'recipe_system.adcc',
              'recipe_system.adcc.servers',
              'recipe_system.cal_service',
              'recipe_system.reduction'
             ]
# IQTOOL NOT FOR PUBLIC RELEASE, YET
#IQTOOL_MODULES = ['iqtool',
#                  'iqtool.iq',
#                  'iqtool.gemplotlib']
ADLIB_PACKAGES = ['FITS','Gemini']   # This is be used to form 'astrodata_Gemini' and 'astrodata_FITS'
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
SUBMODULES.extend(ASTRODATA_MODULES)
SUBMODULES.extend(RS_MODULES)
SUBMODULES.extend(GEMPY_MODULES)
SUBMODULES.extend(RECIPE_MODULES)
SUBMODULES.extend(PIF_MODULES)
SUBMODULES.extend(ADCONFIG_MODULES)
#SUBMODULES.extend(IQTOOL_MODULES)


PACKAGES = []
PACKAGES.extend(SUBMODULES)
for p in ADLIB_PACKAGES:
    PACKAGES.append('astrodata_'+p)
PACKAGE_DIRS = {}
PACKAGE_DIRS[''] = '.'

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

    if os.path.isdir(os.path.join('astrodata_'+p,'RECIPES_'+p)):
        PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p].append('recipe.*')
        PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p].append(os.path.join('subrecipes','recipe.*'))
        PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p].append(os.path.join('demos','recipe.*'))
        PACKAGE_DATA['astrodata_'+p+'.RECIPES_'+p].append(os.path.join('tests','recipe.*'))

astrodatadir = re.compile('astrodata/')
for root, dirs, files in os.walk(os.path.join('astrodata','tests')):
    if not svndir.search(root) and len(files) > 0:
        PACKAGE_DATA['astrodata'].extend( map((lambda f: os.path.join(astrodatadir.sub('',root), f)), files) )

rsdir = re.compile('recipe_system/')
for root, dirs, files in os.walk(os.path.join('recipe_system','adcc','client','adcc_faceplate')):
    if not svndir.search(root) and len(files) > 0:
        PACKAGE_DATA['recipe_system'].extend( map((lambda f: os.path.join(rsdir.sub('',root), f)), files) )

gempydir = re.compile('gempy/')
for root, dirs, files in os.walk(os.path.join('gempy','tests')):
    if not svndir.search(root) and len(files) > 0:
        PACKAGE_DATA['gempy'].extend( map((lambda f: os.path.join(gempydir.sub('',root), f)), files) )

    
# DATA_DIRS and DATA_FILES
DATA_FILES = []

ASTRODATADOC_DIR = os.path.join('share','astrodata')
for root, dirs, files in os.walk(os.path.join('astrodata','doc')):
    if not svndir.search(root) and len(files) > 0:
        dest = root.split('/',2)[2] if len(root.split('/',2)) > 2 else ""
        DOC_FILES = map((lambda f: os.path.join(root,f)), files)
        DATA_FILES.append( (os.path.join(ASTRODATADOC_DIR,dest), DOC_FILES) )

for root, dirs, files in os.walk(os.path.join('package_sample')):
    if not svndir.search(root) and len(files) > 0:
        dest = root.split('/',1)[1] if len(root.split('/',1)) > 1 else ""
        DOC_FILES = map((lambda f: os.path.join(root,f)), files)
        DATA_FILES.append( (os.path.join(ASTRODATADOC_DIR,dest), DOC_FILES) )
        
# GEMINIDOC_DIR = os.path.join('share','astrodata_Gemini')
# for root, dirs, files in os.walk(os.path.join('astrodata_Gemini','doc')):
#     if not svndir.search(root) and len(files) > 0:
#         dest = root.split('/',2)[2] if len(root.split('/',2)) > 2 else ""
#         DOC_FILES = map((lambda f: os.path.join(root,f)), files)      
#         DATA_FILES.append( (os.path.join(GEMINIDOC_DIR,dest), DOC_FILES) )

#FITSDOC_DIR = os.path.join('share','astrodata_FITS')
#for root, dirs, files in os.walk(os.path.join('astrodata_FITS','doc')):
#    if not svndir.search(root) and len(files) > 0:
#        dest = root.split('/',2)[2] if len(root.split('/',2)) > 2 else ""
#        DOC_FILES = map((lambda f: os.path.join(root,f)), files)      
#        DATA_FILES.append( (os.path.join(FITSDOC_DIR,dest), DOC_FILES) )

GEMPYDOC_DIR = os.path.join('share','gempy')
for root, dirs, files in os.walk(os.path.join('gempy','doc')):
    if not svndir.search(root) and len(files) > 0:
        dest = root.split('/',2)[2] if len(root.split('/',2)) > 2 else ""
        DOC_FILES = map((lambda f: os.path.join(root,f)), files)      
        DATA_FILES.append( (os.path.join(GEMPYDOC_DIR,dest), DOC_FILES) )
for root, dirs, files in os.walk(os.path.join('gempy','doc-local')):
    if not svndir.search(root) and len(files) > 0:
        dest = root.split('/',2)[2] if len(root.split('/',2)) > 2 else ""
        DOC_FILES = map((lambda f: os.path.join(root,f)), files)      
        DATA_FILES.append( (os.path.join(GEMPYDOC_DIR,dest), DOC_FILES) )


# SCRIPTS
ASTRODATA_SCRIPTS = [ os.path.join('astrodata','scripts','mkCalculatorInterface'),
                     os.path.join('astrodata','scripts','showd'),
                     os.path.join('astrodata','scripts','superclean'),
                     os.path.join('astrodata','scripts','typewalk'),
                     os.path.join('astrodata','scripts','rsifaces','pif2prim','mkPIF'),
                    ]
RS_SCRIPTS = [ os.path.join('recipe_system','apps','adcc'),
               os.path.join('recipe_system','apps','listprimitives'),
               os.path.join('recipe_system','apps','reduce')
             ]

GEMPY_SCRIPTS = [ #os.path.join('gempy','scripts','cleanir.py')
                  os.path.join('gempy','scripts','fwhm_histogram'),
                  os.path.join('gempy','scripts','profile_all_obj'),
                  os.path.join('gempy','scripts','psf_plot'),
                  os.path.join('gempy','scripts','zp_histogram')
                 ]
#IQTOOL_SCRIPTS = [ os.path.join('iqtool','iqtool.py')]

if "sdist" in sys.argv:
    #GEMPY_SCRIPTS and ASTRODATA_SCRIPTS contain the name of the links which might not be dereferenced during sdist
    #Therefore, here we package the .py those links point to.  During "install" the links are
    #dereferenced, always, as far as I can tell, so there's no need for the .py then.
    PYFILES = []
    dotpy = re.compile(".py$")
    for script in GEMPY_SCRIPTS:
        if not dotpy.match(script):
            PYFILES.append(''.join([script,'.py']))
    GEMPY_SCRIPTS.extend(PYFILES)
    for script in ASTRODATA_SCRIPTS:
        if not dotpy.match(script):
            PYFILES.append(''.join([script,'.py']))
    ASTRODATA_SCRIPTS.extend(PYFILES)
    for script in RS_SCRIPTS:
        if not dotpy.match(script):
            PYFILES.append(''.join([script,'.py']))
    RS_SCRIPTS.extend(PYFILES)


SCRIPTS = []
SCRIPTS.extend(ASTRODATA_SCRIPTS)
SCRIPTS.extend(RS_SCRIPTS)
SCRIPTS.extend(GEMPY_SCRIPTS)
#SCRIPTS.extend(IQTOOL_SCRIPTS)

EXTENSIONS = None

setup ( name='gemini_python',
        version='0.2.0',
        description='Gemini Data Processing Python Package',
        author='Gemini Data Processing Software Group',
        author_email='klabrie@gemini.edu',
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
