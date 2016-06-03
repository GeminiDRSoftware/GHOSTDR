#!/bin/bash

# This script is designed to be called by the jenkins job 'ghost_cicada', when a release is requested.
# In particular it makes use of variables that are set by jenkins and won't usually be in your environment.
#

# Exit if any command returns an error
set -e

# Create the message
GITMSG="GHOST data reduction software release ${RELEASE_VERSION}"

# make the new web pages from the latest rst files
make -C ${WORKSPACE}/hg/rtd SPHINXBUILD=/usr/local/python-2.7.1/bin/sphinx-build html

# remove old software release web pages and install new ones
rm -rf /priv/www/mssso/ghostdr/*
cp -r ${WORKSPACE}/hg/rtd/_build/html/* /priv/www/mssso/ghostdr

# copy everything into the GHOSTDR_github directory
rsync -a --exclude "externals" --exclude ".hg*" ${WORKSPACE}/hg/ ${JENKINS_HOME}/workspace/GHOSTDR_github/

# Push the changes back into the github repository
cd ${JENKINS_HOME}/workspace/GHOSTDR_github
git add .
git commit -m "${GITMSG}"
git tag -a -m "${GITMSG}" ${RELEASE_VERSION}

# We rely on the GHOSTDR_github project to do the actual push, since it has
# the relevant github credentials in it.  And that project is triggered
# automatically when this job completes.

# send out a new notification
mail -s "GHOST data reduction software package release ${RELEASE_VERSION}" ghostdr-release@mso.anu.edu.au <<< "GHOST data reduction software package ${RELEASE_VERSION} has been released.  You can find more information at http://www.mso.anu.edu.au/ghostdr/"
