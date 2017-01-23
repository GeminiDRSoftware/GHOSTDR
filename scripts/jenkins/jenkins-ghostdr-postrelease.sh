#!/bin/bash

# This script is designed to be called by the jenkins job 'ghost_cicada', when a release is requested.
# In particular it makes use of variables that are set by jenkins and won't usually be in your environment.
#

# Exit if any command returns an error
set -e

# Create the message
GITMSG="GHOST data reduction software release ${RELEASE_VERSION}"

# make the new web pages from the latest rst files
RTD=${WORKSPACE}/hg/rtd
PROGMAN=${WORKSPACE}/hg/astrodata_GHOST/docs/progmanuals/GHOST_ProgManual
USERMAN=${WORKSPACE}/hg/astrodata_GHOST/docs/usermanuals/GHOST_UsersManual
DEST=/priv/www/mssso/ghostdr

make -C ${RTD} SPHINXBUILD=/usr/local/python-2.7.1/bin/sphinx-build html
make -C ${PROGMAN} SPHINXBUILD=/usr/local/python-2.7.1/bin/sphinx-build html
make -C ${USERMAN} SPHINXBUILD=/usr/local/python-2.7.1/bin/sphinx-build html

# remove old software release web pages and install new ones
rm -rf ${DEST}/*
cp -r ${RTD}/_build/html/* ${DEST}
cp -r ${RTD}/web/* ${DEST}
mkdir -p ${DEST}/progman
cp -r ${PROGMAN}/_build/html/* ${DEST}/progman
mkdir -p ${DEST}/userman
cp -r ${USERMAN}/_build/html/* ${DEST}/userman

# copy everything into the GHOSTDR_github directory
rsync -a --exclude "externals" --exclude ".hg*" ${WORKSPACE}/hg/ ${JENKINS_HOME}/workspace/GHOSTDR_github/

# Push the changes back into the github repository
cd ${JENKINS_HOME}/workspace/GHOSTDR_github
export HOME=${JENKINS_HOME}
git add .
git commit -m "${GITMSG}"
git tag -a -m "${GITMSG}" ${RELEASE_VERSION}
git push origin ${RELEASE_VERSION}

# send out a new notification
mail -s "GHOST data reduction software package release ${RELEASE_VERSION}" ghostdr-release@mso.anu.edu.au <<< "GHOST data reduction software package ${RELEASE_VERSION} has been released.  You can find more information at http://www.mso.anu.edu.au/ghostdr/"
