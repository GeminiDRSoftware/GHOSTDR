#!/bin/bash

# This script is designed to be called by the jenkins job 'ghostdr'.
# In particular it makes use of variables that are set by jenkins and won't usually be in your environment.
#

# Exit if any command returns an error
set -e

# If this is a release version then use a real version number, otherwise make up a fake one to identify the build.
cd ${WORKSPACE}/hg

# Run pylint over our codebase
/usr/local/python-2.7.1/bin/pylint simulator/pyghost/pyghost

# Now we build everything
# make
