#!/bin/bash

# This script is designed to be called by the jenkins job 'ghostdr'.
# In particular it makes use of variables that are set by jenkins and won't usually be in your environment.
#

# Exit if any command returns an error
set -e

# Change to the appropriate directory
cd ${WORKSPACE}/hg

# Run pylint over our codebase
/usr/local/python-2.7.1/bin/pylint simulator/pyghost/pyghost || exit 0
