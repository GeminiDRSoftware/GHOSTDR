#!/bin/bash
#
# This script downloads test data from ${DATA_URL} and installs it into
# ${GEMINI_DIR}.
#

# These are customizable, if required
GEMINI_DIR="${HOME}/gemini"
DATA_DIR="${GEMINI_DIR}/data_for_ad_unit_tests"
DATA_URL="http://www.gemini.edu/sciops/data/software/gemini_python/gemini_python_datapkg-X1.tar.gz"

# Don't edit anything below here
DATA_FILE=`basename $DATA_URL`

#
# A function to get the given url, using wget or curl whichever is
# available.
#
my_get() {
    if command -v wget >/dev/null 2>&1; then
        wget "$@"
    elif command -v curl >/dev/null 2>&1; then
        curl -O "$@"
    else
        echo "Could not find curl or wget. Aborting."
	exit 1;
    fi
}

# Go to the gemini directory
cd ${GEMINI_DIR}
rm -rf ${DATA_DIR}
# Get the data
my_get ${DATA_URL}
# Extract it
gunzip -c ${DATA_FILE} | tar xf -
