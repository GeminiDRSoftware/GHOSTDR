#!/bin/bash
#
# This script downloads Ureka from ${UREKA_URL} and installs it into
# ${UREKA_DIR}.
#

# These are customizable, if required
GEMINI_DIR="${HOME}/gemini"
UREKA_DIR="${GEMINI_DIR}/Ureka"
UREKA_URL="http://ssb.stsci.edu/ureka/1.5.2/install_ureka_1.5.2"

# Don't edit anything below here
UREKA_SCRIPT=`basename $UREKA_URL`

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

# Remove any old installer script
rm -f ${UREKA_SCRIPT}
# Get the new installer script
my_get ${UREKA_URL}
# Remember where we were (so we can run the script later)
SCRIPT_DIR=`pwd`
# Make the destination directory if required
mkdir -p ${GEMINI_DIR}
# Clean up any old Ureka installation
rm -rf ${UREKA_DIR}
# Go to the install directory (required by the Ureka install script)
cd ${GEMINI_DIR}
# And run the Ureka install script
bash ${SCRIPT_DIR}/${UREKA_SCRIPT}
