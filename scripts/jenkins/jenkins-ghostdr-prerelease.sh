#!/bin/bash

# This script is designed to be called by the jenkins job 'ghostdr', when a release is requested.
# In particular it makes use of variables that are set by jenkins and won't usually be in your environment.
#

# Exit if any command returns an error
set -e

echo "This is the jenkins pre-release script"
echo "RELEASE_VERSION = ${RELEASE_VERSION}"
echo "RELEASE_REVISION = ${RELEASE_REVISION}"
echo "RELEASE_DESCRIPTION = ${RELEASE_DESCRIPTION}"

# Prepend the details of the new release to the file that contains the list of releases
RELEASE_STR="${RELEASE_VERSION}
  ${RELEASE_DESCRIPTION}

"
# Only mess with the repository from one host, not all the build hosts.  Might as well be the master.
if [ "${NODE_NAME}" == "master" ]; then
  # Start in the repository
  cd ${WORKSPACE}/hg

  # Update to the specified revision
  hg update ${RELEASE_REVISION}

  # Edit the releases file
  echo "3a
${RELEASE_STR}
.
w" | ed ${WORKSPACE}/hg/rtd/releases.rst >& /dev/null

  # Commit the updated versions of the file
  hg commit -m "Jenkins prerelease ${RELEASE_VERSION}" ${WORKSPACE}/hg/rtd/releases.rst

  # Tag the repository with this updated file
  hg tag ${RELEASE_VERSION}

  # Now push the tag and the update, but don't run the usual hook
  hg push --config hooks.incoming.jenkins=
fi
