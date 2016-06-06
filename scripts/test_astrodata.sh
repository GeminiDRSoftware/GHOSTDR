#/bin/bash
#
# This script runs the astrodata tests.
#
export ADTESTDATA=$HOME/gemini/gemini_python_datapkg-X1/data_for_ad_unit_tests
cd ../externals/gemini_python/astrodata/tests
py.test -rxX
