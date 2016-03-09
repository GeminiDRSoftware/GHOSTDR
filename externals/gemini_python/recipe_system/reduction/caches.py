#
#                                                                  gemini_python
#
#                                                                      caches.py
# ------------------------------------------------------------------------------
# $Id: caches.py 5142 2015-02-17 21:39:45Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5142 $'[11:-2]
__version_date__ = '$Date: 2015-02-17 11:39:45 -1000 (Tue, 17 Feb 2015) $'[7:-3]
# ------------------------------------------------------------------------------
from os.path import join

# GLOBAL/CONSTANTS (could be exported to config file)
# [DEFAULT]
cals = "calibrations"

# [caches]
cachedirs = [".reducecache", cals, 
             join(cals,"storedcals"), 
             join(cals,"retrievedcals")]

CALDIR     =  join(cals,"storedcals")
adatadir   = "./recipedata/"
calindfile = "./.reducecache/calindex.pkl"
stkindfile = "./.reducecache/stkindex.pkl"

#".reducecache/storedcals/storedbiases",
#".reducecache/storedcals/storeddarks",
#".reducecache/storedcals/storedflats",
#".reducecache/storedcals/storedfringes",
#".reducecache/storedcals/retrievedbiases",
#".reducecache/storedcals/retrieveddarks",
#".reducecache/storedcals/retrievedflats",
#".reducecache/storedcals/retrievedfringes",
