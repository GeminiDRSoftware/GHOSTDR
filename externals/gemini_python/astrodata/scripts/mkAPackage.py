#!/usr/bin/env python
#
#                                                                  gemini_python
#
#                                                              astrodata/scripts
#                                                                  mkAPackage.py
# ------------------------------------------------------------------------------
# $Id: mkAPackage.py 5178 2015-03-26 20:42:38Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5178 $'[11:-2]
__version_date__ = '$Date: 2015-03-26 10:42:38 -1000 (Thu, 26 Mar 2015) $'[7:-2]
# ------------------------------------------------------------------------------
import os
import sys
import shutil

from argparse import ArgumentParser
# ------------------------------------------------------------------------------
tempPack = "astrodata_Sample"
tempPackTag = "_Sample"

parser = ArgumentParser()
parser.add_argument("-w", dest="whereout", default=".", 
                    help="Path to new package. Default is '.'")
parser.add_argument('package_name', help="Name of new package")
args = parser.parse_args()

outpackTag = "_"+args.package_name
outdir = os.path.abspath(args.whereout)

dirname = os.path.dirname(__file__)
dirname = os.path.abspath(os.path.join(dirname, "../../package_sample/"))
templatePack = os.path.join(dirname, tempPack)
outPack = os.path.abspath(os.path.join(outdir, "astrodata"+outpackTag))

print ("cloning\n\t%s \nto\n\t%s\n" %(templatePack, outPack))
shutil.copytree(templatePack, outPack)
print ("...copied, renaming subdirectories")

movelist = []
print "outpack",outPack

for root, dirs, files in os.walk(outPack):
    if (".svn" in root):
        shutil.rmtree(root)

for root, dirs, files in os.walk(outPack):
    print "in %s" % root,
    gotone = False
    for fil in dirs:
        if tempPackTag in fil:
            if gotone == False:
                print ""
            gotone = True
            newfil = fil.replace(tempPackTag, outpackTag)
            print "\tplanning to move %s -> %s" % ( fil, newfil )
            fullfil = os.path.join(root, fil)
            fullnewfil = os.path.join(root, newfil)
            movelist.insert(0,(fullfil,fullnewfil))

    for fil in files:
        if tempPackTag in fil:
            if gotone == False:
                print ""
            gotone = True
            newfil = fil.replace(tempPackTag, outpackTag)
            print "\tplanning to move %s -> %s" % ( fil, newfil )
            fullfil = os.path.join(root, fil)
            fullnewfil = os.path.join(root, newfil)
            movelist.insert(0,(fullfil,fullnewfil))              

    if not gotone:
        print "... nothing to change"

infostr = "Renaming %d Directories" % len(movelist)
print "="*len(infostr)
print infostr
print "-"*len(infostr)

for mov in movelist:
        print "Moving %s" % mov[0]
        print "   --> %s" % mov[1]
        shutil.move(mov[0], mov[1])
