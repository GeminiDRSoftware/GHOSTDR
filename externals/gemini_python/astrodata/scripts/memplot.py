#!/usr/bin/env python
import matplotlib
from pylab import *
figure (figsize=(20,10))

import sys,os
import json

print len(sys.argv), sys.argv

def filize(name):
    charlist = []
    for x in name:
        charlist.append( x if x.isalnum() else "_")    
    clean ="".join(charlist)
    return clean
    
fname = sys.argv[1]

outpng = "%s.png" % fname
label = fname

if (len(sys.argv)>=3):
    label = sys.argv[2]
    outpng = "%s.png" % filize ( label )
    

mifile = open(fname)

mifary = mifile.readlines()

# make a memtrack file with better name
miname = "%s.mem" % filize(label)

num = 1
while os.path.exists(miname):
    miname = "%s.%d.mem" % (filize(label), num)
    num += 1

micpy = open(miname, "w")
micpy.writelines(mifary)
micpy.close()



mifile.close()

rssvals = []
stepnames = []
diffitems = []
titems = []
lastitem = None;
for line in mifary:
    item = json.loads(line)
    item["rss"] = float(item["rss"])
    titems.append(item)
    item["index"] = titems.index(item)
    
    rss = item["rss"]/1e6
    rssvals.append(rss)
    stepnames.append(item["msg"])
    if lastitem:
        print "#%d current = %f, last = %f" % (item["index"], item["rss"],lastitem["rss"])
        diff = item["rss"]-lastitem["rss"]
        print diff, diff>20e6
        if abs(diff)>5e6:
            print type(diff),type(1e6)
            diffstr =  "%2.2f MB" % (float(diff)/1000000.0)
            diffitems.append({
                            "item0": lastitem,
                            "item1": item,
                            "diff": diff,
                            "diffstr": diffstr
                            })
    lastitem = item
                
from numpy import arange

#yvals = frange(0,len(rssvals))
yvals = arange(0.,len(rssvals))

grid("on")

plot(yvals, rssvals)
ymax = max(rssvals)
ymin = min(rssvals)
ymid = ymin + (ymax-ymin)/2

align = None
yspot = ymax
#if False: 

# diff background
for diff in diffitems:
    i0 = diff["item0"]   
    i1 = diff["item1"]
    aty = i0["rss"]   + ((i1["rss"]-i0["rss"]) / 2.0)
    aty = aty/1e6
    atx = i0["index"] + ((i1["index"]-i0["index"]) / 2.0)
    print "process diffitem: ", atx, aty, diff["diffstr"]
    
   
    if diff["diff"] > 0:
        axvspan( xmin = i0["index"], xmax = i1["index"],
            color ="MistyRose"
            )
    else:
        axvspan( xmin = i0["index"], xmax = i1["index"],
            color ="#e4ffe1"
            )
                    
for val in yvals:
    val = int(val)
    color = "blue"
    item = titems[val]
    rss = item["rss"]
    rssmb = rss/1e6
    if rssmb > ymid:
        yspot = 25.0 #rssmb - 25.
        align = "bottom"
    else:
        yspot = 25.0 #rssmb + 25.
        align = "bottom"
    if "END" in stepnames[val]:
        color = "red"
    if "START" in stepnames[val]:
        color = "green"
    print "min, mid, max", ymin, ymid, ymax
    print "primitive step label:", val, yspot, align
    text(val,yspot, stepnames[val], 
            horizontalalignment="center",
            verticalalignment=align,
            rotation="vertical",
            backgroundcolor = (1,1,1,.5),
            color = color)
 
# diff labels diffstr is prably XX MB
for diff in diffitems:
    i0 = diff["item0"]   
    i1 = diff["item1"]
    aty = i0["rss"]   + ((i1["rss"]-i0["rss"]) / 2.0)
    aty = aty/1e6
    atx = float(i0["index"]) + float((i1["index"]-i0["index"]) / 2.0)
    print "process diffitem: ", atx, aty, diff["diffstr"]
    
    #text( 10, 200.0e6, "label")
    text(atx, aty, diff["diffstr"],
        color = "black",
        rotation = "vertical",
        horizontalalignment = "center",
        verticalalignment = "center"
        )
          
title ( label)
xlabel ( "recipe step #" )
ylabel ( "real memory in Millions of Bytes")


savefig(outpng, dpi=42)
show()