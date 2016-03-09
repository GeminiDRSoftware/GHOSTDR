#!/usr/bin/env python
import os
import json
import fcntl
import socket
import struct
import urllib
import threading
from time import time
from pprint import pformat
from optparse import OptionParser
from urlparse import urlparse, parse_qs

import astrodata
from astrodata.utils import jsutil
from astrodata.interface import AstroDataType
from astrodata.utils.ConfigSpace import cs
from astrodata.interface.AstroDataType import ClassificationLibrary

from recipe_system.reduction.priminspect import PrimInspect
from recipe_system.reduction..RecipeManager import RecipeLibrary

recipelibrary = RecipeLibrary()
classificationlib = ClassificationLibrary()

# ------------------------------------------------------------------------------
def get_ip_address(ifname):
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    return socket.inet_ntoa(fcntl.ioctl(
        s.fileno(),
        0x8915,  # SIOCGIFADDR
        struct.pack('256s', ifname[:15])
    )[20:24])

# ------------------------------------------------------------------------------
# set up commandline args and options
parser = OptionParser()
parser.set_description( "Gemini Observatory Pipeline Development Kit Tool"
                       "(v_1.0 2011)")
parser.add_option("-d", "--server", action="store_true", dest="server",
                  default=True,
                  help="start a server for interactive use via browser")

parser.add_option("-p", "--port", dest = "port", default = 8888, type = "int",
                  help="ADKTOOL service provided on <port> (http://localhost:<port>)")
       
parser.add_option("--nonlocal", action="store_true", dest="nonlocal",
                help="used to put the server in non-local mode (no editing)")
                
(options,  args) = parser.parse_args()
server = options.server

# ------------------------------------------------------------------------------
class EditThread(threading.Thread):
    def __init__(self, fname, lnum):
        self.fname = fname
        self.lnum = lnum
        threading.Thread.__init__(self)

    def run(self):
        from subprocess import call
        if "PDK_EDITOR" in os.environ:
           PDK_EDITOR = os.environ["PDK_EDITOR"]
        else:
            PDK_EDITOR = "nedit +%(lnum)s %(filename)s"
        editcmd = PDK_EDITOR % {
                                "lnum":self.lnum,
                                "filename":self.fname
                               }
        
        print "pt41:", editcmd
        call ( editcmd, shell=True)


class RCat():
    def __init__(self, trd):
        self.cat = trd
                    
    def makeCats(self,catlist, curcat = None):
        if curcat == None:
            curcat = self.cat
        if len(catlist) ==0:
            return curcat
        if "subcats" not in curcat:
            curcat.update({"subcats":[]})
        subcats = curcat["subcats"]
        goodcat = None
        for cat in subcats:
            if cat["cat_name"] == catlist[0]:
                goodcat = cat
        if goodcat == None:
            goodcat = {"cat_name":catlist[0]}
            subcats.append(goodcat)
        subcat = goodcat
        return self.makeCats(catlist[1:], curcat = subcat)

    def add(self, nd):
        path = nd["path"]
        category = os.path.dirname(path)
        recipename = os.path.basename(path)
        cats = category.split("/")
        nd.update({"cat_name":cats[-1]})
        if nd["package"] != "*":
            cats = cats[cats.index(nd["package"])+1:]
        else: # package does equal "*", the wildcard, set when typeinfo == True
            for i in range(0,len(cats)):
                if "astrodata_" in cats[i]:
                    break
            cats =cats[i+1:]
            correctcat = self.makeCats(cats)
            if "members" not in correctcat:
                correctcat.update({"members":[]})
            members = correctcat["members"]
            members.append(nd)

    def printDict(self,pdict = None, indent = ""):
        if pdict == None:
            pdict = self.cat
        print indent,"p-------t"
                    
        if "name" in pdict:
            print indent,"name",pdict["name"]
        if "cat_name" in pdict:
            print indent,"cat_name", pdict["cat_name"]
        if "members" in pdict:
            membs = pdict["members"]
            memnames = [ mem["name"] for mem in membs]
            print indent, ", ".join(memnames)
        if "subcats" in pdict:
            for cat in pdict["subcats"]:
                self.printDict(cat, indent=indent+" ")
                    
if server:
    import BaseHTTPServer
    
    HOST_NAME = ""
    PORT_NUMBER = options.port
    
    class IfaceHandler(BaseHTTPServer.BaseHTTPRequestHandler):
                
        pi = None;     
        
        def address_string(self):
            host, port = self.client_address[:2]
            return host   
        def do_HEAD(s):
            s.send_response(200)
            s.send_header("Content-type", "text/html")
            s.end_headers()
       
        def do_GET(s):
            """Respond to a GET request."""
            try:
                myIP = get_ip_address('eth0')
            except:
                print "Can't Get Local Host Addr on 'etho0'."
            myIP = None
            self = s
            clientIP = self.client_address[0]
            print "at159: host", myIP,"client", clientIP
            if clientIP == myIP or clientIP == "127.0.0.1":
                localClient = True
            else:
                localClient = False            
            
            print "at165:",localClient,options.nonlocal
            if localClient == True and options.nonlocal == True:
                localClient = False
                              
            self = s
            # If someone went to "http://something.somewhere.net/foo/bar/",
            # then s.path equals "/foo/bar/".
            #s.wfile.write("<p>You accessed path: %s</p>" % s.path)
            pres = urlparse(s.path)
            urlparms = parse_qs(pres.query)
            
            cmds = pres.path.split("/")[1:]
            command = scope = subscope = None
            if len(cmds)>0:
                command = cmds[0]
            if len(cmds)>1:
                scope  = cmds[1]
            if len(cmds)>2:
                subscope = cmds[2]
            pparams = pres.query
            print "186:", command, scope, subscope, "urlparms:", repr(urlparms)           
            
            if command=="qap" or command=="docs":
                if ".." in self.path:
                    self.send_response(200)
                    self.send_header('Content-type', 'text/html')
                    self.end_headers()
                    data = "<b>bad path error</b>"
                    self.wfile.write(data)
                dirname = os.path.dirname(astrodata.__file__)
                print "lP114:", dirname
                if command == "qap":
                    joinlist = [dirname, "scripts/adcc_faceplate/"]
                elif command == "docs":
                    joinlist = [dirname, "doc/docscripts/build/"]
                    
                # Split out any parameters in the URL
                self.path = self.path.split("?")[0]
                print "adk204:", repr(self.path)
                #append any further directory info.
                joinlist.append( self.path[len(command)+2:])
                
                # print "ppw790:", repr(joinlist), self.path
                
                fname = os.path.join(*joinlist)
                # print "pt97: QAP IF: trying to open\n\t%s" % fname
                responsecode = 200
                try:
                    f = open(fname, "r")
                    data = f.read()
                    f.close()
                except IOError:
                    print "lP131: no such resource as",fname
                    data = "<b>NO SUCH RESOURCE AVAILABLE</b>"
                    responsecode = 404
                #print "pt100:",data
                self.send_response(responsecode)
                if  self.path.endswith(".js"):
                    self.send_header('Content-type', 'text/javascript')
                    self.send_header("Expires", "157785000")
                elif self.path.endswith(".css"):
                    self.send_header("Content-type", "text/css")
                elif fname.endswith(".png"):
                    self.send_header('Content-type', "image/png")
                else:
                    self.send_header('Content-type', 'text/html')
                self.end_headers()
                self.wfile.write(data)
                return 
                
            print "233: at second if:(%s)" % command
            # preprocess command for typedict arg    
            graphtype = "NICI" #"GMOS_MOS"      
            showprims = None
            if "showprims" in urlparms:
                from recipe_system.reduction.RecipeManager import centralPrimitivesIndex
                showprims = centralPrimitivesIndex
                
                
            if command == "typedict.py":            
                command = ""
                if "astrotype" in urlparms:
                    graphtype = urlparms["astrotype"][0]
            
            if command == "" or command.isdigit():
                s.send_response(200)
                s.send_header("Content-type", "text/html")
                s.end_headers()
            
                dirname = os.path.dirname(__file__)
                templname = os.path.join(dirname, "adcc_faceplate/adktool%s.html"% str(command))
                try:
                    template = open(templname)
                except IOError:
                    s.wfile.write("CAN'T FIND TEMPLATE: %s" % templname)
                    return
                template = template.read()
                template_args = {}
                
                selopts = []
                for pkpath in cs.package_paths:
                    packname = os.path.basename(pkpath)
                    path = pkpath
                    link = "/pkinfo/"+packname
                    selopts.append( { "name":packname,
                                      "path":path,
                                      "url":link
                                     })
                
                accordfeature = jsutil.JSAccord()
                template_args.update({"accordian_div":accordfeature.div()})

                sel = jsutil.JSPackageSelection(selopts)
                template_args.update({"select_options": sel.div()})
                
                typdiv = jsutil.JSTypes()
                template_args.update({"types_div": typdiv.div()})

                desdiv = jsutil.JSDescriptors()
                template_args.update({"descriptors_div": desdiv.div()})
     
                recdiv = jsutil.JSRecipeSystem()
                template_args.update({"recipes_div": recdiv.div()})

                from astrodata.interface.AstroDataTypes import get_classification_library
                classlib = get_classification_library()
                typeobj = classlib.get_type_obj("GEMINI")
                print "akdtool281:",pformat(typeobj.json_typetree())
                typestree_json = typeobj.json_typetree()
                
                
                # dot = classlib.gviz_doc(astrotype=graphtype, assign_dict = showprims)
                
                template_args.update(  
                            {
                                "local_client": "true" if localClient else "false",
                                # "types_tree_dot": dot,
                                "typestree_json":typestree_json,
                                
                            })              
                pagetext = template % template_args                          

#                s.wfile.write('<div id="pk_content" style="padding:2px;margin:2px;float:left;border:solid black 1px"></div>')

                s.wfile.write(pagetext)
                             
            elif command == "old": # main interface
                s.send_response(200)
                s.send_header("Content-type", "text/html")
                s.end_headers()
                s.wfile.write("""
<html><head><title>Astrodata Package Viewer</title>
<link href="/qap/js/jquery-ui/css/ui-lightness/jquery-ui-1.8.20.custom.css" rel=
"stylesheet">
<script type="text/javascript" src="/qap/js/jquery-ui/js/jquery.js"></script>
<script type="text/javascript" src="/qap/js/jquery-ui/js/jquery-ui.js"></script>

<script type="text/javascript">
function ajaxLink(link)
{
    $.ajax({url:link,
            type:"post",
            });
}


function SERVEReditFile(link,lnum)
{
    if (typeof(lnum) == "undefined")
    {
        lnum = 0;
    }
    $.ajax({url:"/edit",
            type:"post",
            data: { lnum:lnum,
                    filename: link
                   }
           });
    
}

function aceditFile(link,lnum)
{
    if (typeof(lnum) == "undefined")
    {
        lnum = 0;
    }
    
    var ddict = {   lnum:lnum,
                    filename:link
                };
    var normal = true;
    if (normal)
    {
        var acewin = window.open("/acedit?lnum="+lnum+"&filename="+encodeURI(link),
                                "_blank",
                                "width=600,height=600");
        // acewin.moveTo(window.screenX, window.screenY);
        acewin.moveTo(200,200);
        //acewin.resizeTo(600,500);
    }
    
    else
    {
        $.ajax({url:"/acedit",
                type:"get",
                data: { lnum:lnum,
                        filename: link
                       },
                success:function(data)
                {
                    var acewin = window.open("","","width=600,height=500");
                    // acewin.moveTo(window.screenX, window.screenY);
                    acewin.moveTo(0,0);
                    ac
                    acewin.document.write(data);
                    acewin.focus();
                }
               });
     }   
}

editFile = aceditFile;

function empty_pdk_focus(){
    $($(".pdk_focus")[0]).slideUp(500, function () {
                     $($(".pdk_focus")[0]).empty();
                     });
}
                                        
function getKeys(obj)
{
    var keys = []
    for (var key in obj)
    {
        keys[keys.length]=key;
    }
    keys.sort();
    return keys;
}
localClient = %(local_client)s;
</script>                    
                    </head>
                    """ % {"local_client": "true" if localClient else "false" }
                    )
                s.wfile.write("<body>")
                s.wfile.write("<a href='/docs'>Astrodata Manual</a><br/>")
                ###                
               
                selopts = []
                for pkpath in cs.package_paths:
                    packname = os.path.basename(pkpath)
                    path = pkpath
                    link = "/pkinfo/"+packname
                    selopts.append( { "name":packname,
                                      "path":path,
                                      "url":link
                                     })
                accordfeature = jsutil.JSAccord()
                s.wfile.write(accordfeature.div())
                sel = jsutil.JSPackageSelection(selopts)
                s.wfile.write(sel.div())
                
                #act = jsutil.JSAce()
                #s.wfile.write(act.div())
                
                s.wfile.write("""<div class="pdk_focus">
                                </div>
                                """)                
                
                typdiv = jsutil.JSTypes()
                s.wfile.write(typdiv.div())               

                desdiv = jsutil.JSDescriptors()
                s.wfile.write(desdiv.div())               
     
                recdiv = jsutil.JSRecipeSystem()
                s.wfile.write(recdiv.div())
                #typdiv = jsutil.JSTypes()
                #s.wfile.write(typdiv.div())               
                
#                s.wfile.write('<div id="pk_content" style="padding:2px;margin:2px;float:left;border:solid black 1px"></div>')
                s.wfile.write("</body></html>")
                            
            elif command == "acedit":
                st = time()
                print "got ace edit", st
                
                s.send_response(200)
                print "c1", time() -st, " elapsed"           
                s.send_header("Content-type", "text/html")
                print "c1", time() -st, " elapsed"           
                s.end_headers()
                
                print "before open file ace edit", urlparms["filename"], time() -st, " elapsed"           
                
                sendfile = open(urlparms["filename"][0])
                print "before read file ace edit", time() -st, " elapsed"           
                code = sendfile.read()
                print "read file ace edit", time() -st, " elapsed"           
                sendfile.close()
                code = code.replace("<","&lt;")
                code = code.replace(">","&rt;")
                print "pt353:",localClient
                if "target" in urlparms:
                    target = urlparms["target"][0]
                else:
                    target = "unknown"
                ace = jsutil.JSAce(code, 
                        int(urlparms["lnum"][0]), 
                        local_client = localClient,
                        filename = urlparms["filename"][0],
                        target = target
                        );
                s.wfile.write(ace.page())
                # print "pt303:", ace.page()     
                print "returned ace edit", time() -st, " elapsed"           
                return
                
                
                
            elif command == "pkinfo" or command == "typeinfo":
                packinfo   = False
                typeinfo = False
                astrotype = None    
                packname = None           
                if command == "pkinfo":
                    packinfo = True
                    packname = scope                
                if command == "typeinfo":
                    typeinfo = True
                    if "astrotype" in urlparms:
                        astrotype = urlparms["astrotype"]
                    else:
                        astrotype = scope
                
                cl = AstroDataType.ClassificationLibrary.get_classification_library()
                s.send_response(200)
                s.send_header("Content-type", "application/json")
                s.end_headers()

                rd = {"package":scope, "types":[]}
                retprimsets = rettypes = retdescriptors = retrecipes = retrecipesystem = False
                if packname:
                    if (subscope == None or subscope == "all"):
                        print "sub", subscope
                        
                        rettypes = True
                        retdescriptors = True
                        retrecipesystem = True
                        retprimsets = True
                        
                    if subscope == "types":
                        rettypes = True
                    if subscope == "descriptors":
                        retdescriptors = True
                    if subscope == "recipesystem":
                        retrecipesystem = True
                    if subscope == "primsets":
                        retrecipesystem = True
                if astrotype:
                    if (subscope == None or subscope == "all"):
                        print "sub", subscope
                        
                        rettypes = True
                        retdescriptors =    True
                        retrecipesystem =   True #True
                        retprimsets =       True #True
                        
                    if subscope == "types":
                        rettypes = True
                    if subscope == "descriptors":
                        retdescriptors = True
                    if subscope == "recipesystem":
                        retrecipesystem = True
                    if subscope == "primsets":
                        retrecipesystem = True
                
                # returning types    
                if rettypes:    
                    cns = []
                    if "type_meta" not in rd:
                        rd.update({"type_meta": {}})

                    if typeinfo:
                        typenames = [astrotype]
                    else: # packinfo
                        typenames = cl.typesDict.keys()
                        
                    for cn in typenames:
                        co = cl.typesDict[cn]
                        print "lp184:", co.fullpath
                        editlink = "/edit?%s" % urllib.urlencode(
                                                      {"filename":co.fullpath,
                                                      "lnum":0})
                                                      
                        if typeinfo or  ( packname and packname in co.fullpath ):
                            rd["types"].append(cn)
                            rd["type_meta"].update({cn:
                                                {"type":cn, 
                                                "fullpath":co.fullpath,
                                                "edit_link": editlink}
                                                   })
                        rd["types"].sort()
                if retdescriptors:
                    from astrodata.interface import Descriptors
                    das = [];
                    if packinfo:
                        for dv in Descriptors.centralCalculatorIndex:
                            mod = Descriptors.centralCalculatorIndex[dv].split(".")[0]
                            exec("import %s" % mod)
                            module = eval(mod)
                            path = module.__file__
                            # print "173:",scope, path
                            if scope in path:
                                di = {"type":dv,
                                        "descriptor_class": Descriptors.centralCalculatorIndex[dv],
                                        "path": path                             
                                     }
                                das.append(di)
                                
                    if typeinfo:
                        if astrotype in Descriptors.centralCalculatorIndex.keys():
                            mod = Descriptors.centralCalculatorIndex[astrotype].split(".")[0]
                            exec ("import %s" % mod)
                            module = eval(mod)
                            path = module.__file__
                            di = { "type":astrotype,
                                    "descriptor_class": Descriptors.centralCalculatorIndex[astrotype],
                                    "path":path
                                  }
                            das.append(di)                            
                            
                    das.sort(key= lambda da: da["type"])
                    rd.update({"descriptors":das})

                if retrecipesystem:
                    
                            
                    from recipe_system.reduction import RecipeManager
                    ri = RecipeManager.centralRecipeIndex
                    recdict = {}      
                    if "recipes" not in rd:
                        rd.update({"recipes":{}})
                    bigrecdict = rd["recipes"]
                    recipe_cat = RCat({ "cat_name": scope,
                                      })
                    
                    
                    for rec in ri.keys():
                        if typeinfo:
                                relevant = ri[rec].endswith(astrotype)
                        elif packinfo:
                            relevant = scope in ri[rec]
                            
                        print "adk625:", relevant, astrotype, ri[rec]
                        #if relevant == True:
                        #    print "pt252:", relevant, rec, scope, ri[rec]
                        if packinfo:
                            package = scope   
                        elif typeinfo:
                            package = "*"
                        if relevant == True:
                            #print "pt254: TRUE"
                            # print "pt252:", rec, scope, ri[rec]
                            trd = {"name":rec,
                                   "path":ri[rec],
                                   "package":package}
                                  
                            bigrecdict.update({rec: trd})
                            recipe_cat.add(trd)
                    rd.update({"recipe_cat":recipe_cat.cat})
                                
                if retprimsets:    
                    from recipe_system.reduction import RecipeManager
                    from recipe_system.reduction.RecipeManager import RecipeLibrary, centralPrimitivesCatalog
                    rl = RecipeLibrary()
                    if packinfo:
                        primsetkeys = RecipeManager.centralPrimitivesIndex.keys()
                    else:
                        primsetkeys = [astrotype]
                    
                    print "adkt649:", primsetkeys
                                            
                    cpi = RecipeManager.centralPrimitivesIndex
                    primsetkeys.sort()
                    primitive_cat = {}
                    for key in primsetkeys:
                        if key in cpi:
                            for entry in cpi[key]:                        
                                print "pt345:",key,entry
                                impname = os.path.splitext(entry[0])[0]
                                #[0]
                                #exec("import "+impname)
                                #fname = eval(impname+".__file__")
                                fname = "astrodata_Gemini"                            
                                #print "fname",fname
                                
                                ps = { "astrotype":key,
                                       "module": entry[0],
                                       "class":entry[1],                                                                 
                                     }
                                cpc = centralPrimitivesCatalog
                                pmd = cpc.get_primcat_dict(entry)                             
                                if pmd:
                                    ps.update({"index_path":pmd["path"],
                                                "package":pmd["package"]
                                               })
                                    # print "pt382:", pformat(ps)
                                if scope == pmd["package"] or typeinfo:  
                                    print "at527:", key                            
                                    primitive_cat.update({key:ps})
                    rd.update({"primitives_cat":primitive_cat})
                # print "lp189:", pprint.pformat(rd["recipes"]) 
                
                rd["astrotype"] = astrotype
                                     
                s.wfile.write(json.dumps(rd))                
                return
            elif command == "calculator_by_type":
                from astrodata.interface.Descriptors import centralCalculatorIndex
                from astrodata.interface.Descriptors import calculatorPackageMap
                import astrodata.utils.adinspect as adi
                print "538:"
                s.send_response(200)
                s.send_header("Content-type", "application/json")
                s.end_headers()
                    
                s.wfile.write(json.dumps(adi.get_descriptors(astrotype=scope), sort_keys=True, indent=4))
                return                 
            elif command == "primset_by_type":
                from recipe_system.reduction.priminspect import PrimInspect
                from inspect import getsourcelines, getsourcefile
                import astrodata.utils.adinspect as adi
                from recipe_system.reduction.RecipeManager import centralPrimitivesIndex as cpi
                # print "at529:", pformat(cpi[scope])   

                primsets = recipelibrary.retrieve_primitive_set(astrotype=scope)
                
                # print "at536:", pformat(primsets)
                             
                # if self.pi == None:
                #    self.pi = PrimInspect()
                #pi = self.pi 
                #mdict = pi.master_dict[scope]
                s.send_response(200)
                s.send_header("Content-type", "application/json")
                s.end_headers()
                # dstr = json.dumps(pidict["inheritance"]["GMOS"].keys(), sort_keys=True, indent=4)     
                # psObject = mdict["instance"]
                # classname = mdict["class"].__name__
                parsed = urlparse(self.path)
                parms = parse_qs(parsed.query)
                # print "at550:", pformat(parms)
                pclass = parms["primclass"][0]
                
                oprims = []
                tprimset = None
                for primset in primsets:
                    if primset.__class__.__name__ == pclass:
                        tprimset = primset
                    else:
                        oprims.append(
                            {   "primclass":pclass,
                                "astrotype":scope,
                                "link_frag": '''also assigned
                                            <i><b>%(primclass)s</b></i>
                                            (<a href="javascript:void(0)" 
                                            onclick="showPrimset('%(primclass)s','%(astrotype)s')"
                                            >detail</a>)
                                        ''' % {"primclass":primset.__class__.__name__,
                                                "astrotype":scope
                                               }
                             }
                           )
                # print "at572:", pformat(oprims)
                # which object?  one that has right class
                psObject = tprimset
                classname = psObject.__class__.__name__
                
                primslist =  adi.get_primitives(psObject) # mdict["primitives"].keys()
                exec("import "+psObject.__module__)
                tpath = eval(psObject.__module__+".__file__")
                if tpath[-4:] == ".pyc":
                    tpath = tpath[:-1]
                tpath = os.path.splitext(tpath)[0]+".py"               
                # print "405:", tpath  
                pd = {}   
                for prim in primslist:
                    co = eval("psObject."+prim)
                    gsl = getsourcelines(co)
                    fpath = getsourcefile(co)
                    if fpath[-4:] ==".pyc":
                        fpath = fpath[:-1] 
                     
                    # print "pt420:", fpath
                    
                    pd.update({prim: {"lnum":gsl[1],
                                        "path": fpath,
                                        "basename": os.path.basename(fpath),
                                        "dirname":os.path.dirname(fpath)}})   
                retdict = {"class":classname,
                           "path":tpath,
                           "prims":primslist,
                           "prim_info":pd,
                           "other_primsets":oprims
                          }
                    # print "pt549:", pformat(retdict)
                s.wfile.write(json.dumps(retdict))
                return            
            elif command == "primset_info":
                from recipe_system.reduction.priminspect import PrimInspect
                from inspect import getsourcelines
                if self.pi == None:
                    self.pi = PrimInspect()
                pi = self.pi 
                # print "388:", pformat(urlparms)
                # print "400:",pformat(pi.class2instance)
                classname = urlparms["class"][0]
                psObject = pi.class2instance[classname]
                primslist = pi.primsdict[psObject]
                exec("import "+psObject.__module__)
                tpath = eval(psObject.__module__+".__file__")                
                # print "405:", tpath     
                       
                s.send_response(200)
                s.send_header("Content-type", "application/json")
                s.end_headers()
                
                tpath = os.path.splitext(tpath)
                tpath = tpath[0]+".py"
                pd = {}
                for prim in primslist:
                    co = eval("psObject."+prim)
                    gsl = getsourcelines(co)
                    # print "pt420:", pformat(gsl)
                    pd.update({prim: {"lnum":gsl[1]}})   
                retdict = {"class":classname,
                           "path":tpath,
                           "prims":primslist,
                           "prim_info":pd
                          }
                # print "pt549:", pformat(retdict)
                s.wfile.write(json.dumps(retdict))
                return      
                
                
        def do_POST(s):
            """Respond to a POST request."""
            import cgi
            self = s
            try:
                ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))
                if ctype == 'multipart/form-data':
                    postvars = cgi.parse_multipart(self.rfile, pdict)
                elif ctype == 'application/x-www-form-urlencoded':
                    length = int(self.headers.getheader('content-length'))
                    postvars = cgi.parse_qs(self.rfile.read(length), keep_blank_values=1)
                else:
                    postvars = {}
            except:
                postvars = {}
            print "350:", pformat(postvars)
            self = s
            # If someone went to "http://something.somewhere.net/foo/bar/",
            # then s.path equals "/foo/bar/".
            #s.wfile.write("<p>You accessed path: %s</p>" % s.path)
            pres = urlparse(s.path)
            urlparms = parse_qs(pres.query)
            cmds = pres.path.split("/")[1:]
            command = scope = subscope = None
            if len(cmds)>0:
                command = cmds[0]
            if len(cmds)>1:
                scope  = cmds[1]
            if len(cmds)>2:
                subscope = cmds[2]    
            postvars.update(urlparms)
            if command == "edit":
                
                print "pt165:",s.path
                print "pt169:",postvars["filename"]
                fn = postvars["filename"][0]
                lnum = postvars["lnum"][0]
                ethread = EditThread(fn, lnum)
                ethread.start()

                s.send_response(200)
                s.send_header("Content-type", "application/json")
                s.end_headers()
                s.wfile.write(json.dumps({"status":"success"}))                
                return    
            
    server_class = BaseHTTPServer.HTTPServer
    httpd = server_class((HOST_NAME, PORT_NUMBER), IfaceHandler)
    try:
        print "Starting HTTP Interface at %s:%d" % (HOST_NAME, PORT_NUMBER)
        print "Opening Window in Default Browser"
        if False:
            # this opens a tab/window in the default browser
            import webbrowser
            webbrowser.open("http://localhost:%d/" % PORT_NUMBER, 0, True)
        
        
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    httpd.server_close()

      
      
