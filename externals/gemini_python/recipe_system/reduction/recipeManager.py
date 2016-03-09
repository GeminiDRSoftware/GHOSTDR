#
#                                                                  gemini_python
#
#                                                        recipe_system.reduction
#                                                               RecipeManager.py
# ------------------------------------------------------------------------------
# $Id: recipeManager.py 5201 2015-04-06 16:28:32Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5201 $'[11:-2]
__version_date__ = '$Date: 2015-04-06 06:28:32 -1000 (Mon, 06 Apr 2015) $'[7:-2]
# ------------------------------------------------------------------------------
# This module operates like a singleton
import os
import re
import new
import sys
import traceback

from copy import copy
from datetime import datetime

from astrodata import AstroData
from astrodata.utils import logutils
from astrodata.utils import ConfigSpace

from astrodata.utils.Errors import PrimitiveError
from astrodata.utils.ConfigSpace import RECIPEMARKER

from astrodata.utils.gdpgutil import inherit_index
from astrodata.utils.gdpgutil import pick_config

from .primitivescat import PrimitivesCatalog
from .reductionObjects import ReductionObject
from .reductionObjects import ReductionError

#------------------------------------------------------------------------------
centralPrimitivesIndex = {}
centralRecipeIndex = {}
centralRecipeInfo = {}
centralReductionMap = {}
centralAstroTypeRecipeIndex = {}
centralParametersIndex = {}
centralAstroTypeParametersIndex = {}
centralPrimitivesCatalog = PrimitivesCatalog()
#------------------------------------------------------------------------------
class RecipeError(Exception):
    pass

class SettingFixedParam(RecipeError):
    pass

class RCBadParmValue(RecipeError):
    pass

#------------------------------------------------------------------------------
def get_recipe_info(recname, filename):
    rd = {}
    paths = filename.split("/")
    recipefile = paths[-1]
    for d in paths:
        if RECIPEMARKER in d:
            ind = paths.index(d)
            paths = paths[ind:-1]
            break
    recinfo = { "recipe_name":recname,
                "fullpath":filename,
                "basename":recipefile,
                "category_list":repr(paths),
                "package_path":"/".join(paths)
               }
    return recinfo

def open_if_name(dataset):
    """
    Utility function to handle accepting datasets as AstroData
    instances or string filenames. Works in conjunction with close_if_name.
    The way it works, open_if_name opens returns an AstroData instance.
    """    
    bNeedsClosing = False    
    if type(dataset) == str:
        bNeedsClosing = True
        gd = AstroData(dataset)
    elif isinstance(dataset, AstroData):
        bNeedsClosing = False
        gd = dataset
    else:
        raise RecipeError("BadArgument in recipe utility function: "+
                           "open_if_name(..)\n MUST be filename (string) "+
                           "or AstroData instance")
    return (gd, bNeedsClosing)
    
def close_if_name(dataset, b_needs_closing):
    """
    Utility function to handle accepting datasets as AstroData
    instances or string filenames. Works in conjunction with open_if_name.
    """
    if b_needs_closing == True:
        dataset.close()
    return

#------------------------------------------------------------------------------ 
class RecipeLibrary(object):
    prim_load_times = {}
    
    def add_load_time(self, source, start, end):
        key = datetime.now()
        pair = {key: {"source":source, "start":start, "end":end}}
        self.prim_load_times.update(pair)
        return

    def discover_correct_prim_type(self, context):
        ref = context.get_reference_image()
        if ref == None:
            return None
        val = pick_config(ref, centralPrimitivesIndex, "leaves")
        k = val.keys()
        if False: # Allow this. Recipes and multiple packages len(k) != 1:
            errstr = "Can't discover correct primtype for %s, more than one (%s)"
            raise RecipeError(errstr % (ref.filename, repr(k)))
        return k[0]
        
    def report_history(self):
        self.report_load_times()
        return
        
    def report_load_times(self):
        skeys = self.prim_load_times.keys()
        skeys.sort()
        
        for key in skeys:
            primrecord = self.prim_load_times[key]
            source = primrecord["source"]
            start = primrecord["start"]
            end = primrecord["end"]
            duration = end - start
            
            pargs = {   "module":source,
                        "duration":duration,
                        }
            print "Module '%(module)s took %(duration)s to load'" % pargs
        return

    def bind_recipe(self, redobj, name, recipefunc):
        rprimset = redobj.new_primitive_set(redobj.curPrimType, btype="RECIPE")
        bindstr = "rprimset.%s = new.instancemethod(recipefunc, redobj, None)" % name
        exec(bindstr)
        redobj.add_prim_set(rprimset, add_to_front = False)
        return redobj

    def load_and_bind_recipe(self, ro, name, dataset=None, astrotype=None, src=None):
        """
        Will load a single recipe, compile and bind it to the given reduction 
        objects. If src is set, dataset and astrotype are ignored (no recipe lookup)
        """
        re_msg1 = "NAME CONFLICT: ASSIGNING RECIPE %s BUT EXISTS AS PRIMITIVE:\n\t%s"
        re_msg2 = "Recipe Source for '%s' Not Found\n\ttype=%s, "+ \
                  "instruction_name=%s, src=%s"
        if src != None:
            rec = src
            # compose to python source
            prec = self.compose_recipe(name, rec)
            # compile to unbound function
            rfunc = self.compile_recipe(name, prec)
            ro = self.bind_recipe(ro, name, rfunc)
        elif astrotype != None:
            rec = self.retrieve_recipe(name, astrotype=astrotype)
            try:
                ps = ro.get_prim_set(name)
                if ps:
                    if rec == None:
                        return         #not a recipe; exists as primitive
                    else:
                        raise RecipeError(re_msg1 % rec, repr(ps))
            except ReductionError:
                pass # no primset, that function throws
                
            if rec:
                prec = self.compose_recipe(name, rec)
                # compile to unbound function
                rfunc = self.compile_recipe(name, prec)
                ro = self.bind_recipe(ro, name, rfunc)
            else:
                raise RecipeError(re_msg2 % (name, astrotype, name, src))

        elif dataset != None:
            gd, bnc = open_if_name(dataset)
            types = gd.type()
            rec = None
            for typ in types:
                rec = self.retrieve_recipe(name, astrotype=typ, inherit=False)
                if rec:
                    prec  = self.compose_recipe(name, rec)
                    rfunc = self.compile_recipe(name, prec)
                    ro = self.bind_recipe(ro, name, rfunc)
            # no recipe, see if there is a generic one
            if rec == None:
                rec = self.retrieve_recipe(name)
                if rec:
                    prec = self.compose_recipe(name, rec)
                    rfunc = self.compile_recipe(name, prec)
                    ro = self.bind_recipe(ro, name, rfunc)
            close_if_name(gd, bnc)
        return

    def get_applicable_recipes(self, dataset=None, 
                                     astrotype=None, 
                                     collate=False,
                                     prune=False):
        """
        Get list of recipes associated with all the types that apply to this 
        dataset.
        """
        if dataset != None and astrotype != None:
            raise RecipeError("dataset AND astrotype set not allowed")
        if dataset == None and astrotype == None:
            raise RecipeError("must pass dataset OR explicit astrotype")
        byfname = False
        if dataset:
            if  type(dataset) == str:
                astrod = AstroData(dataset)
                byfname = True
            elif type(dataset) == AstroData:
                byfname = False
                astrod = dataset
            else:
                raise BadArgument()
            # get the types
            types = astrod.type(prune=True)
        else:
            types = [astrotype]
        # look up recipes, fill list
        reclist = []
        recdict = {}

        for typ in types:
            if False:
                if typ in centralAstroTypeRecipeIndex.keys():
                    recnames = centralAstroTypeRecipeIndex[typ]

            recnames = inherit_index(typ, centralAstroTypeRecipeIndex)

            if recnames:
                reclist.extend(recnames[1])
                recdict.update({recnames[0]: recnames[1]})
        reclist = list(set(reclist))

        if byfname:
            astrod.close()
        
        if collate == False:
            return reclist
        else:
            return recdict
        
    def recipe_index(self, as_xml=False):
        cri = centralRecipeIndex
        
        if as_xml == False:
            return copy(cri)
        else:
            rs  = '<?xml version="1.0" encoding="UTF-8" ?>\n'
            rs += "<recipe_index>\n"
            for typ in cri.keys():
                recipe = cri[typ]
                rs += '\t<recipeAssignment type="%s" recipe="%s"/>\n' % \
                      (typ, recipe)
            rs += "</recipe_index>\n"
            return rs

    def list_recipes(self, as_xml=False):
        cri = centralRecipeIndex
        recipelist = cri.keys()
        if as_xml==True:
            retxml  = '<?xml version="1.0" encoding="UTF-8" ?>\n'
            retxml += "<recipes>\n"
            for recipe in recipelist:
                retxml += """\t<recipe name="%s" path="%s"/>\n""" % \
                (recipe, cri[recipe])
            retxml += "</recipes>\n"
            return retxml
        else:
            return recipelist

    def retrieve_recipe(self, name, astrotype=None, inherit=True):
        cri = centralRecipeIndex
        if astrotype:
            akey = name + "." + astrotype
            key = name 
        else:
            key = name
            akey = name + ".None"

        bdefRecipe = key in cri
        bastroRecipe = akey in cri
        
        fname = None
        if bastroRecipe:
            fname = cri[akey]
        elif bdefRecipe:
            if astrotype == None:
                fname = cri[key]
            else:
                # @@NOTE: OLD WAY: User must SPECIFY none to get the generic recipe
                # return None
                # @@....: new way: inherit generic recipe!
                if inherit == True:
                    fname = cri[key]
                else:
                    return None        
        else:
            return None

        rfile = file(fname, "r")
        rtext = rfile.read()
        # print "RM1433:", rtext
        return rtext
            
    def retrieve_reduction_object(self, dataset=None, astrotype=None):
        a = datetime.now()
        
        # If astrotpye is None, but dataset is set, then we need to get the 
        # astrotype from the dataset.  For reduction objects, there can be 
        # only one assigned to a real object if there are multiple reduction 
        # objects associated with type we must find out through inheritance 
        # relationships which one applies. E.g. if a dataset is GMOS_SPEC and
        # GMOS_IFU, then an inheritance relationship is sought, and the child 
        # type has priority. 
        # If they cannot be resolved, because there are unrelated types or 
        # through multiple inheritance multiple ROs may apply, then we raise 
        # an exceptions, this is a configuration problem.
        
        ro = ReductionObject()
        primsetlist = self.retrieve_primitive_set(dataset=dataset, 
                                                  astrotype=astrotype)
        ro.recipeLib = self

        if primsetlist:
            ro.curPrimType = primsetlist[0].astrotype
        else:
            return None
        for primset in primsetlist:
            ro.add_prim_set(primset)
        
        b = datetime.now()
        if astrotype != None:
            source = "TYPE: " + astrotype
        elif dataset != None:
            source = "FILE: " + str(dataset)
        else:
            source = "UNKNOWN"

        self.add_load_time(source, a, b)
        return ro
        
    def retrieve_primitive_set(self, dataset=None, astrotype=None):
        val = None
        primset = None
        primlist = []
        k = []

        if (astrotype == None) and (dataset != None):
            val = pick_config(dataset, centralPrimitivesIndex, style="leaves")
            k = val.keys()

        if (astrotype != None):
            k = [astrotype]

        for astrotype in k:
            blmsg =  "##### PRIMITIVE SET IMPORT ERROR: SKIPPING %(impname)s\n" * 3
            if (astrotype != None) and (astrotype in centralPrimitivesIndex):
                primdeflist = centralPrimitivesIndex[astrotype]
                for primdef in primdeflist:
                    rfilename = primdef[0]              # primset file
                    rpathname = centralReductionMap[rfilename]
                    rootpath = os.path.dirname(rpathname)
                    importname = os.path.splitext(rfilename)[0]
                    a = datetime.now()
                    try:
                        exec ("import " + importname)
                    except KeyboardInterrupt:
                        log = logutils.get_logger(__name__)
                        log.error("Primitive import interrupted")
                        log.error("retrieve_primitive_set() received Ctrl-C event.")
                        raise KeyboardInterrupt
                    except NameError:
                        log = logutils.get_logger(__name__)
                        blmsg = "##### PRIMITIVE SET IMPORT ERROR: SKIPPING "\
                                "%(impname)s\n" * 3
                        blmsg = blmsg % {"impname": importname}
                        msg = blmsg + traceback.format_exc() + blmsg
                        if log:
                            log.error(msg)
                        else:                    
                            print "PRINTED, not logged:\n"+msg
                        raise PrimitiveError

                    try:
                        primset = eval(importname + "." + primdef[1] + "()")
                    except NameError:
                        traceback.print_exc()
                        print "NOTE "*15
                        print "NOTE: If you have trouble importing a "
                        print "      primitive set, you may need to "
                        print "      create login.cl.  With IRAF installed,"
                        print "      this can be done with the 'mkiraf' command."
                        print "NOTE "*15
                        raise NameError
                    except Exception:
                        print
                        print ("!@"*40)
                        print "PROBLEM CREATING PRIMITIVE SET"
                        print (("!@"*40))
                        traceback.print_exc()
                        print ("!@"*40)
                        print
                        raise PrimitiveError

                    primset.astrotype = astrotype
                    primset.acquire_param_dict()
                    primlist.append(primset)

        if len(primlist):
            return primlist
        else:
            return None
        
    def compose_recipe(self, name, recipebuffer):
        templ = """
def %(name)s(self,cfgObj):
    #print "${BOLD}RECIPE BEGINS: %(name)s${NORMAL}" #$$$$$$$$$$$$$$$$$$$$$$$$$$$
    recipeLocalParms = cfgObj.localparms
%(lines)s
    #print "${BOLD}RECIPE ENDS:   %(name)s${NORMAL}" #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    yield cfgObj
"""
        newl_str = """
            
    if "%(line)s" in recipeLocalParms:
        dostep = (str(recipeLocalParms["%(line)s"]).lower() != "false")
    else:
        dostep = True
    if dostep:
        cfgObj.localparms = eval('''%(parms)s''')
        #cfgObj.localparms.update(recipeLocalParms)
        # add parms specified
        for pkey in cfgObj.localparms:
            val = cfgObj.localparms[pkey]
            if val[0]=="[" and val[-1]=="]":
                vkey = val[1:-1]
                if vkey in recipeLocalParms:
                    cfgObj.localparms[pkey] = recipeLocalParms[vkey]
        for co in self.substeps('%(line)s', cfgObj):
            if (co.is_finished()):
                break
            yield co
    yield co""" 

        recipelines = recipebuffer.splitlines()
        lines = ""
        
        for line in recipelines:
            line = re.sub("#.*?$", "",line)
            line = line.strip()
            # PARSE PRIMITIVE ARGUMENT LIST
            # take parenthesis off, make arg dict with it
            m = re.match("(?P<prim>.*?)\((?P<args>.*?)\)$", line)
            d = {}
            if m:
                prim = m.group("prim").strip()
                args = m.group("args")
                elems = args.split(",")
                for elem in elems:
                    selem = elem.strip()
                    if "=" in selem:
                        parmname, parmval = elem.split("=")
                        parmname = parmname.strip()
                        parmval = parmval.strip()
                        if parmval[0] == '"' or parmval[0] == "'":
                            parmval = parmval[1:]
                        if parmval[-1] == '"' or parmval [-1] == "'":
                            parmval = parmval[:-1]
                        d.update({parmname:parmval})
                    else:
                        if len(selem)>0:
                            d.update({selem:True})
                line = prim
            # need to add dictionary to context
            
            if line == "" or line[0] == "#":
                continue

            newl = newl_str % {"parms":repr(d),"line":line}
            lines += newl
            
        rets = templ % { "name" : name,
                         "lines" : lines,
                       }
        return rets
        
    def compile_recipe(self, name, recipeinpython):
        exec(recipeinpython)
        func = eval(name)
        return func
            
    def check_method(self, redobj, primitivename):
        ps = redobj.get_prim_set(primitivename)
        if ps == None:
            return False
        else:
            return True
        
    def check_and_bind(self, redobj, name, context=None):
        if self.check_method(redobj, name):
            return False
        else:
            self.load_and_bind_recipe(redobj, name, astrotype = redobj.curPrimType)
            return True

# ------------------------------------------------------------------------------
# CODE THAT RUNS ON IMPORT
# THIS MODULE ACTS AS A SINGLETON FOR RECIPE FEATURES

# NOTE: The issue of a central service for recipes implies a need for
# a singleton as with ClassificationLibrary and Descriptors.py.
# I have adopted the module-as-singleton approach for Structures, which 
# does not involve the message try-instantiate-except block used in the 
# ClassificationLibrary. 

#: recipeIndexREMask used to identify which files by filename
#: are those with tables relating type names to structure types

primitivesIndexREMask = r"primitivesIndex\.(?P<modname>.*?)\.py$"
parameterIndexREMask  = r"parametersIndex\.(?P<modname>.*?)\.py$"
recipeIndexREMask     = r"recipeIndex\.(?P<modname>.*?)\.py$"

#theorectically could be automatically correlated by modname

reductionObjREMask    = r"primitives_(?P<redname>.*?)\.py$"
recipeAstroTypeREMask = r"(?P<recipename>.*?)\.(?P<astrotype>.*?)$"

parameterREMask = r"parameters\.(?P<recipename>.*?)\.py$"
recipeREMask    = r"recipe\.(?P<recipename>.*?)$"

# was firstrun logic... python interpreter makes sure this module 
# only runs once already

if True: 
    # WALK the directory structure
    # add each directory to the sytem path (from which import can be done)
    # and exec the structureIndex.***.py files
    # These indices are meant to append to the centralDescriptorIndex
            
    for root, dirn, files in ConfigSpace.config_walk("recipes"):
        root = os.path.abspath(root)
        #print "RM2193:", root
        sys.path.append(root)
        curdir = root
        curpack = ConfigSpace.from_which(curdir)
        for sfilename in files:
            curpath = os.path.join(curdir, sfilename)
            m   = re.match(recipeREMask, sfilename)
            mpI = re.match(primitivesIndexREMask, sfilename)
            mri = re.match(recipeIndexREMask, sfilename)
            mro = re.match(reductionObjREMask, sfilename) 
            mpa = re.match(parameterREMask, sfilename)
            mpaI = re.match(parameterIndexREMask, sfilename)

            fullpath = os.path.join(root, sfilename)

            if m:
                # this is a recipe file
                recname = m.group("recipename")
                # For duplicate recipe names, add extras.
                if centralRecipeIndex.has_key(recname):
                    # check if the paths are really the same file
                    if os.path.abspath(fullpath) != os.path.abspath(centralRecipeIndex[recname]):

                        print "-" * 35 + " WARNING " + "-" * 35
                        print "There are two recipes with the same name."
                        print "The duplicate:"
                        print fullpath
                        print "The Original:"
                        print centralRecipeIndex[recname]
                        print
                        
                        # @@TODO: eventually continue, don't raise!
                        # This makes bad recipe packages halt the whole package!
                        # raise now because this should NEVER happen.
                        raise RecipeError("Two Recipes with the same name.")

                centralRecipeIndex.update({recname: fullpath})
                recinfo = get_recipe_info(recname, fullpath)
                centralRecipeInfo.update({recname: recinfo})               
                am = re.match(recipeAstroTypeREMask, m.group("recipename"))

            elif mpI: # this is an primitives index
                efile = open(fullpath, "r")
                exec (efile)
                efile.close()
                cpis = set(centralPrimitivesIndex.keys())
                cpi = centralPrimitivesIndex
                try:
                    lpis = set(localPrimitiveIndex.keys())
                    lpi = localPrimitiveIndex
                except NameError:
                    print "WARNING: localPrimitiveIndex not found in %s" % fullpath
                    continue

                intersect = cpis & lpis
                if intersect:
                    for typ in intersect:

                        # Allow this
                        # @@NOTE: there may be a conflict, in which case order 
                        # is used to give preference we should have a tool to 
                        # check this, because really it's only OK if none of 
                        # the members of the primitive set have the same name
                        # which we don't know until later, if we actually load 
                        # and use the primtiveset
                        if False:
                            rs = "Multiple Primitive Sets Found for Type %s" % typ
                            rs += "\n  Primitive Index Entry from %s" % fullpath
                            rs += "\n  adds ... %s" % repr(localPrimitiveIndex[typ])
                            rs += "\n  conflicts with present setting ... %s" \
                                  % repr(centralPrimitivesIndex[typ])
                            print "WARNING:\n" + rs
                        pass

                for key in lpis:
                    if key not in cpis:
                        centralPrimitivesIndex.update({key:[]})
                    
                    plist = centralPrimitivesIndex[key]
                    val = lpi[key]
                    if type(val) == tuple:
                        plist.append(localPrimitiveIndex[key])
                        centralPrimitivesCatalog.add_primitive_set(
                            package = curpack, 
                            primsetEntry = val,
                            primsetPath = curpath)
                    else:
                        plist.extend(val)
                           
            elif mro: # reduction object file ... contains primitives as members
                centralReductionMap.update({sfilename: fullpath})

            elif mri:                        # this is a recipe index
                efile = open(fullpath, "r")
                exec efile
                efile.close()

                for key in localAstroTypeRecipeIndex.keys():
                    if centralRecipeIndex.has_key(key):
                        curl = centralRecipeIndex[key]
                        curl.append(localAstroTypeRecipeIndex[key])
                        localAstroTypeRecipeIndex.update({key: curl})
                    if key in centralAstroTypeRecipeIndex:
                        ls = centralAstroTypeRecipeIndex[key]
                    else:
                        ls = []
                        centralAstroTypeRecipeIndex.update({key:ls})
                    ls.extend(localAstroTypeRecipeIndex[key])

            elif mpa:                        # Parameter file
                efile = open(fullpath, "r")
                exec(efile)
                efile.close()
                recname = mpa.group("recipename")
                centralParametersIndex.update({recname:localParameterIndex})

            elif mpaI:                       # ParameterIndex file
                efile = open(fullpath, "r")
                exec(efile)
                efile.close()
                centralAstroTypeParametersIndex.update(localparameterTypeIndex)
