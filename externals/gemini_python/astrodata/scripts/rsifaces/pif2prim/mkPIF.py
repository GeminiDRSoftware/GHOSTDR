#!/usr/bin/env python
import os
import re

class PIFFunction():
    PIF_TEMPLATE = """
def %(pifname)s(*args, **argv):
    ro = mkRO(astrotype="%(astrotype)s", copy_input=%(copy_input)s, 
              args=args, argv=argv)
    ro.runstep("%(primname)s", ro.context)
    outputs = ro.context.get_outputs(style="AD")
    if len(outputs)==0:
        return None
    elif len(outputs)==1:
        return outputs[0]
    else:
        return outputs
    
    """
    
    output_filename = None
    pifdict = None
    modulestr = None
    primname = None
    basedir = ""
    def __init__(self, primname, pifdict):
        # print "25:", repr(pifdict)
        self.primname = primname
        self.pifdict = pifdict
        if "basedir" in pifdict.keys():
            self.basedir = pifdict["basedir"]
        self.prepare_output_location()
        
    @classmethod
    def camelTo_pep8(cls, name):
        # function uses "Ironic Naming" PEP8i approved
        pep8 = ""
        up = 0
        for i in range(len(name)):
            char = name[i]
            if char.isupper():
                if up == 0:
                    pep8 += "_"+char.lower()
                else:
                    pep8 += char.lower()
                up += 1
            else:
                pep8 += char
                up = 0
        return pep8
           
    def prepare_output_location(self):
        if "module" in self.pifdict:
            modulestr = self.pifdict["module"]
        else:
            modulestr = ""
        
        dirs = list(os.path.split(self.basedir))
        basedirlen = len(dirs)
        dirs.extend(modulestr.split("."))
        self.modulestr = modulestr
        for i in range(basedirlen-1, len(dirs)):
            curpath = os.path.join(*dirs[0:i+1])
            if i > 0:
                parpath = os.path.join(*dirs[0:i])
            else:
                parpath = "dirs[0]"
            # print "24:", repr(dirs), repr(curpath)
            if not os.path.exists(curpath):
                if curpath != "":
                    print "PIFFunc making %s" % curpath
                    os.mkdir(curpath);
            initpy = os.path.join(curpath,"__init__.py")
            if not os.path.exists(initpy):
                ipf = open(initpy, "w+")
                ipf.write('#test\nfrom recipe_system.reduction.mkro import *\n\n')
                ipf.close()
                
                pipf = open(os.path.join(parpath, "__init__.py"), "a+")
                pipf.write("\nimport %s\n" % dirs[i])
                pipf.close()
        self.output_filename = initpy
        
        
    def write_pif(self):
        if "pep8name" in self.pifdict.keys():
            pifname = self.pifdict["pep8name"]
        else:
            pifname = self.camelTo_pep8(self.primname)
        if "astrotype" in self.pifdict.keys():
            astrotype = self.pifdict["astrotype"]
        else:
            astrotype = "GEMINI"
        if "copy_input" in self.pifdict.keys():
            copy_input = self.pifdict["copy_input"]
        else:
            copy_input = True
                        
        self.prepare_output_location()
        print "=----> WRITING PIF %s to %s" % (self.primname, self.output_filename)
        outfile = open(self.output_filename, "a+")
        primname = self.primname
        outfile.write(self.PIF_TEMPLATE % 
                        {"primname":primname,
                         "astrotype":astrotype,
                         "pifname" :pifname,
                         "copy_input":copy_input,
                         }
                     )
        outfile.close()
        
            
class PIFGenerator():
    output_dir = "./pif"
    primdicts_dir = "primdicts"
    primdict_re = r"pif2primDict(.*?)\.py$"
    primdict = None
    primlist = None
    pif_funcs = None # to be in lists keyed by module name
               
    def __init__(self):
        self.load_primdict()
        self.pif_funcs = {}
        
            
    @classmethod
    def camelTo_pep8(cls, name):
        # function uses "Ironic Naming" PEP8i approved
        pep8 = ""
        up = 0
        for i in range(len(name)):
            char = name[i]
            if char.isupper():
                if up == 0:
                    pep8 += "_"+char.lower()
                else:
                    pep8 += char.lower()
                up += 1
            else:
                pep8 += char
                up = 0
        return pep8
   
    
    
    
    def generate(self):
        if os.path.exists(self.output_dir):
            print "\n\n"+"OUTPUT DIRECTORY EXISTS, please remove manually: %s, ./__init__.py\n" % self.output_dir *3 + "\n"
            raise self
        for prim in self.primdict:
            print "planning to generate PIF for: %s" % prim
            pifdict = self.primdict[prim]
            pifdict.update({"basedir":self.output_dir})
            pif = PIFFunction(prim, pifdict)
            if (pif.modulestr not in self.pif_funcs):
                self.pif_funcs.update({pif.modulestr:[]})
                
            pifs = self.pif_funcs[pif.modulestr]
            pifs = pifs.append(pif)
            
        for module in self.pif_funcs:
            print "MODULE:", module
            pifs = self.pif_funcs[module]
            print "PIFS:"+repr(pifs)
            for pif in pifs:
                pif.write_pif()
                
    def load_primdict(self):
        primdict = {}
        for root, dirs, files in os.walk(self.primdicts_dir):
            for potfile in files:
                if re.match(self.primdict_re, potfile):
                    
                    pdf = open(os.path.join(root, potfile))
                    pd  = eval(pdf.read())
                    print "adding %s to primdict" % potfile
                    primdict.update(pd)
                    
        self.primdict = primdict
        self.primlist = primdict.keys()
        print "loadpd: 173", self.primdict
        return 
pifgen = PIFGenerator()
pifgen.generate()
