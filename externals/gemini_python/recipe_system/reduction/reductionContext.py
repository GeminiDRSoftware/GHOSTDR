#
#                                                                  gemini_python
#
#                                                        recipe_system.reduction
#                                                            reductionContext.py
# ------------------------------------------------------------------------------
# $Id: reductionContext.py 5418 2015-12-02 13:26:18Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5418 $'[11:-2]
__version_date__ = '$Date: 2015-12-02 03:26:18 -1000 (Wed, 02 Dec 2015) $'[7:-2]
# ------------------------------------------------------------------------------
import os
import re
import sys
import socket
import pickle

from copy import copy
from datetime import datetime

from astrodata import AstroData
from astrodata.utils import gdpgutil
from astrodata.utils.Errors import ReduceError

import IDFactory
import eventsManager

from .recipeManager import RecipeError
from .recipeManager import RCBadParmValue
from .recipeManager import SettingFixedParam

from .reductionContextRecords import AstroDataRecord
from .reductionContextRecords import CalibrationRecord
from .reductionObjectRequests import UpdateStackableRequest
from .reductionObjectRequests import GetStackableRequest
from .reductionObjectRequests import DisplayRequest
from .reductionObjectRequests import ImageQualityRequest

from .StackKeeper import StackKeeper
from .StackKeeper import FringeKeeper

from ..cal_service.CalibrationDefinitionLibrary import CalibrationDefinitionLibrary
#------------------------------------------------------------------------------
MAINSTREAM = "main"
#------------------------------------------------------------------------------

class ReductionContext(dict):
    """The ReductionContext is used by primitives and recipies, implicitely in
    the later case, to get input and report output. This allows primitives to be
    controlled in many different running environments, from pipelines to command
    line interactive reduction.
    
    :sort: __init__,__contains__,__getitem__,__str__,
        _*,a*,b*,c*,d*,e*,f*,g*,h*,i*,j*,k*,l*,m*,n*,o*,p*,q*,r*,
        s*,t*,u*,v*,w*,x*,y*,z*
    """
    status      = "EXTANT"
    reason      = "EXTANT"
    cmd_request = "NONE"

    ro = None       # What is this for?
    
    def __init__(self, adcc_mode="start_early"):
        """
        Constructor creates empty dictionaries and lists, 
        members set to None in the class.
        """
        self._adcc_mode      = adcc_mode
        self._localparms     = None
        self._current_stream = MAINSTREAM
        self._output_streams = []
        self._metricEvents = eventsManager.EventsManager(rc = self)
        self._nonstandard_stream = []
        self._running_contexts = None

        self.arguments    = []

        self.cache_files  = {}
        self.calibrations = {}
        self.callbacks    = {}
        self.calindfile   = None
        self.cmd_history  = []
        self.cmd_index    = {}

        self.display_id   = None
        self.display_mode = None
        self.display_name = None

        self.hostname     = socket.gethostname()

        self.indent       = 0
        self.inputs       = []
        self.inputs_history = []

        self.irafstdout   = None
        self.irafstderr   = None

        self.outputs      = {"main":[]}
        self.proxy_id     = 1            # used to ensure uniqueness

        self.rorqs        = []           # terrible name
        self.stephistory  = {}

        self.user_params  = None         # UserParams() instance

        # return behavior
        self._return_command = None
        self._stream_meta    = {}

        # TESTING
        self.cdl = CalibrationDefinitionLibrary()
        self.fringes = FringeKeeper()

        # Stack Keep is a resource for all RecipeManager functions... 
        # one shared StackKeeper to simulate the shared ObservationService
        # used in PRS mode.
        self.stackKeep = StackKeeper(local=True if 
                                     adcc_mode == "start_lazy" 
                                     else False)
        
    def __getitem__(self, arg):
        """Note, the ReductionContext version of __getitem__ returns None
        instead of throwing a KeyError.
        """
        if self.localparms and arg in self.localparms:
            value = self.localparms[arg]
        else:
            try:
                value = dict.__getitem__(self, arg)
            except KeyError:
                return None
        if value is None:
            retval = None
        else:
            retval = self.convert_parm_to_val(arg, value)
        return retval
    
    def __contains__(self, thing):
        """
        The ``__contains__`` function implements the Python ``in`` operator. 
        The ``ReductionContext`` is a subclass of a ``dict``, but it also has 
        a secondary dict of "local parameters" which are available to the 
        current primitive only, which are also tested by the ``__contains__(..)`` 
        member. These parameters will generally be those passed in as arguments
        to a primitive call from a recipe.

        :param thing: key to check for presence in the Reduction Context
        :type thing: <str>

        """
        if thing in self._localparms:
            return True
        return dict.__contains__(self, thing)
        
    def get_context(self):
        return ":".join(self._running_contexts)
    getContext = get_context
    
    context = property(getContext)
    
    def in_context(self, context):
        context = context.lower()
        return context in self._running_contexts
    inContext = in_context
    
    def add_context(self, context):
        context = context.lower()
        if context not in self._running_contexts:
            self._running_contexts.append(context)
    addContext = add_context        
    
    def set_context(self, context):
        if type(context) == list:
            context = [cstr.lower() for cstr in context]
            self._running_contexts = context
        else:
            self._running_contexts = [context.lower()]
    setContext = set_context
    
    def clear_context(self, context = None):
        if context is None:
            self._running_contexts = []
        else:
            context = context.lower()
            if context in self._running_contexts:
                self._running_contexts.remove(context)
    clearContext = clear_context
    
    def convert_parm_to_val(self, parmname, value):
        legalvartypes = ["bool", "int", "str", "float", None]
        vartype = self.ro.parameter_prop( parmname, prop="type")
        if vartype not in legalvartypes:
            mes =  "TEMPORARY EXCEPTION: illegal type in parameter defintions"
            mes += " for %s." % str(value)
            raise RCBadParmValue(mes)
            return value
        if vartype:
            # bool needs special handling
            if vartype == "bool":
                if type(value) == str:
                    if (value.lower() == "true"):
                        value = True
                    elif (value.lower() == "false"):
                        value = False
                    else:
                        mes = "%s is not legal boolean setting " % value
                        mes += 'for "boolean %s"' % parmname
                        raise RCBadParmValue(mes)
            retval = eval("%s(value)"%(vartype))
        else:
            retval = value
        return retval
    
    def parm_dict_by_tag(self, primname, tag, **otherkv):
        rd = self.ro.parm_dict_by_tag(primname, tag)
        rd.update(otherkv)
        return rd

    def add_cal(self, data, caltyp, calname, timestamp=None):
        '''
        :param data: Path or AstroData to which calibration will be applied.
        :type data: str or AstroData instance
        
        :param caltyp: calibration type. Eg., 'bias' and 'flat'.
        :type caltyp: str
        
        :param calname: URI of  MEF calibration file.
        :type calname: str
        
        :param timestamp: Timestamp for when calibration was added.
            for is that of datetime.datetime.
        :type timestamp: str
        
        Add a calibration to the calibration index with a key related to the
        dataset's "datalabel", so it will apply, generally to later, processed
        versions of the dataset, and thus allow retrieval of the same calibration.
        
        '''
        adID = IDFactory.generate_astro_data_id(data)
        calname = os.path.abspath(calname)
        
        if timestamp is None:
            timestamp = datetime.now()
        else:
            timestamp = timestamp
        if self.calibrations is None:
            self.calibrations = {}
        if isinstance(data, AstroData):
            filename = data.filename
        else:
            filename = data
        calrec = CalibrationRecord(filename, calname, caltyp, timestamp)
        key = (adID, caltyp)
        self.calibrations.update({key: calrec})
        return
    
    def add_callback(self, name, function):
        callbacks = self.callbacks
        if name in callbacks:
            l = callbacks[name]
        else:
            l = []
            callbacks.update({name:l})
        l.append(function)
        return
    
    def clear_input(self):
        self.inputs = []
        return

    def add_input(self, filenames):
        """
        Add input to be processed the next batch around. If this is the first
        input being added, it is also added to original_inputs.
        
        :parameter filenames: Inputs you want added.
        :type filenames: <list>, <str>, <AstroData>

        """
        if type(filenames) != list:
            filenames = [filenames]
        
        for filename in filenames:
            if filename is None:
                continue
            elif type(filename) == str:
                filename = AstroDataRecord(filename)  
            elif type(filename) == AstroData:
                filename = AstroDataRecord(filename)
            elif type(filename) == AstroDataRecord:
                pass
            else:
                m = "BadArgument: '%(name)s' is an invalid type '%(type)s'." \
                    % {'name':str(filename), 'type':str(type(filename))}
                m += "Should be str, AstroData, AstroDataRecord."
                raise ReduceError(m) 
            
            #@@CONFUSING: filename here is an AstroDataRecord!
            if filename not in self.inputs:
                self.inputs.append(filename)
        return

    def add_rq(self, rq):
        """
        Add a request to be evaluated by the control loop.
        
        @param rq: The request.
        @type rq: ReductionObjectRequests instance

        """
        self.rorqs.append(rq)
        return
        
    def begin(self, stepname):
        key = datetime.now()
        # value = dictionary
        val = self.step_moment(stepname, "begin")
        self.indent += 1
        self.stephistory.update({key: val}) 
        self.lastBeginDt = key
        self.initialize_inputs()
        return self
        
    def get_begin_mark(self, stepname, indent=None):
        for time in self.stephistory.keys():
            if self.stephistory[time]["stepname"] == stepname and \
            self.stephistory[time]["mark"] == "begin":
                if indent != None:
                    if self.stephistory[time]["indent"] == indent:
                        return (time, self.stephistory[time])
                else:
                    return (time, self.stephistory[time])    
        return None
    
    def cal_filename(self, caltype):
        """
        Returns a local filename for a retrieved calibration

        """
        retl = {}

        for inp in self.get_inputs_as_astrodata():
            key = (idFac.generate_astro_data_id(inp.ad), caltype)
            calfile = self.calibrations[key].filename
            infile = os.path.basename(inp.filename)
            if retl.has_key(calfile):
                retl.update({calfile:retl[calfile] + [infile]})
            else:
                retl.update({calfile:[infile]})
        return retl
                     
    def call_callbacks(self, name, **params):
        callbacks = self.callbacks
        if name in callbacks:
            for f in callbacks[name]:
                f(**params)
        return
                    
    def cal_summary(self, mode="text"):
        rets = ""
        for key in self.calibrations.keys():
            rets += str(key)
            rets += str(self.calibrations[key])
        return rets
    
    def check_control(self):
        return self.cmd_request
    
    def clear_rqs(self, rtype=None):
        '''
        Clear all requests.
        '''
        if rtype is None:
            self.rorqs = []
        else:
            rql = copy(self.rorqs)
            for rq in rql:
                if type(rq) == type(rtype):
                    self.rorqs.remove(rq)
        return
    
    def control(self, cmd="NONE"):
        self.cmd_request = cmd
        return
    
    def end(self, stepname):
        key = datetime.now()
        self.indent -= 1
        val = self.step_moment(stepname, "end")
        # this step saves inputs
        self.stephistory.update({key: val})
        # this step moves outputs[MAINSTREAM] to inputs and clears outputs
        self.finalize_outputs()
        self.localparms = None
        self._output_streams = []
        return self
    
    def finalize_outputs(self):
        """ 
        This function means there are no more outputs, generally called
        in a control loop when a generator function primitive ends.  Standard
        outputs become the new inputs. Calibrations and non-standard output
        is not affected.

        """
        # only push if outputs is filled
        if len(self.outputs[self._current_stream]) != 0:
            newinputlist = []
            # this code to be executed in initialize_inputs
            
    def initialize_inputs(self):
        newinputlist = []
        # print "initialize_inputs"
        for out in self.outputs[self._current_stream]:
            if type(out) == AstroDataRecord:
                newinputlist.append(out)
            else:
                mes = "Bad Argument: Wrong Type '%(val)s' '%(typ)s'." \
                    % {'val':str(out), 'typ':str(type(out))}
                raise RuntimeError(mes)
        self.inputs = newinputlist
    
    # finish and is_finished is combined using property
    def is_finished(self, arg=None):
        if arg is None:
            return self.status == "FINISHED"
        else:
            if arg == True:
                self.status = "FINISHED"
            elif self.status != "FINISHED":
                mes = "Attempt to change status from %s to FINISHED" % \
                    self.status
                raise ReduceError(mes)
        return self.is_finished()

    def finish(self):
        self.is_finished(True)
    finished = property(is_finished, is_finished)
    
    def get_cal(self, data, caltype):
        """
        Retrieve calibration.
        
        :param data: File for which calibration must be retrieved.
        :type data: string or AstroData instance
        
        :param caltype: The type of calibration (ex.'bias', 'flat').
        :type caltype: string
        
        :return: The URI of the currently stored calibration or None.
        :rtype: string or None.

        """
        adID = IDFactory.generate_astro_data_id(data)
        key = (adID, caltype)
        if key in self.calibrations.keys():
            return self.calibrations[(adID, caltype)].filename
        return None
    
    def get_end_mark(self, stepname, indent=None):
        for time in self.stephistory.keys():
            if     self.stephistory[time]["stepname"] == stepname \
               and self.stephistory[time]["mark"] == "end":
                if indent != None:
                    if self.stephistory[time]["indent"] == indent:
                        return (time, self.stephistory[time])
                else:
                    return (time, self.stephistory[time])
        return None    
    
    def get_inputs(self, style=None):
        """
        :param style: Controls the type of return value. Supported values are 
        "AD" and "FN" for ``AstroData`` and ``string`` filenames respectively.
        :type style: <str>

        :return: a list of ``AstroData`` instances or ``string`` filenames
        :rtype: <list>
        
        ``get_inputs(..)`` gets the current input datasets from the current 
        stream. You cannot choose the stream, use ``get_stream(..)`` for that.
        To report modified datasets back to the stream use ``report_output(..)``.
        """
        if style is None:
            return self.inputs
        elif style == "AD": #@@HARDCODED: means "as AstroData instances"
            retl = []
            for inp in self.inputs:
                if inp.ad is None:
                    inp.load()
                retl.append(inp.ad)
            return retl
        elif style == "FN": #@@HARDCODED: means "as Filenames"
            retl = [inp.filename for inp in self.inputs]
            return retl
        else:
            return None # this should not happen, but given a mispelled style arg

    def get_outputs(self, style = None):
        return self.get_stream(style = style, stream = MAINSTREAM, empty = False)
        
    def get_stream(self, stream=MAINSTREAM, empty=False, style = None):
        """
        returns a list of AstroData instances in the specified stream.

        :param stream: A string name for the stream in question.  
                       To use the standard stream do not set.
        :type stream: <str>

        :param empty: Controls if the stream is emptied, defaults to "False".
        :type empty: <bool>

        :param style: controls the type of output. "AD" directs the function
            to return a list of AstroData instances. "FN" directs it to return 
            a list of filenames. If left blank or set to ``None``, the 
            AstroDataRecord structures used by the Reduction Context will be 
            returned.

        :returns: a list of <AstroDataRecord> objects, <AstroData> objects or 
                  filenames.
        :rtype: <list>

        """
        if stream in self.outputs:
            outputs = self.outputs[stream]
        else:
            return None
        if empty:
            self.outputs.update({stream:[]})
        
        if style is None:
            return outputs
        elif style == "AD":
            retl = []
            for adrec in outputs:
                if not adrec.is_loaded():
                    adrec.load()
                retl.append(adrec.ad)
            return retl
        elif style == "FN":
            retl = [ad.filename for ad in outputs]
            return retl
        else:
            raise ReduceError("get_outputs: BAD STYLE ARGUMENT")
            
    def get_inputs_as_astrodata(self):
        """
        This function is equivalent to get_inputs(style="AD")
        """
        return self.get_inputs(style="AD")
    get_inputs_as_astro_data = get_inputs_as_astrodata
    
    def get_inputs_as_filenames(self):
        """
        This function is equivalent to get_inputs(style="FN")
        """
        return self.get_inputs(style="FN")
           
    def get_iraf_stderr(self):
        if self.irafstderr != None:
            return self.irafstderr
        else:
            return sys.stderr
        
    def get_iraf_stdout(self):
        if self.irafstdout != None:
            return self.irafstdout
        else:
            return sys.stdout
        
    def get_reference_image(self):
        """
        This function returns the current reference image.  At the moment
        this is simply the first dataset in the current inputs.  However,
        use of this function allows us to evolve our concept of reference
        image for more complicated cases where the choice of a "reference" image
        may need to be different (e.g. require some data analysis to determine).

        """
        if len(self.inputs) == 0:
            return None
        if self.inputs[0].ad is None:
            return None
        return self.inputs[0].ad
    
    def get_stack_ids(self):
        cachefile = self.get_cache_file("stackIndexFile")
        retval = self.stackKeep.get_stack_ids(cachefile )
        return retval
    
    def populate_stream(self, infiles, stream=None, load = True):
        self.report_output(infiles, stream = stream, load = load)
        
        if stream is None:
            stream = self._current_stream
        if stream in self._output_streams:
            self._output_streams.remove(stream)
        return
    
    def get_stack(self, purpose=""):
        sidset = set()
        purpose = self["purpose"]
        if purpose is None:
            purpose = ""
        
        # Get ID for all inputs
        for inp in self.inputs:
            sidset.add(purpose+IDFactory.generate_stackable_id(inp.ad))

        wholelist = []
        for sid in sidset:
            stacklist = self.get_list(sid)  #.filelist
            wholelist.extend(stacklist)
        return wholelist
        
    def get_list(self, id):
        """
        :param id: Lists are associated with arbitrary identifiers,
            passed as strings.  See ``IDFactory`` for IDs built from
            standard ``astrodata`` characteristics.
        :type id: str
        
        The list functionality allows storing dataset names in a list
        which is shared by all instances of reduce running in a given
        directory.  The list is kept by an ``adcc`` instance in charge of that
        sub-directory.  The ``get_list(..)`` function retrieves a list that has
        already been requested via ``rq_stack_get(..)`` which initiates the
        interprocess request.
        
        This function does not block, and if the stack was not requested
        prior to a yeild, prior to this call, then None or an out of date
        version of this list will be retrieved.
        
        :note: "get_stack" calls get_list but takes a "purpose" to which it adds
               a stackingID as a suffix to the list identifier.
        
        """
        cachefile = self.get_cache_file("stackIndexFile")
        retval = self.stackKeep.get(id, cachefile )
        return retval
    
    def inputs_as_str(self, strippath=True):
        if self.inputs is None:
            return ""
        else:
            inputlist = []
            for inp in self.inputs:
                if inp.ad != None:
                    inputlist.append(inp.ad.filename)
                else:
                    inputlist.append(inp.filename)
            if strippath == False:
                return ",".join(inputlist)                
            else:
                return ",".join([os.path.basename(path) for path in inputlist])

    def localparms_set(self, lpd):
        self._localparms = lpd
        
    def localparms_get(self):
        if self._localparms is None:
            self._localparms = {}
        return self._localparms 
    localparms = property(localparms_get, localparms_set)
    
    def make_inlist_file(self, filename, filelist):
        try:
            fh = open(filename, 'w')
            for item in filelist:
                fh.writelines(item + '\n')
        except:
            raise "Could not write inlist file for stacking." 
        finally:
            fh.close()
        return "@" + filename
                
    def parameter_collate(self, astrotype, primset, primname):
        """
        This method looks at the default primset paramaters for primname
        and sets the localparms member.
        """
        if primname in primset.param_dict:
            correctUPD = None
            if self.user_params and not self.user_params.is_empty():
                correctUPD = self.user_params.get_user_param(astrotype, primname)
                if correctUPD != None:
                    for param in correctUPD.keys():
                        if param in self.localparms:
                            exs  = "User attempting to override parameter set "
                            exs += "in recipe\n\tastrotype = %s\n" % astrotype
                            exs += "\tprimitive = %s\n" % primname
                            exs += "\tparameter = %s\n" % str(param)
                            exs += "\t\tattempt to set to = %s\n" % \
                                correctUPD[param]
                            exs += "\t\trecipe setting = %s\n" % \
                                self.localparms[param]
                            raise SettingFixedParam(exs)
                            
            # use primset.param_dict to update self.localparms
            for param in primset.param_dict[primname].keys():
                if param in self.localparms: #  or param in self:
                    repOvrd = ("recipeOverride" not in primset.param_dict[primname][param])\
                                 or primset.param_dict[primname][param]["recipeOverride"]
                    if not repOvrd:
                        exs =  "Recipe attempts to set fixed parameter\n"
                        exs += "\tastrotype = %s\n" % astrotype
                        exs += "\tprimitive = %s\n" % primname
                        exs += "\tparameter = %s\n" % str(param)
                        exs += "\t\tattempt to set to = %s\n" % self[param]
                        exs += "\t\tfixed setting = %s\n" %  \
                            primset.param_dict[primname][param]["default"]
                        raise SettingFixedParam(exs)
                if param not in self.localparms and param not in self:
                    if "default" in primset.param_dict[primname][param]:
                        self.localparms.update(
                            {param:primset.param_dict[primname][param]["default"]}
                        )

            # about to add user paramets... some of which may be in the global
            # context (and not in correct UPD), strictly speaking these may 
            # not have been added by the user but we consider it user space
            # and at any rate expect it to not be overrided by ANY means (we
            # may want a different flag than userOverride
            for param in primset.param_dict[primname].keys():
                userOvrd = ("userOverride" not in primset.param_dict[primname][param])\
                             or primset.param_dict[primname][param]["userOverride"]
                
                if dict.__contains__(self, param):
                    
                    # note: if it's in self.localparms, that's due to legal
                    # behavior above... primitives parameters (as passed in
                    # recipes) are always added to the localparms space
                    # thus, if a value is in the main context, it MUST be
                    # userOverridable
                    if not userOvrd:
                        exs =  "\nParm set in context when userOverride is False\n"
                        exs += "\tastrotype = %s\n" % astrotype
                        exs += "\tprimitive = %s\n" % primname
                        exs += "\tparameter = %s\n" % str(param)
                        exs += "\t\tattempt to set to = %s\n" % self[param]
                        exs += "\t\tfixed setting = %s\n" % \
                            primset.param_dict[primname][param]["default"]
                        raise SettingFixedParam(exs, astrotype=astrotype)
            
            # users override everything else if  it gets here... and is allowed
            if correctUPD:
                for param in correctUPD:
                    userOvrd = ("userOverride" not in primset.param_dict[primname][param])\
                                 or primset.param_dict[primname][param]["userOverride"]
                    if param in self.localparms or param in self:
                        
                        if not userOvrd:
                            exs =  "User attempted to set fixed parameter\n"
                            exs += "\tastrotype = %s\n" % astrotype
                            exs += "\tprimitive = %s\n" % primname
                            exs += "\tparameter = %s\n" % str(param)
                            exs += "\t\tattempt to set to = %s\n" % correctUPD[param]
                            exs += "\t\tfixed setting = %s\n" % \
                                primset.param_dict[primname][param]["default"]
                            raise SettingFixedParam(exs, astrotype = astrotype)
                        else:
                            self.localparms.update({param:correctUPD[param]})
                    else:
                        self.localparms.update({param:correctUPD[param]})    
        return
    
    def param_names(self, subset = None):
        if subset == "local":
            return self.localparms.keys()
        else:
            lpkeys = set(self.localparms.keys())
            rckeys = set(self.keys())
            retl = list(lpkeys | rckeys)
            return retl
                                                                
    def outputs_as_str(self, strippath=True):
        if self.outputs is None:
            return ""
        else:
            outputlist = []
            for inp in self.outputs['main']: 
                outputlist.append(inp.filename)

            if strippath == False:
                return ", ".join(outputlist)
            else:
                return ", ".join([os.path.basename(path) for path in outputlist])
    
    def run(self, stepname):
        """ 
        :param stepname: The primitive or recipe name to run. Note: this is 
                actually compiled as a recipe. Proxy recipe names may appear
                in the logs.
        :type stepname: <str>
            
        The ``run(..)`` function allows a primitive to use the reduction
        context to execute another recipe or primitive.

        """
        a = stepname.split()
        cleanname = ""
        for line in a:
            cleanname = re.sub(r'\(.*?\).*?$', '', line)
            cleanname = re.sub(r'#.*?$', '', cleanname)
            if line != "":
                break;
        # cleanname not used!
        name = "proxy_recipe%d" % self.proxy_id
        self.proxy_id += 1
        self.ro.recipeLib.load_and_bind_recipe(self.ro, name, src=stepname)
        ret = self.ro.runstep(name, self)
        self.initialize_inputs()
        return ret

    # ------------------ PAUSE -------------------------------------------------
    def is_paused(self, bpaused=None):
        if bpaused is None:
            return self.status == "PAUSED"
        else:
            if bpaused:
                self.status = "PAUSED"
            else:
                self.status = "RUNNING"
        
        return self.is_paused()
    
    def pause(self):
        self.call_callbacks("pause")
        self.is_paused(True)
    
    def unpause (self):
        self.is_paused(False)
    paused = property(is_paused, is_paused)
    
    def request_pause(self):
        self.control("pause") 
    
    def pause_requested(self):
        return self.cmd_request == "pause"
    # ------------------ PAUSE -------------------------------------------------

    def report(self, report_history=False, internal_dict=False, context_vars=False,
                report_inputs=False, report_parameters=False, showall=False):
        """
        Prints out a report of the contents of the context object 
        
        @return: The formatted message for all the current parameters.
        @rtype: str
        """
        if showall:
            context_vars      = True
            internal_dict     = True
            report_history    = True
            report_inputs     = True
            report_parameters = True
    
        rets = "\n\n" + "-" * 50  + "\n\n\n"
        rets += " "*11 + "C O N T E X T  R E P O R T\n\n\n"
        rets += "-" * 50 + "\n" 
        varlist = [
            "inputs_history", "calibrations", "rorqs", "status", 
            "reason", "cmd_request", "hostname", "display_name",
            "stackKeep", "calindfile", "display_mode", "display_id", 
            "callbacks", "arguments", "_nonstandard_stream",
            "_current_stream", "proxy_id", "cmd_history","cmd_index" 
        ]
        
        # add in vars to show they are not there
        if not self._localparms:
            varlist.append("_localparms")
        if not self.stephistory:
            varlist.append("stephistory")
        if not self.user_params:
            varlist.append("user_params")
        if not self.inputs:
            varlist.append("inputs")

        if report_inputs:
            # inputs
            if self.inputs:
                rets += "\nInput (self.inputs)"
                for rcr in self.inputs:
                    if isinstance(rcr, AstroDataRecord):
                        rets += "\n    AstroDataRecord:"
                        rets += str(rcr)
                        rets += rcr.ad._infostr()
        
        if context_vars:
            # context vars
            rets += "\nContext Variables (self.<var>):"
            varlist.sort()
            for var in varlist:
                rets += "\n    %-20s = %s" % (var, eval("self.%s" % var ))
        
        if report_parameters:
            # _localparms
            if self._localparms:
                rets += "\n\nLocal Parameters (self._localparms)"
                pkeys = self._localparms.keys()
                pkeys.sort
                for pkey in pkeys:
                    rets += "\n    %-13s : %s" % \
                        (str(pkey), str(self._localparms[pkey]))
            
            # user params (from original varlist)
            if self.user_params:
                rets += "User Parameters:"
                rets += repr(self.user_params.user_param_dict)
            rets += "\n"

        if self.stephistory and report_history == True:
            # stephistory
            
            rets += "\n\nStep History (self.stephistory):\n"
            rets += "    " + "-"*41 + "\n\n"
            rets += "Feature deprecated until memory issue is resolved \
(callen@gemini.edu)"
            
            """
            shkeys = self.stephistory.keys()
            shkeys.sort()
            count = 0
            for key in shkeys:
                sh_dict = self.stephistory[key]
                rets += "\n" + "          S T E P " + str(count+1) 
                rets += ": " + sh_dict["stepname"] 
                rets += "\n\n\n    " + "-"*41
                rets += "\n    " + str(key) + ":"
                sh_dictkeys = sh_dict.keys()
                sh_dictkeys.sort()
                if sh_dict.has_key("inputs"):
                    rets += "\n%s%-10s : " % (" "*8, "inputs (self.inputs)")
                    for rcr in sh_dict["inputs"]:
                        if isinstance(rcr, \
                recipe_system.reduction.reductionContextRecords.AstroDataRecord):
                            rets += "'%s':\n\n    %s" % \
                (str(jkey), "reductionContextRecords.AstroDataRecord:")
                            rets += str(rcr)
                        else:
                            rets += str(rcr)
                    sh_dictkeys.remove("inputs")
                for ikey in sh_dictkeys:
                    if ikey != "outputs":
                        rets += "\n%s%-10s : %s" % \
                           (" "*8, str(ikey), sh_dict[ikey])
                if sh_dict.has_key("outputs"):
                    rets += "\n%s%-10s : " % (" "*8, "outputs (self.outputs)")
                    outputs_dict = sh_dict["outputs"]
                    outputs_dictkeys = outputs_dict.keys()
                    outputs_dictkeys.sort()
                    for jkey in outputs_dictkeys:
                        for rcr in outputs_dict[jkey]:
                            if isinstance(rcr, \
                    recipe_system.reduction.reductionContextRecords.AstroDataRecord):
                                rets += "'%s':\n\n    %s" % \
                    (str(jkey), "reductionContextRecords.AstroDataRecord:")
                                rets += str(rcr)
                                rets += "\n    OUTPUT AD.INFO (rcr.ad._infostr())"
                                rets += rcr.ad._infostr() + "\n"
                            else:
                                rets += str(rcr)
                rets += "    " + "-"*41 + "\n\n"
                count += 1"""

        if internal_dict:
            # internal dictionary contents
            cokeys = self.keys()
            rets += "\n       I N T E R N A L  D I C T I O N A R Y\n"
            loglist = []
            cachedirs = []
            others = []
            for key in cokeys:
                if key  == "cachedict":
                    rets += "\nCached Files (self[cachedict:{}]):\n"
                    cache_dict = self[key]
                    cdkeys = cache_dict.keys()
                    cdkeys.remove("storedcals")
                    cdkeys.remove("reducecache")
                    cdkeys.sort()
                    for ikey in cdkeys:
                        dirfiles = os.listdir(cache_dict[ikey])
                        if len(dirfiles) == 0:
                            dirfiles = "None"
                        rets += "    %-20s : %s\n" % (ikey, dirfiles)
                elif key[:3] == "log":
                    loglist.append(key)
                elif key == "reducecache" or key[:9] == "retrieved" or \
                    key[:6] == "stored":
                    cachedirs.append(key)
                else:
                    others.append(key)

            rets += "\nCache Directories (self[<dir>]):"
            for dir_ in cachedirs:
                rets +="\n    %-20s : %s" % (dir_, str(self[dir_]))
            rets += "\n\nLogger Info (self[<log...>]):"
            for l in loglist:
                rets +="\n    %-20s : %s" % (l, str(self[l]))
            if len(others) > 0:
                rets += "\nOther (self[<Other>]):\n"
                for o in others:
                    rets +="\n    %-20s : %s" % (o, str(self[o]))
                
        rets += "\n\n" + "-" * 50  + "\n"
        return rets
    
    def persist_cal_index(self, filename = None, newindex = None):
        if newindex != None:
            self.calibrations = newindex
        try:
            pickle.dump(self.calibrations, open(filename, "w"))
            self.calindfile = filename
        except:
            print "Could not persist the calibration cache."
            raise 
    
    def persist_fringe_index(self, filename):
        try:
            pickle.dump(self.fringes.stack_lists, open(filename, "w"))
        except:
            raise 'Could not persist the fringe cache.'
        return
            
    def persist_stk_index(self, filename):
        self.stackKeep.persist(filename)
        return
    
    def prepend_names(self, prepend, current_dir=True, filepaths=None):
        '''
        :param prepend: The string to be put at the front of the file.
        :type prepend: string
        
        :param current_dir: Used if the filename (astrodata filename) is in the
                            current working directory.
        :type current_dir: boolean
        
        :param filepaths: If present, these file paths will be modified, otherwise
                          the current inputs are modified.
        :type filepaths:
        
        :return: List of new prepended paths.
        :rtype: list  
        
        Prepends a prefix string to either the inputs or the given list of filenames.
        
        '''
        retlist = []
        if filepaths is None:
            dataset = self.inputs
        else:
            
            dataset = filepaths
            
        for data in dataset:
            parent = None
            if type(data) == AstroData:
                filename = data.filename
            elif type(data) == str:
                filename = data
            elif type(data) == AstroDataRecord:
                filename = data.filename
                parent = data.parent
            else:
                raise RecipeError("BAD ARGUMENT: '%(data)s'->'%(type)s'" % 
                                  {'data':str(data), 'type':str(type(data))})
               
            if current_dir == True:
                root = os.getcwd()
            else:
                root = os.path.dirname(filename)

            bname = os.path.basename(filename)
            prependfile = os.path.join(root, prepend + bname)
            if parent is None:
                retlist.append(prependfile)
            else:
                retlist.append((prependfile, parent))
        
        return retlist
    
    def print_headers(self):
        for inp in self.inputs:
            if type(inp) == str:
                ad = AstroData(inp)
            elif type(inp) == AstroData:
                ad = inp
            try:
                outfile = open(os.path.basename(ad.filename) + ".headers", 'w')
                for ext in ad.hdulist:
                    outfile.write("\n" + "*" * 80 + "\n")
                    outfile.write(str(ext.header))
            except:
                raise "Error writing headers for '%{name}s'." % {'name':ad.filename}
            finally:
                outfile.close()
        return
    
    def process_cmd_req(self):
        if self.cmd_request == "pause":
            self.cmd_request = "NONE"
            self.pause()
        return
            
    def remove_callback(self, name, function):
        if name in self.callbacks:
            if function in self.callbackp[name]:
                self.callbacks[name].remove(function)
        else:
            return
    
    def report_history(self):
        sh = self.stephistory
        ks = self.stephistory.keys()
        ks.sort()
        
        lastdt = None
        startdt = None
        enddt = None

        retstr = "RUNNING TIMES\n"
        retstr += "-------------\n"
        for dt in ks:
            indent = sh[dt]["indent"]
            indentstr = "".join(["  " for i in range(0, indent)])
            mark = sh[dt]["mark"]
            if mark == "begin":
                elapsed = ""
                format = "%(indent)s%(stepname)s begin at %(time)s"
            elif mark == "end":
                elapsed = "(" + str(dt - lastdt) + ") "
                format = "\x1b[1m%(indent)s%(stepname)s %(elapsed)s \x1b[22mends at %(time)s"
            else:
                elapsed = ""
                format = "%(indent)s%(stepname)s %(elapsed)s%(mark)s at %(time)s"
                
            lastdt = dtpostpend
            if startdt is None:
                startdt = dt

            pargs = {  "indent":indentstr,
                        "stepname":str(sh[dt]['stepname']),
                        "mark":str(sh[dt]['mark']),
                        "inputs":str(",".join(sh[dt]['inputs'])),
                        "outputs":str(sh[dt]['outputs']),
                        "time":str(dt),
                        "elapsed":elapsed,
                        "runtime":str(dt - startdt),
                    }
            retstr += format % pargs + "\n"
            retstr += "%(indent)sTOTAL RUNNING TIME: %(runtime)s (MM:SS:ms)" % pargs + "\n"
       
        startdt = None
        lastdt = None
        enddt = None
        wide = 75
        retstr += "\n\n"
        retstr += "SHOW IO".center(wide) + "\n"
        retstr += "-------".center(wide) + "\n"
        retstr += "\n"
        for dt in ks: # self.stephistory.keys():
            indent = sh[dt]["indent"]
            indentstr = "".join(["  " for i in range(0, indent)])
            
            mark = sh[dt]["mark"]
            if mark == "begin":
                elapsed = ""
            elif mark == "end":
                elapsed = "(" + str(dt - lastdt) + ") "
                
            pargs = {  "indent":indentstr,
                        "stepname":str(sh[dt]['stepname']),
                        "mark":str(sh[dt]['mark']),
                        "inputs":str(",".join(sh[dt]['inputs'])),
                        "outputs":str(",".join(sh[dt]['outputs']['main'])),
                        "time":str(dt),
                        "elapsed":elapsed,
                    }
            if startdt is None:
                retstr += ("%(inputs)s" % pargs).center(wide) + "\n"

            if (pargs["mark"] == "end"):
                retstr += " | ".center(wide) + "\n"
                retstr += "\|/".center(wide) + "\n"
                retstr += " ' ".center(wide) + "\n"
                
                line = ("%(stepname)s" % pargs).center(wide)
                line = "\x1b[1m" + line + "\x1b[22m" + "\n"
                retstr += line
                
            if len(sh[dt]["outputs"][MAINSTREAM]) != 0:
                retstr += " | ".center(wide) + "\n"
                retstr += "\|/".center(wide) + "\n"
                retstr += " ' ".center(wide) + "\n"
                retstr += ("%(outputs)s" % pargs).center(wide) + "\n"
                
                
            lastdt = dt
            if startdt is None:
                startdt = dt
        
        return retstr
        
    def report_output(self, inp, stream=None, load=True):
        """
        This function, along with ``get_inputs(..)`` allows a primitive to
        interact with the datastream in which it was invoked (or access
        other streams).

        :param inp: The inputs to report (add to the given or current stream).
            Input can be a string (filename), an AstroData instance, or a list of
            strings and/or AstroData instances.  Each individual dataset is
            wrapped in an AstroDataRecord and stored in the current stream.
        :type inp: str, AstroData instance, or list

        :param stream: If not specified the default ("main") stream is used.
            When specified the named stream is created if necessary.
        :type stream: str

        :param load: A boolean (default: True) which specifies whether string
            arguments (pathnames) should be loaded into AstroData instances
            or if it should be kept as a filename, unloaded.  This argument
            has no effect when "report"
            ``AstroData`` instances already in memory.
        :type load: <bool>

        """
        ##@@TODO: Read the new way code is done.
        if stream is None:
            stream = self._current_stream
        # this clause saves the output stream so we know when to 
        # the first report happens so we can clear the set at that time.
        if stream not in self._output_streams:
            self._output_streams.append(stream)
            self.outputs.update({stream:[]})
            
        # this clause makes sure there is a list in self.outputs
        if stream not in self.outputs:
            self.outputs.update({stream:[]})
            
        if type(inp) == str:
            self.outputs[stream].append(
                AstroDataRecord(inp, self.display_id, load=load)
            )
        elif isinstance(inp, AstroData):
            self.outputs[stream].append(AstroDataRecord(inp))
        elif type(inp) == list:
            for temp in inp:
                # This is a good way to check if IRAF failed.
                
                if type(temp) == tuple:
                    #@@CHECK: seems bad to assume a tuple means it is from 
                    #@@.....: a primitive that needs it's output checked!
                    if not os.path.exists(temp[0]):
                        raise "LAST PRIMITIVE FAILED: %s does not exist" % temp[0]
                    orecord = AstroDataRecord(temp[0], self.display_id, 
                                              parent=temp[1], load=load)
                elif isinstance(temp, AstroData):
                    orecord = AstroDataRecord(temp)
                elif isinstance(temp, AstroDataRecord):
                    orecord = temp
                elif type(temp) == str:
                    if not os.path.exists(temp):
                        raise "LAST PRIMITIVE FAILED."
                    orecord = AstroDataRecord(temp, self.display_id , load=load)
                else:
                    raise "RM292 type: " + str(type(temp))

                if stream not in self.outputs:
                    self.outputs.update({stream:[]})
                self.outputs[stream].append(orecord)
    
    def restore_fringe_index(self, filename):
        if os.path.exists(filename):
            self.fringes.stack_lists = pickle.load(open(filename, 'r'))
        else:
            pickle.dump({}, open(filename, 'w'))
        return

    def rm_cal(self, data, caltype):
        """
        Remove a calibration. This is used in command line argument (rmcal). 
        This may end up being used for some sort of TTL thing for cals in the 
        future.
        
        @param data: Images who desire their cals to be removed.
        @type data: str, list or AstroData instance.
        
        @param caltype: Calibration type (e.g. 'bias').
        @type caltype: str

        """
        datalist = gdpgutil.check_data_set(data)
        for dat in datalist:
            datid = IDFactory.generate_astro_data_id(data)
            key = (datid, caltype)
            if key in self.calibrations.keys():
                self.calibrations.pop(key)
            else:
                print "'%(tup)s', was not registered in the calibrations."
        return
    
    def rq_cal(self, caltype, inputs=None, source="all"):
        """
        Create calibration requests based on raw inputs.
        
        :param caltype: The type of calibration. For example, 'bias' and 'flat'.
        :type caltype: str
        
        :param inputs: The datasets for which to find calibrations, if not present
                        or ``None`` current "inputs" are used.
        :type inputs: list of AstroData instances
        
        :param source: Directs what calibration service to contact, for future
                        compatibility, currently only "all" is supported.
        :type source: <str>

        """
        if type(caltype) != str:
            raise RecipeError("caltype not string, type = " + str( type(caltype)))
        if inputs is None:
            # note: this was using original inputs!
            addToCmdQueue = self.cdl.get_cal_req(self.get_inputs_as_astrodata(),
                                                 caltype)
        else:
            addToCmdQueue = self.cdl.get_cal_req(inputs, caltype)
        
        for re in addToCmdQueue:
            re.calurl_dict = self["calurl_dict"]
            re.source = source
            self.add_rq(re)
        return
        
    def return_from_recipe(self):
        self._return_command = "return_from_recipe"
    
    def terminate_primitive(self):
        self._return_command = "terminate_primitive"
    
    def pop_return_command(self):
        rc = self._return_command
        self._return_command = None
        return rc
        
    def report_qametric(self, ad=None, name=None, metric_report=None, metadata=None):
        self._metricEvents.append_event(ad, name, metric_report, metadata=metadata)
        return

    def report_status(self, status):
        ad = status['adinput']
        current = { "current": status['current'],
                    "logfile": status['logfile'] }
        self._metricEvents.append_event(ad=ad, name="status", mdict=current, 
                                        msgtype='reduce_status')
        return

    def get_metric_list(self, clear=False):
        ml = self._metricEvents.get_list()
        if clear:
            self._metricEvents.clear_list()
        return ml
    
    def rq_display(self, display_id=None):
        """
        self, filename = None
        if None use self.inputs
        
        Create requests to display inputs.

        """
        ver = "1_0"
        displayObject = DisplayRequest()
        if display_id:
            Did = display_id
        else:
            Did = IDFactory.generate_display_id(self.inputs[0].filename, ver)
        displayObject.disID = Did
        displayObject.disList = self.inputs
        self.add_rq(displayObject)
        return
    
    def rq_iq(self, ad, e_m, e_s, f_m, f_s):
        iqReq = ImageQualityRequest(ad, e_m, e_s, f_m, f_s)
        self.add_rq(iqReq)
        return
    rq_iqput = rq_iq
        
    def rq_stack_get(self, purpose = ""):
        """
        The stackingID (see IDFactory module) is used to identify the list.
        The first input in the rc.inputs list is used as the reference image 
        to generate  
        the stackingID portion of the list identifier.
        
        The stackingID function in IDFactory is meant to produce identical
        stacking identifiers for different images which can/should be stacked 
        together, e.g. based
        on program ID and/or other details.  Again, see IDFactory for the
        particular algorithm in use.
        
        :note: a versioning system is latent within the code, and is added
            to the id to allow adaptation in the future if identifer construction
            methods change.

        :param purpose: The purpose is a string prepended to the stackingID
                        used to identify the list (see ``get_list(..)``).
        :type purpose: string

        """
        ver = "1_0"
        for orig in self.get_inputs_as_astrodata():
            Sid = purpose + IDFactory.generate_stackable_id(orig.ad, ver)
            stackUEv = GetStackableRequest()
            stackUEv.stk_id = Sid
            self.add_rq(stackUEv)
        return
                
    def rq_stack_update(self, purpose = None):
        """
        This function creates requests to update a stack list with the files
        in the current rc.inputs list.  Each will go in a stack based on its
        own stackingID (prepended with "purpose").
        
        :note: this function places a message on an outbound message queue
            which will not be sent until the next "yield", allowing the
            ReductionObject command clause to execute.

        :param purpose: The purpose argument is a string prefixed to the
            generated stackingID.  This allows two images which would
            produce identical stackingIDs to go in different lists,
            i.e. such as a fringe frame which might be prepended with
            "fringe" as the purpose.
        :type purpose: str
 
        """
        if purpose is None:
            purpose = ""
        ver = "1_0"
        inputs = self.get_inputs_as_astrodata()
        for inp in inputs:
            stackUEv = UpdateStackableRequest()
            Sid = purpose + IDFactory.generate_stackable_id(inp, ver)
            stackUEv.stk_id = Sid
            stackUEv.stk_list = inp.filename
            self.add_rq(stackUEv)
        return
    rq_stack_put = rq_stack_update
    
    def set_cache_file(self, key, filename):
        filename = os.path.abspath(filename)
        self.cache_files.update({key:filename})
        return
        
    def get_cache_file(self, key):
        if key in self.cache_files:
            return self.cache_files[key]
        else:
            return None
            
    def set_iraf_stderr(self, so):
        self.irafstderr = so
        return
    
    def set_iraf_stdout(self, so):
        self.irafstdout = so
        return
    
    def list_append(self, id, files, cachefile=None):
        """
        :param id: String identifies list to append the filenames.
        :type id: <str>

        :param files: List of filenames to add to the list.
        :type files: <list>

        :param cachefile: Filename to use to store the list.
        :type cachefile: <str>
        
        The caller is expected to supply ``cachefile``, which in principle
        a value of "None" could mean the "default cachefile" this is not
        supported by the ``adcc`` as of yet. The desired behavior is for
        reduce instances running in the same directory to cooperate, and those
        running in separate directories be kept separate, and this is 
        implemented by providing an argument for ``cachefile`` which is in a 
        generated subdirectory (hidden) based on the startup directory
        for the reduce process.  
        
        The adcc will negotiate all contention and race conditions regarding
        multiple applications manipulating a list simultaneously in separate
        process.
        """
        self.stackKeep.add(id, files, cachefile)
        return
    stack_append = list_append
        
    def list_inputs_as_str(self, id):
        """
        :param id: The identifier of the list to return as a comma separated 
        string w/ no whitespace
        :type id: str
        
        This is used to provide the list of names as a single string.
        """
        #pass back the stack files as strings
        stack = self.stackKeep.get(id)
        return ",".join(stack.filelist)
    stack_inputs_as_str = list_inputs_as_str

    def step_moment(self, stepname, mark):
        val = { "stepname"  : stepname,
                "indent"    : self.indent,
                "mark"      : mark,
                "inputs"    : [inp.filename for inp in self.inputs],
                "outputs"   : None,
                "processed" : False
                }
        return val
    
    def suffix_names(self, suffix, current_dir=True):
        '''
        
        '''
        newlist = []
        for nam in self.inputs:
            if current_dir == True:
                path = os.getcwd()
            else:
                path = os.path.dirname(nam.filename)
            
            fn = os.path.basename(nam.filename)
            finame, ext = os.path.splitext(fn)
            fn = finame + "_" + suffix + ext
            newpath = os.path.join(path, fn) 
            newlist.append(newpath)
        return newlist
        
    def switch_stream(self, switch_to = None):
        """
        :param switch_to: The string name of the stream to switch to. The 
            named stream must already exist.
        :type switch_to: str
        
        :note: This function is used by the infrastructure (in an application
            such as reduce and in the ReductionContext) to switch the stream
            being used. Reported output then goes to the specified stream.
        """
        if switch_to in self._output_streams:
            self._output_streams.remove(switch_to)

        if switch_to not in self.outputs:
            self.outputs.update({switch_to:[]})

        if switch_to not in self._stream_meta:
            self._stream_meta[switch_to] = {}

        switch_meta = self._stream_meta[switch_to]
        switch_meta["active_stream"] = True
        switch_meta["outputs_reported"] = False    
        
        self._current_stream = switch_to
        self._nonstandard_stream.append(switch_to)
        for ad in self.outputs[switch_to]:
            self.add_input(ad)

        return switch_to
        
    def restore_stream(self, from_stream=None):
        """
        Revert to the last stream prior to previous switch_stream(..) call.

        :param from_stream: This is the stream being reverted from. It does 
          not need to be passed in but can be used to ensure it is the same
          stream the rc thinks it is  popping off.
        :type from_stream: <str>

        """        
        if len(self._nonstandard_stream) > 0:
            prevstream = self._nonstandard_stream.pop()
            if from_stream and prevstream != from_stream:
                raise ReduceError("from_stream does not match last stream")
                
            if len(self._nonstandard_stream)>0:
                self._current_stream = self._nonstandard_stream[-1]
            else:
                self._current_stream = MAINSTREAM
        else:
            errmsg = "Cannot revert stream. No stream on stream list."
            errmsg += "The switch_stream(..) function not called."
            raise ReduceError(errmsg)

        return
                    
    def bad_call(self, arg=None):
        raise ReduceError("Forbidden: DO NOT USE ORIGINAL INPUTS")
        return
    original_inputs = property(bad_call, bad_call)
