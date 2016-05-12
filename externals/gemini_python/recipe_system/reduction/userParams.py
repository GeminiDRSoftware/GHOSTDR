#
#                                                                  gemini_python
#
#
#                                                                  userParams.py
# ------------------------------------------------------------------------------
# $Id: userParams.py 5142 2015-02-17 21:39:45Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5142 $'[11:-3]
__version_date__ = '$Date: 2015-02-18 08:39:45 +1100 (Wed, 18 Feb 2015) $'[7:-3]
# ------------------------------------------------------------------------------
""" Module provides user parameter classes for structured user parameters as
may be passed through the reduce command line interface.

Classes:
    UserParam  -- a single parameter
    UserParams -- a set of UserParam instances
"""
from .recipeManager import RecipeError
# ------------------------------------------------------------------------------
class UserParam(object):

    def __init__(self, astrotype, primname, param, value):
        self.astrotype = astrotype
        self.primname  = primname
        self.param     = param
        self.value     = value
    
    def __repr__(self):
        ret = "UserParam: adtype=%s primname=%s %s=%s" % (repr(self.astrotype),
                                                          repr(self.primname),
                                                          repr(self.param),
                                                          repr(self.value),
                                                      )
        return ret

# ------------------------------------------------------------------------------
class UserParams(object):

    def __init__(self):
        self.user_param_dict = None

    def __repr__(self):
        ret = "UserParams: "
        ret += repr(self.user_param_dict)
        return ret
        
    def __contains__(self, arg):
        if self.user_param_dict:
            return self.user_param_dict.__contains__(arg0)
    
    def is_empty(self):
        if self.user_param_dict:
            return False
        else:
            return True
            
    def get_user_param(self, astrotype, primname):
        if self.user_param_dict == None:
            return None
        if astrotype not in self.user_param_dict:
            return None
        if primname not in self.user_param_dict[astrotype]:
            return None
        return self.user_param_dict[astrotype][primname]
        
    def add_user_param(self, userparam):
        up = userparam
        if userparam == None:
            return
        if self.user_param_dict == None:
            self.user_param_dict = {}
            
        if up.astrotype not in self.user_param_dict:
            self.user_param_dict.update({up.astrotype: {}})
        
        if up.primname not in self.user_param_dict[up.astrotype]:
            self.user_param_dict[up.astrotype].update({up.primname: {}})
            
        if up.param in self.user_param_dict[up.astrotype][up.primname]:
            raise RecipeError("Parameter (%s.%s%s) already set by user" % \
            (up.astrotype, up.primname, up.param))
        else:
            self.user_param_dict[up.astrotype][up.primname].update({up.param:up.value})
        return
