#
#                                                                  gemini_python
#
#                                                        recipe_system.reduction
#                                                                   IDFactory.py
# ------------------------------------------------------------------------------
# $Id: IDFactory.py 5142 2015-02-17 21:39:45Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5142 $'[11:-2]
__version_date__ = '$Date: 2015-02-18 08:39:45 +1100 (Wed, 18 Feb 2015) $'[7:-2]
# ------------------------------------------------------------------------------
import re
import hashlib
import pyfits as pf
#------------------------------------------------------------------------------ 
version_index = {"stackID":"1_0", "recipeID":"1_0", "display_id":"1_0"}

def generate_md5_file(filename, version="1_0"):
    f = open(filename);
    block_size = 2**18 # 262 144 byte chunks
    md5 = hashlib.md5()
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()

    
def generate_fingerprint(dataset, version="1_0"):
    h = hashlib.md5()
    fullid = repr(dataset.types)+repr(dataset.descriptors)
    h.update(fullid)
    return h.hexdigest()
    
    
def generate_stackable_id(dataset, version="1_0"):
    """
    Generate an ID from which all similar stackable data will have in common.
    
    @param dataset: Input AstroData instance or fits filename.
    @type dataset: AstroData instances or str
    
    @param version: The version from which to run.
    @type version: string
    
    @return: A stackable id.
    @rtype: string  
    """
    from astrodata import AstroData
    if version != version_index['stackID']:
        try:
            # designed to call generateStackableID_
            idFunc = getattr(globals()['IDFactory'], 
                             'generateStackableID_' + version )
        except:
            raise "Version: '" + version + "' is either invalid or not supported." 
        
        return idFunc(inputf, version)
    
    """
        shaObj = hashlib.sha1()
        phu = pf.getheader( inputf[0], 0 )
        shaObj.update( phu['OBSID'] )
        shaObj.update( phu['OBJECT'] )
    """
    try:
        if isinstance(dataset, str):
            ad = AstroData(dataset)
            ID = version + str(ad.group_id())
        elif isinstance(dataset, AstroData):
            ID = version + str(dataset.group_id())
    except:
        print "Filename:", dataset.filename
        
    ID = make_id_safe_for_filename(ID)
    return ID
  
def generate_display_id(dataset, version):
    '''
    Generate an ID from which all similar stackable data will have in common.
    
    @param dataset: Input AstroData or fits filename
    @type dataset: list of AstroData instance
    
    @param version: The version from which to run.
    @type version: string
    
    @return: A display id.
    @rtype: string  
    '''
    from astrodata import AstroData
    if version != version_index['display_id']:
        try:
            # designed to call generateStackableID_
            idFunc = getattr( globals()['IDFactory'], 
                              'generateDisplayID_' + version )
        except:
            raise "Version: '" + version + "' is either invalid or not supported." 
        
        return idFunc(inputf, version)
    
    """
        shaObj = hashlib.sha1()
        phu = pf.getheader( inputf[0], 0 )
        shaObj.update( phu['OBSID'] )
        shaObj.update( phu['OBJECT'] )
    """
    
    if isinstance(dataset, str):
        phu = pf.getheader(dataset)
        ID = version + "_" + phu['OBSID'] + "_" + phu['OBJECT'] 
    elif isinstance(dataset, AstroData):
        ID = version + "_" + dataset.phu_get_key_value('OBSID') + "_" + \
             dataset.phu_get_key_value('OBJECT')
    return ID

def generate_astro_data_id(dataset, version="1_0"):
    """
    An id to be used to identify AstroData types. This is used for:
    
    1) Calibrations:
    
    Let's say a recipe performs
    
    getProcessedBias

    prepare
    biasCorrect
    
    Because of the prepare step, the calibration key determined at 
    getProcessedBias will not match biasCorrect because (N2009..., bias) will 
    not match (gN2009..., bias). By using an astro_id, you can avoid this issue 
    as you will have (DATALAB, bias). So, any steps inbetween getProcessedBias and
    biasCorrect will have no impact.
    
    2) Fringe:
    
    Fringe uses this as a FringeID, which is based off the first input of the list.
    
    @param dataset: Input AstroData instance or fits filename.
    @type dataset: AstroData instances or str
    
    @param version: The version from which to run.
    @type version: string
    
    @return: An astrodata id.
    @rtype: string  
    """
    from astrodata import AstroData
    if isinstance(dataset, str):
        ad = AstroData(dataset)
        return ad.data_label().as_pytype()
    elif isinstance(dataset, AstroData):
        return dataset.data_label().as_pytype()
    else:
        raise "BAD ARGUMENT TYPE: "+type(dataset)
    
def generate_fringe_list_id(dataset, version='1_0'):
    '''
    Generate an id from which all similar stackable data will have in common.
    
    @param dataset: Input AstroData instance or fits filename.
    @type dataset: AstroData instances or str
    
    @param version: The version from which to run.
    @type version: string
    
    @return: A stackable id.
    @rtype: string
    '''
    return generate_stackable_id(dataset, version)

def make_id_safe_for_filename(id):
    return re.sub(" ", "_", id)
