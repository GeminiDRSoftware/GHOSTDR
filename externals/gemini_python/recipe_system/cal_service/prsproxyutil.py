#
#                                                                     QAP Gemini
#
#                                                                prsproxyutil.py
#                                                                        07-2013
# ------------------------------------------------------------------------------
# $Id: prsproxyutil.py 5331 2015-09-14 23:48:57Z phirst $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5331 $'[11:-2]
__version_date__ = '$Date: 2015-09-14 13:48:57 -1000 (Mon, 14 Sep 2015) $'[7:-2]
# ------------------------------------------------------------------------------
from os.path import join, basename
from xml.dom import minidom
from pprint  import pformat

from astrodata.utils.Lookups import get_lookup_table
# ------------------------------------------------------------------------------
CALURL_DICT   = get_lookup_table("Gemini/calurl_dict","calurl_dict")
UPLOADPROCCAL = CALURL_DICT["UPLOADPROCCAL"]
UPLOADCOOKIE = CALURL_DICT["UPLOADCOOKIE"]
_CALMGR       = CALMGR      = CALURL_DICT["CALMGR"]
_LOCALCALMGR  = LOCALCALMGR = CALURL_DICT["LOCALCALMGR"]
# ------------------------------------------------------------------------------
CALTYPEDICT = { "arc": "arc",
                "bias": "bias",
                "dark": "dark",
                "flat": "flat",
                "processed_arc":   "processed_arc",
                "processed_bias":   "processed_bias",
                "processed_dark":   "processed_dark",
                "processed_flat":   "processed_flat",
                "processed_fringe": "processed_fringe"}
# -----------------------------------------------------------------------------
RESPONSESTR = """%%%%%%%%%%%%%%%% Request Data BEGIN:
%(sequence)s
%%%%%%%%%% Request Data END

########## Calibration Server Response BEGIN:
%(response)s
########## Calibration Server Response END

########## Nones Report (descriptors that returned None):
%(nones)s
########## Note: all descriptors shown above, scroll up.
        """
# -----------------------------------------------------------------------------

def upload_calibration(filename):
    """Uploads a calibration file to the FITS Store.

    parameters: <string>, the file to be uploaded.
    return:     <void>
    """
    import urllib2

    fn  = basename(filename)
    url = join(UPLOADPROCCAL, fn)
    postdata = open(filename).read()

    #"UPLOADCOOKIE": "qap_upload_processed_cal_ok",

    try:
        rq = urllib2.Request(url)
        rq.add_header('Content-Length', '%d' % len(postdata))
        rq.add_header('Content-Type', 'application/octet-stream')
        rq.add_header('Cookie', 'gemini_fits_upload_auth=%s' % UPLOADCOOKIE)
        u = urllib2.urlopen(rq, postdata)
        response = u.read()
    except urllib2.HTTPError, error:
        contents = error.read()
        raise
    return


def calibration_search(rq, fullResult=False):
    import urllib, urllib2
    import datetime

    calserv_msg = None
    print "\n@ppu074: calibration_search() ..."

    if "descriptors" in rq and "ut_datetime" in rq["descriptors"]:
        utc = rq["descriptors"]["ut_datetime"]
        pyutc = datetime.datetime.strptime(utc.value, "%Y%m%dT%H:%M:%S")
        print "@ppu079: OBS UT Date Time:", pyutc
        rq["descriptors"].update({"ut_datetime":pyutc} )
    
    # if rq has "calurl_dict" use it!!!
    if "calurl_dict" in rq :
        CALMGR = rq["calurl_dict"]["CALMGR"]
        LOCALCALMGR = rq["calurl_dict"]["LOCALCALMGR"]
    else:
        CALMGR = _CALMGR
        LOCALCALMGR = _LOCALCALMGR   
    
    if "source" not in rq:
        source = "central"
    else:
        source = rq["source"]
    
    token = ""             # used for GETs, we're using the post method
    rqurl = None

    print "@ppu098: CALMGR:     ", CALMGR
    print "@ppu099: LOCALCALMGR:", LOCALCALMGR

    if source == "central" or source == "all":
        rqurl = join(CALMGR, CALTYPEDICT[rq['caltype']])
        print "@ppu103: CENTRAL SEARCH: " + rqurl

    if source == 'local' or (rqurl == None and source == "all"):
        rqurl = LOCALCALMGR % { "httpport": 8777,
                                "caltype":  CALTYPEDICT[rq['caltype']],
                                } # "tokenstr":tokenstr}
        print "@ppu116: LOCAL SEARCH: " + rqurl

    rqurl = rqurl + "/%s" % rq["filename"]

    ### send request
    sequence = [("descriptors", rq["descriptors"]), ("types", rq["types"])]
    postdata = urllib.urlencode(sequence)
    response = "CALIBRATION_NOT_FOUND"

    try:
        print "@ppu119: REQUEST: ", rqurl
        calRQ = urllib2.Request(rqurl)
        if source == "local":
            u = urllib2.urlopen(calRQ, postdata)
        else:
            u = urllib2.urlopen(calRQ, postdata)
        response = u.read()
        print "@ppu126: RESPONSE: ", response
    except urllib2.HTTPError, error:
        calserv_msg = error.read()
        print "@ppu129: HTTPError- server returns:", error.read()
        import traceback
        traceback.print_exc()
        return (None, calserv_msg)

    #response = urllib.urlopen(rqurl).read()
    print "prs141:", response

    if fullResult:
        return response

    nones = []
    descripts = rq["descriptors"]
    for desc in descripts:
        if descripts[desc] == None:
            nones.append(desc)

    preerr = RESPONSESTR % { "sequence": pformat(sequence),
                             "response": response.strip(),
                             "nones"   : ", ".join(nones) \
                             if len(nones) > 0 else "No Nones Sent" }

    try:
        dom = minidom.parseString(response)
        calel = dom.getElementsByTagName("calibration")
        calurlel = dom.getElementsByTagName('url')[0].childNodes[0]
        calurlmd5 = dom.getElementsByTagName('md5')[0].childNodes[0]
    except IndexError:
        print "No url for calibration in response, calibration not found"
        return (None, preerr)
    except:
        return (None, preerr)
        
    #print "prs70:", calurlel.data
    
    #@@TODO: test only 
    print "@ppu165: ", repr(calurlel.data)
    return (calurlel.data, calurlmd5.data)