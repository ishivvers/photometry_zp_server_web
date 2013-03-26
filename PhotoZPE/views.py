########################################################################
from flask import redirect, request, session, url_for, render_template, Response, send_from_directory # helper functions
from PhotoZPE import app #the flask object itself, created by __init__.py
import numpy as np
from time import time, strftime
import json
import re, os
import pymongo as pm

from sys import platform
if 'linux' in platform.lower():
    import xmlrpclib
    # startup the xml connection to beast
    gs = xmlrpclib.ServerProxy('http://beast.berkeley.edu:5555', allow_none=True)
else:
    # import the local version
    from my_code import get_SEDs as gs


########################################################################
# The below can stay as global variables, since they don't change across threads
ALL_FILTERS = ['u','g','r','i','z','y','B','V','R','I','J','H','K']
# band: (central wavelength (AA), zeropoint (erg/s/cm^2/AA), catalog index)
FILTER_PARAMS =  {'u': (3551., 8.6387e-9, 0), 'g': (4686., 4.9607e-9, 1),
                  'r': (6165., 2.8660e-9, 2), 'i': (7481., 1.9464e-9, 3),
                  'z': (8931., 1.3657e-9, 4), 'y': (10091., 1.0696e-9, 5),
                  'B': (4400., 6.6000e-9, 6), 'V': (5490., 3.6100e-9, 7),
                  'R':(6500., 2.1900e-9, 7),  'I': (7885., 1.1900e-9, 9),
                  'J':(12350., 3.1353e-10, 10), 'H':(16620., 1.1121e-10, 11),
                  'K':(21590., 4.2909e-11, 12)}

# note, cannot put spectra folder in an app.config entry because 
#  the app hasn't yet been configured when this stuff is imported
spectra_folder = '/static/spectra/'
try:
    MODELS = np.load( app.root_path+spectra_folder+'all_models.npy' )
except:
    raise IOError('cannot find models file')
# convert the MODELS np array into a dictionary of arrays, so we can call by index (faster)
MODELS_DICT = {}
for model in MODELS[1:]:
    MODELS_DICT[model[0]] = model[1:]
del(MODELS) #just to free memory

try:
    SPEC_TYPES = np.loadtxt( app.root_path+spectra_folder+'pickles_types.txt', dtype=str )
except:
    raise IOError('cannot find spectral types file')

# SPEC_DICT contains all Pickles spectra in FLAM, as well as one array of lambda (Angstrom)
SPEC_DICT = {}
try:
    spec_files = os.listdir( app.root_path+spectra_folder )
    for fn in spec_files:
        try:
            n = int( re.findall('\d+', fn)[0] )
        except:
            continue
        data = np.load( app.root_path+spectra_folder+fn )
        SPEC_DICT[n] = data[1]
        if n == 1:
            SPEC_DICT['lam'] = data[0]
except:
    raise IOError('cannot find spectra files')
    

# initialize the database
DB = pm.MongoClient().PZserver

# the maximum number of sources to show on the spectrum page (does not affect catalog download)
MAX_DISP = 500
# the maximum field size to allow
MAX_FIELD = 3600.


# define a special class to handle xmlrpc weirdness
class special_dict(dict):
    '''
    A class that offers dictionary lookup through class attributes,
     so that interaction over the xmlrpc server is the same as through
     the local class.
    '''
    def __getattr__(self,name):
        return self[name]
    def __setattr__(self,name,value):
        self[name] = value

#######################################################################
@app.route('/upload', methods=['GET', 'POST'])
def show_upload():
    '''
    Uses upload.html (and base.html).
    When file uploaded, save it (date/time string), and then pass a pointer
    to the file name along to results! Done!
    '''
    if request.method == 'GET':
        # if this is a redirect, include any feedback (error messages)
        try:
            feedback = session['feedback']
            session['feedback'] = ''
        except:
            # if no feedback found, will simply fall to here
            feedback = ''
        # serve up the upload interface
        return render_template( "upload.html", feedback=feedback )
    else:
        # see which mode we're in
        mode = int(request.form['mode'])
        if mode == 1:
            # produce catalog for single field
            #  first test whether ra,dec, and fs pass muster:
            try:
                ra = round(parse_ra( request.form['RA'] ), 6)
                dec = round(parse_dec( request.form['DEC'] ), 6)
                fs = float( request.form['FS'] )
            except:
                return render_template( "upload.html", feedback="Could not interpret input - please enter valid coordinates "+\
                                                                "in decimal degrees or sexagesimal (HH:MM:SS.S, DD:MM:SS.S), "+\
                                                                "and make sure the requested field is {} arcseconds or smaller.".format(MAX_FIELD))
            if not (0. <= ra <= 360.) or not (-90. <= dec <= 90.) or not (0. <= fs <= MAX_FIELD):
                return render_template( "upload.html", feedback="Coordinates or field size beyond limits.")
            # if all's good, create the collection and populate it
            data = ra, dec, fs
            coll = create_collection()
            coll.insert( {"entry":"mode", "mode":mode} )
            coll.insert( {"entry":"search_field", "ra":ra, "dec":dec, "fs":fs} )
            # and take them straight to the results
            return redirect(url_for('show_results'))
        
        else:
            # parse the uploaded file
            source_file = request.files["source_file"]
            if source_file and allowed_file(source_file.filename):
                # try to parse with numpy
                try:
                    fn = app.root_path + '/tmp/' + str(np.random.randint(9999)) + '.txt'
                    source_file.save( fn )
                    data = np.loadtxt( fn ) #[:1000] #only accept first 1000 sources
                    # check to see whether we exceed the 1-degree limit
                    center, size = gs.find_field( data[:,:2].tolist() )
                    if max(size) > MAX_FIELD:
                        return render_template( "upload.html", feedback="Requested field exceeds size limit!")
                except:
                   return render_template( "upload.html", feedback="File upload failed! Make sure the file " + \
                                                            "is a properly-formatted ASCII file.")
            else:
                return render_template( "upload.html", feedback="File upload failed! Make sure the file " + \
                                                            "is a properly-formatted ASCII file that ends in .txt.")
            
            if mode == 2:
                # create matched catalog
                # if all's good, create the collection and populate it
                coll = create_collection()
                coll.insert( {"entry":"mode", "mode":mode} )
                for row in data:
                    coll["requested_sources"].insert( {"ra":row[0], "dec":row[1] })
                # have the user check that the file uploaded correctly
                return render_template( "upload.html", mode=mode, data=data[:5] )
            elif mode == 3:
                # produce matched catalog and zeropoint estimate
                # if all's good, create the collection and populate it
                band = request.form["band"]
                coll = create_collection()
                coll.insert( {"entry":"mode", "mode":mode} )
                coll.insert( {"entry":"passband", "passband":band} )
                for row in data:
                    coll["requested_sources"].insert( {"ra":row[0], "dec":row[1], "inst_mag":row[2] })
                # have the user check that the file uploaded correctly
                return render_template( "upload.html", mode=mode, data=data[:5], band=band )


## show_upload() helper functions
def create_collection():
    # create a new collection in the database, using the unix time and a random integer
    #  to ensure a datable (yet unique) collection.  Remove these with cronjob!'
    sid = str(time()).split('.')[0]+'_'+str(np.random.randint(9999))
    session['sid'] = sid #I use the flask.sessions interface to keep track of data across requests
    return DB[sid]


ALLOWED_EXTENSIONS = set(['txt','dat', 'cat'])
def allowed_file(filename):
    '''
    For now, simply return True.  Change if you want to constrain input filetypes.
    '''
    return True
    #return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


def parse_ra( inn ):
    '''
    Parse input RA string, either decimal degrees or sexagesimal HH:MM:SS.SS (or similar variants).
    Returns decimal degrees.
    '''
    # if simple float, assume decimal degrees
    try:
        ra = float(inn)
        return ra
    except:
        # try to parse with phmsdms:
        res = phmsdms(inn)
        ra = 15.*( res['vals'][0] + res['vals'][1]/60. + res['vals'][2]/3600. )
        return ra


def parse_dec( inn ):
    '''
    Parse input Dec string, either decimal degrees or sexagesimal DD:MM:SS.SS (or similar variants).
    Returns decimal degrees.
    '''
    # if simple float, assume decimal degrees
    try:
        dec = float(inn)
        return dec
    except:
        # try to parse with phmsdms:
        res = phmsdms(inn)
        dec = res['sign']*( res['vals'][0] + res['vals'][1]/60. + res['vals'][2]/3600. )
        return dec


def phmsdms(hmsdms):
    """
    +++ Pulled from python package 'angles' +++
    Parse a string containing a sexageismal number.
    
    This can handle several types of delimiters and will process
    reasonably valid strings. See examples.
    
    Parameters
    ----------
    hmsdms : str
        String containing a sexagesimal number.
        
    Returns
    -------
    d : dict
    
        parts : a 3 element list of floats
            The three parts of the sexagesimal number that were
            identified.
        vals : 3 element list of floats
            The numerical values of the three parts of the sexagesimal
            number.
        sign : int
            Sign of the sexagesimal number; 1 for positive and -1 for
            negative.
        units : {"degrees", "hours"}
            The units of the sexagesimal number. This is infered from
            the characters present in the string. If it a pure number
            then units is "degrees".
    """
    units = None
    sign = None
    # Floating point regex:
    # http://www.regular-expressions.info/floatingpoint.html
    #
    # pattern1: find a decimal number (int or float) and any
    # characters following it upto the next decimal number.  [^0-9\-+]*
    # => keep gathering elements until we get to a digit, a - or a
    # +. These three indicates the possible start of the next number.
    pattern1 = re.compile(r"([-+]?[0-9]*\.?[0-9]+[^0-9\-+]*)")
    # pattern2: find decimal number (int or float) in string.
    pattern2 = re.compile(r"([-+]?[0-9]*\.?[0-9]+)")
    hmsdms = hmsdms.lower()
    hdlist = pattern1.findall(hmsdms)
    parts = [None, None, None]
    
    def _fill_right_not_none():
        # Find the pos. where parts is not None. Next value must
        # be inserted to the right of this. If this is 2 then we have
        # already filled seconds part, raise exception. If this is 1
        # then fill 2. If this is 0 fill 1. If none of these then fill
        # 0.
        rp = reversed(parts)
        for i, j in enumerate(rp):
            if j is not None:
                break
        if  i == 0:
            # Seconds part already filled.
            raise ValueError("Invalid string.")
        elif i == 1:
            parts[2] = v
        elif i == 2:
            # Either parts[0] is None so fill it, or it is filled
            # and hence fill parts[1].
            if parts[0] is None:
                parts[0] = v
            else:
                parts[1] = v
                
    for valun in hdlist:
        try:
            # See if this is pure number.
            v = float(valun)
            # Sexagesimal part cannot be determined. So guess it by
            # seeing which all parts have already been identified.
            _fill_right_not_none()
        except ValueError:
            # Not a pure number. Infer sexagesimal part from the
            # suffix.
            if "hh" in valun or "h" in valun:
                m = pattern2.search(valun)
                parts[0] = float(valun[m.start():m.end()])
                units = "hours"
            if "dd" in valun or "d" in valun:
                m = pattern2.search(valun)
                parts[0] = float(valun[m.start():m.end()])
                units = "degrees"
            if "mm" in valun or "m" in valun:
                m = pattern2.search(valun)
                parts[1] = float(valun[m.start():m.end()])
            if "ss" in valun or "s" in valun:
                m = pattern2.search(valun)
                parts[2] = float(valun[m.start():m.end()])
            if "'" in valun:
                m = pattern2.search(valun)
                parts[1] = float(valun[m.start():m.end()])
            if '"' in valun:
                m = pattern2.search(valun)
                parts[2] = float(valun[m.start():m.end()])
            if ":" in valun:
                # Sexagesimal part cannot be determined. So guess it by
                # seeing which all parts have already been identified.
                v = valun.replace(":", "")
                v = float(v)
                _fill_right_not_none()
        if not units:
            units = "degrees"
            
    # Find sign. Only the first identified part can have a -ve sign.
    for i in parts:
        if i and i < 0.0:
            if sign is None:
                sign = -1
            else:
                raise ValueError("Only one number can be negative.")
                
    if sign is None:  # None of these are negative.
        sign = 1
        
    vals = [abs(i) if i is not None else 0.0 for i in parts]
    return dict(sign=sign, units=units, vals=vals, parts=parts)


@app.route('/', methods=['GET', 'POST'])
def home():
    # homepage simply points to upload
    return redirect(url_for('show_upload'))


@app.route('/favicon.ico')
def favicon():
    # quick redirect to show favicon
    return send_from_directory(app.root_path+'/static/img','favicon.ico')


@app.route('/robots.txt')
def robots():
    # quick redirect to serve robots.txt
    return send_from_directory(app.root_path+'/static/','robots.txt')


#######################################################################
@app.route('/results', methods=['GET'])
def show_results():
    '''
    The main results page, uses results*.html (and base.html).
    '''
    # two-tiered try-except clause to figure out whether the user submitted
    #  their own SID and then whether the known SID has a database entry tied to it
    try:
        # see whether a session ID was passed as an argument
        sid = request.args.get('sid')
        coll = DB[ sid ]
        session['sid'] = sid
        # The gymnastics above are neccessary because of strange errors
        #  when accessing the database with a redefined session['sid'].
        #  I don't fully understand what was wrong before, and don't fully
        #  understand why this works.
    except:
        try:
            coll = DB[ session['sid'] ]
        except:
            session['feedback'] = "Error: Could not find query results in database, please submit a new request."
            return redirect(url_for('show_upload'))
    try:
        # if database entry not found, shunt them back to upload
        mode = coll.find_one( {"entry":"mode"} )['mode']
    except:
        session['feedback'] = "Error: Could not find query results in database, please submit a new request."
        return redirect(url_for('show_upload'))
        
    # first, test to see whether we've already built a database, and simply display it
    if coll['data'].find_one() != None:
        curs = coll['data'].find()
        model_indices, coords = [], []
        for i in range(curs.count()):
            obj = curs.next()
            model_indices.append( obj["models"] )
            coords.append( obj["coords"] )
        coords = np.array(coords[:MAX_DISP])
        center, width, offsets = coords_to_offsets( coords )
        send_coords = []
        for i in range(len(coords)):
            send_coords.append( [coords[i][0], coords[i][1], offsets[i][0], offsets[i][1], model_indices[i], i] )
        if (mode == 1) or (mode == 2):
            return render_template( "results12.html", coords=json.dumps(send_coords),\
                                    field_center=json.dumps(center.tolist()), field_width=str(width))
        elif mode == 3:
            band = coll.find_one( {"entry":"passband"} )["passband"]
            zp_est = coll.find_one( {"entry":"zeropoint_estimate"})["zp"]
            zp_mad = coll.find_one( {"entry":"zeropoint_estimate"})["mad"]
            return render_template( "results3.html", coords=json.dumps(send_coords), field_center=json.dumps(center.tolist()),\
                        field_width=str(width), zp=round(zp_est,2), mad=round(zp_mad,2), band=band )
        
    # build the catalog and display it                     
    if mode == 1:
        search_field = coll.find_one( {"entry":"search_field"})
        ra = search_field['ra']
        dec = search_field['dec']
        fs = search_field['fs']
        try:
            cat = gs.catalog( (ra,dec), fs ) #square box of size fs
        except ValueError:
            session['feedback'] = "Error: No good sources found in requested field."
            return redirect(url_for('show_upload'))
        
        # handle differences in beast or local operation
        if not app.config['BEAST_MODE']:
            cat.SEDs = cat.SEDs.tolist()
            cat.full_errors = cat.full_errors.tolist()
            cat.coords = cat.coords.tolist()
        else:
            cat = special_dict(cat)
            
        # put into database
        for i in range(len(cat.SEDs)):
            coll['data'].insert( {"index":i, "sed":cat.SEDs[i], "errors":cat.full_errors[i],\
                                    "mode":cat.modes[i], "coords":cat.coords[i], "models":int(cat.models[i])} )
                                    
        # send to webpage
        coords = cat.coords[:MAX_DISP]
        center, width, offsets = coords_to_offsets( coords )
        send_coords = []
        for i in range(len(coords)):
            send_coords.append( [coords[i][0], coords[i][1], offsets[i][0], offsets[i][1], cat.models[i], i] )
        return render_template( "results12.html", coords=json.dumps(send_coords),\
                                    field_center=json.dumps(center.tolist()), field_width=str(width))
        
    elif mode == 2:
        requested_coords = []
        curs = coll['requested_sources'].find()
        for i in range(curs.count()):
            obj = curs.next()
            requested_coords.append( [obj['ra'], obj['dec']] )
        
        center, size = gs.find_field( requested_coords )
        try:
            cat = gs.catalog( center, max(size), requested_coords )
        except ValueError:
            session['feedback'] = "Error: No good sources found in requested field."
            return redirect(url_for('show_upload'))
        
        # handle differences in beast or local operation
        if not app.config['BEAST_MODE']:
            cat.SEDs = cat.SEDs.tolist()
            cat.full_errors = cat.full_errors.tolist()
            cat.coords = cat.coords.tolist()
        else:
            cat = special_dict(cat)
            
        # match requested coords to returned sources
        matches,tmp = gs.identify_matches( requested_coords, cat.coords )
        
        # put into database
        out_coords, out_model_indices = [],[]
        i = 0
        for j,match in enumerate(matches):
            if match >= 0:
                coll['data'].insert( {"index":i, "sed":cat.SEDs[match], "errors":cat.full_errors[match],\
                                        "mode":cat.modes[match], "coords":requested_coords[j], "models":int(cat.models[match])} )
                out_model_indices.append( cat.models[match] )
                out_coords.append( requested_coords[j] )
                i +=1
        
        # send to webpage
        out_coords = np.array(out_coords[:MAX_DISP])
        center, width, offsets = coords_to_offsets( out_coords )
        send_coords = []
        for i in range(len(out_coords)):
            send_coords.append( [out_coords[i][0], out_coords[i][1], offsets[i][0], offsets[i][1], out_model_indices[i], i] )
        return render_template( "results12.html", coords=json.dumps(send_coords),\
                                            field_center=json.dumps(center.tolist()), field_width=str(width))
        
    elif mode == 3:
        band = coll.find_one( {"entry":"passband"} )["passband"]
        requested_coords, inst_mags = [],[]
        curs = coll['requested_sources'].find()
        for i in range(curs.count()):
            obj = curs.next()
            requested_coords.append( [obj['ra'], obj['dec']] )
            inst_mags.append( obj['inst_mag'] )
        
        center, size = gs.find_field( requested_coords )
        try:
            cat = gs.catalog( center, max(size), requested_coords )
        except ValueError:
            session['feedback'] = "Error: No good sources found in requested field."
            return redirect(url_for('show_upload'))
        
        # handle differences in beast or local operation
        if not app.config['BEAST_MODE']:
            cat.SEDs = cat.SEDs.tolist()
            cat.full_errors = cat.full_errors.tolist()
            cat.coords = cat.coords.tolist()
        else:
            cat = special_dict(cat)
            
        # pull out only the band we care about
        iband = ALL_FILTERS.index(band)
        mod_mags = [ sss[ iband ] for sss in cat.SEDs ]
        
        # calculate zeropoint
        median_zp, mad_zp, matches, all_zps = gs.calc_zeropoint( requested_coords, cat.coords, inst_mags, mod_mags )
        coll.insert( {"entry":"zeropoint_estimate", "zp":median_zp, "mad":mad_zp} )
        for val in all_zps:
            coll['zeropoints'].insert( {"zp":val} )
        
        # put into db
        out_coords, out_model_indices = [],[]
        i = 0
        for j,match in enumerate(matches):
            if match >= 0:
                coll['data'].insert( {"index":i, "sed":cat.SEDs[match], "errors":cat.full_errors[match],\
                                        "mode":cat.modes[match], "coords":requested_coords[j], "models":int(cat.models[match])} )
                out_model_indices.append( cat.models[match] )
                out_coords.append( requested_coords[j] )
                i +=1
        
        # send to webpage
        out_coords = np.array(out_coords[:MAX_DISP])
        center, width, offsets = coords_to_offsets( out_coords )
        send_coords = []
        for i in range(len(out_coords)):
            send_coords.append( [out_coords[i][0], out_coords[i][1], offsets[i][0], offsets[i][1], out_model_indices[i], i] )
            
        return render_template( "results3.html", coords=json.dumps(send_coords), field_center=json.dumps(center.tolist()),\
                        field_width=str(width), zp=round(median_zp,2), mad=round(mad_zp,2), band=band )



def coords_to_offsets( coords ):
    '''
    Converts each set of coordinates to an angular offset from center.
    Used to query images from online DSS image server.
    Returns:
     field_center (ra,dec in decimal degrees)
	 field_width (in arcminutes) and
	 offsets_from_center (in arcminutes)
    '''
    coords = np.array(coords)
    edge = .5 # size to add to field edges in arcminutes
    
    r_c = np.deg2rad(min(coords[:,0]) + ( max(coords[:,0]) - min(coords[:,0]) )/2.)
    d_c = np.deg2rad(min(coords[:,1]) + ( max(coords[:,1]) - min(coords[:,1]) )/2.)
    
    # compute the offsets in radians
    r = np.deg2rad(coords[:,0])
    d = np.deg2rad(coords[:,1])
    offsets = np.empty_like( coords )
    offsets[:,0] = (np.cos(d)*np.sin(r-r_c))/(np.cos(d_c)*np.cos(d)*np.cos(r-r_c) +\
                                              np.sin(d_c)*np.sin(d))
    offsets[:,1] = (np.cos(d_c)*np.sin(d)-np.cos(d)*np.sin(d_c)*np.cos(r-r_c))/ \
                   (np.sin(d_c)*np.sin(d)+np.cos(d_c)*np.cos(d)*np.cos(r-r_c))
    
    # slightly rotate the offsets to mesh with the images from DSS server
    '''
    print abs(d_c) / np.pi
    if abs(d_c)>(np.pi/2):
        theta = (d_c/np.pi)*0.37
    else:
        theta = (d_c/np.pi)*0
    '''
    offsets = rotate_coords( offsets, 0. )
    
    # and finally present them in arcminutes
    offsets = np.rad2deg(offsets) * 60.
    
    # and now the width in arcminutes
    r_fw = max(offsets[:,0]) - min(offsets[:,0])
    d_fw = max(offsets[:,1]) - min(offsets[:,1])
    
    return np.rad2deg([r_c, d_c]).round(6), max([r_fw, d_fw])+2*edge, offsets


def rotate_coords( coords, theta ):
    '''
    Rotates coordinates about center by angle theta (all in radians).
    '''
    R = np.array( [[np.cos(theta), -1*np.sin(theta)], [np.sin(theta), np.cos(theta)]] )
    return np.dot(R, (coords).T).T

    
#######################################################################
@app.route('/info', methods=['GET'])
def show_info():
    '''
    A detailed description of the photometric estimates this produces.
    '''
    return render_template( "info.html" )


@app.route('/servespectrum', methods=['GET'])

def serve_spectrum():
    '''
    Loads spectrum loaded from numpy file and renormalizes it
    (given a spectrum id with url?spec=value&index=value), and then
    returns it in JSON format, as a set of objects with an x (Angstroms) and y (Flam).
    '''
    spec = int(request.args.get('spec',''))
    sed_index = int(request.args.get('index',''))
    # truncate the data below 2500AA    
    wl = SPEC_DICT['lam'][ SPEC_DICT['lam']>2500 ] #Angstroms
    
    # now match the model spectrum to the SED, using the Y-band
    #  to match (since Y will always be modeled)
    coll = DB[ session['sid'] ]
    curs = coll['data'].find_one( {"index":sed_index} )
    sed_mags = curs["sed"]
    sed_flam = mag2flam( sed_mags, ALL_FILTERS )
    # pull out the already-calculated model mag for yband
    y_model_flam = mag2flam( [MODELS_DICT[spec][5]], ['y'] ) #need to be array-like
    #D is the multiplier needed to make the model mesh with the fitted model
    D = sed_flam[5]/y_model_flam
    
    spec = D*SPEC_DICT[spec][ SPEC_DICT['lam']>2500 ]
    
    # push data into a json-able format: a list of dictionaries
    json_list = [{'x': wl[i], 'y': spec[i]} for i in range(len(wl))]
    return Response(json.dumps( json_list, indent=2 ), mimetype='application/json')


@app.route('/serveflams', methods=['GET'])
def serve_sed_flams():
    '''
    Loads magnitudes (obs and modeled) from database created by upload, returns FLAM.
    Needs database index (given as url?index=value&spec=value).
    '''
    sed_index = int(request.args.get('index',''))
    coll = DB[ session['sid'] ]
    curs = coll['data'].find_one( {"index":sed_index} )
    sed_mags = curs["sed"]
    mode = curs["mode"]
    if mode == 1:
        #USNOB+2MASS
        modeled = ['m']*6 + ['o','m']*2 + ['o']*3
    else:
        #SDSS+2MASS
        modeled = ['o']*5 + ['m']*5 + ['o']*3
    sed_flam = mag2flam( sed_mags, ALL_FILTERS )
    
    # push everything into json-able format
    json_list = [{'x': FILTER_PARAMS[ALL_FILTERS[i]][0], 'y': sed_flam[i], 'name': ALL_FILTERS[i],\
                    'modeled':modeled[i]} for i in range(len(sed_flam))]
    return Response(json.dumps( json_list, indent=2 ), mimetype='application/json')


@app.route('/servemags', methods=['GET'])
def serve_sed_mags():
    '''
    Loads magnitudes (obs and modeled) from database created by results page, returns mag.
     and the spectral type of the fit.
    Needs database index and spectrum id (given as url?index=value&spec=value).
    '''
    sed_index = int(request.args.get('index',''))   
    spec_index = int(request.args.get('spec',''))
    spec_type = SPEC_TYPES[spec_index].strip('IV') #remove dwarf/giant classifications
    
    coll = DB[ session['sid'] ]
    curs = coll['data'].find_one( {"index":sed_index} )
    sed_mags = curs["sed"]
    sed_errs = curs["errors"]
    mode = curs["mode"]
    
    if mode == 1:
        #USNOB+2MASS
        modeled = ['m']*6 + ['o','m']*2 + ['o']*3
    else:
        #SDSS+2MASS
        modeled = ['o']*5 + ['m']*5 + ['o']*3
    # push everything into json-able format
    json_list = [{'x': FILTER_PARAMS[ALL_FILTERS[i]][0], 'y': sed_mags[i], 'err': sed_errs[i], \
                 'modeled': modeled[i], 'name': ALL_FILTERS[i]} for i in range(len(sed_mags))]
    json_list += [{'name':'spec_type', 'value':spec_type}]
    return Response(json.dumps( json_list, indent=2 ), mimetype='application/json')


@app.route('/servecatalog', methods=['GET'])
def serve_full_catalog():
    '''
    Returns formatted & human-readable ASCII catalog of all sources.
    '''
    #mags, errs, mods, coords = DATA
    catalog_txt = \
    "# Catalog produced by the Photometric Estimate Server\n"+\
    "# http://classy.astro.berkeley.edu/ \n" +\
    "# Generated: {}\n".format(strftime("%H:%M %B %d, %Y")) +\
    '#  Mode = 0: -> B,V,R,I,y modeled from SDSS and 2-MASS\n' +\
    '#  Mode = 1: -> u,g,r,i,z,y,V,I modeled from USNOB-1 and 2-MASS\n' +\
    "# " + "RA".ljust(8) + "DEC".ljust(10) + "".join([f.ljust(8) for f in ALL_FILTERS]) +\
    "".join([(f+"_err").ljust(8) for f in ALL_FILTERS]) + "Mode\n"
    
    coll = DB[ session['sid'] ]
    curs = coll['data'].find()
    for i in range(curs.count()):
        obj = curs.next()
        catalog_txt += "".join([ s.ljust(10) for s in map(lambda x: "%.6f"%x, obj["coords"]) ])
        catalog_txt += "".join([ s.ljust(8) for s in map(lambda x: "%.3f"%x, obj["sed"]) ])
        catalog_txt += "".join([ s.ljust(8) for s in map(lambda x: "%.4f"%x, obj["errors"]) ])
        catalog_txt += str(obj["mode"])+"\n"
        
    response = Response(catalog_txt, mimetype='text/plain')
    response.headers['Content-Disposition'] = 'attachment; filename=catalog.txt'
    return response

    

@app.route('/servezp', methods=['GET'])
def serve_zeropoints():
    '''
    Loads and returns zeropoint estimates as inserted into the database by results page.
    Needs database key.
    '''
    coll = DB[ session['sid'] ]
    curs = coll['zeropoints'].find()
    json_list = []
    for i in range(curs.count()):
        obj = curs.next()
        json_list.append( {'zp': obj['zp']} )
    return Response(json.dumps( json_list, indent=2 ), mimetype='application/json')


@app.route('/api', methods=['GET','POST'])
def api_handler():
    '''
    Creates and returns (in JSON or ASCII format) a catalog of sources in a field,
     including both observed and modeled mags.  If fed an input file, attempts to 
     process it and provide cross-matched models and zeropoint estimate as well.
    If a 'GET' request, the url keys are: 'ra','dec','size', optional:'response'
    If a 'POST', a properly-formatted file must be uploaded, and a passband url key
     ('band') as well as a mode key ('mode') must be passed along.
    '''
    response_type = 'ascii' #what type of response to give. {'ascii', 'json'}
    try:
        typ = request.args.get('response')
        if 'json' in typ.lower():
            response_type = 'json'
    except:
        pass
    if request.method == 'GET':
        '''
        Produce & return a catalog in a field.  This is mode 1.
        '''
        try:
            ra = parse_ra( request.args.get('ra') )
            dec = parse_dec( request.args.get('dec') )
            size = float( request.args.get('size') )
        except:
            return Response( '{ success:false, message:"Request not formatted properly."}',
                            mimetype='application/json')
        if not (0. < ra < 360.) or not (-90. < dec < 90.) or not (0. < size < MAX_FIELD):
            return Response( '{ success:false, message:"Requested field either too large or entered incorrectly."}',
                            mimetype='application/json')
    
        # create a database entry, so we can view these results through the GUI too
        coll = create_collection()
        coll.insert( {"entry":"mode", "mode":1} )
        coll.insert( {"entry":"search_field", "ra":ra, "dec":dec, "fs":size} )
        try:
            cat = gs.catalog( (ra,dec), size )
        except ValueError:
            return Response( '{ success:false, message:"No suitable sources found in requested field."}',
                            mimetype='application/json')
            
        # handle differences in beast or local operation
        if not app.config['BEAST_MODE']:
            cat.SEDs = cat.SEDs.tolist()
            cat.full_errors = cat.full_errors.tolist()
            cat.coords = cat.coords.tolist()
        else:
            cat = special_dict(cat)
            
        # put the catalog entries both into the database and into a response
        json_list = [ {'query_ID':session['sid']} ]
        for i in range(len(cat.SEDs)):
            coll['data'].insert( {"index":i, "sed":cat.SEDs[i], "errors":cat.full_errors[i],\
                                    "mode":cat.modes[i], "coords":cat.coords[i], "models":int(cat.models[i])} )
            json_list.append({ 'ra':cat.coords[i][0], 'dec':cat.coords[i][1], 'mode':cat.modes[i], 'phot':cat.SEDs[i],\
                               'errors':cat.full_errors[i] })
        
        # return the catalog in either ascii or JSON
        if response_type == 'json':
            return Response(json.dumps( json_list, indent=2 ), mimetype='application/json')
        elif response_type == 'ascii':
            ascii_out = build_ascii( json_list[1:], "Query ID: {}".format(json_list[0]['query_ID']) )
            response = Response(ascii_out, mimetype='text/plain')
            response.headers['Content-Disposition'] = 'attachment; filename=catalog.txt'
            return response
    else:
        '''
        Produce and return a cross-matched catalog, perhaps with zeropoint estimate.  Modes 2, 3.
        '''
        try:
            mode = int( request.args.get('mode') )
            source_file = request.files["source_file"]
        except:
            return Response( '{ success:false, message:"Request not formatted properly."}',
                                mimetype='application/json')
        if source_file and allowed_file(source_file.filename):
            # try to parse with numpy
            try:
                fn = app.root_path + '/tmp/' + str(np.random.randint(9999)) + '.txt' 
                source_file.save( fn )
                data = np.loadtxt( fn ) #[:1000] #only accept first 1000 sources
                # check to see whether we exceed the 1-degree limit
                center, size = gs.find_field( data[:,:2].tolist() )
                if max(size) > MAX_FIELD:
                    return Response( '{ success:false, message:"Requested field size exceeds limit."}',
                                    mimetype='application/json')
            except:
                return Response( '{ success:false, message:"Uploaded file incorrectly formatted."}',
                                mimetype='application/json')
        else:
            return Response( '{ success:false, message:"Uploaded file needs to end in .txt."}',
                            mimetype='application/json')
        
        if mode == 2:
            # produce matched catalog only
            # if all's good, create the collection and populate it
            coll = create_collection()
            coll.insert( {"entry":"mode", "mode":mode} )
            for row in data:
                coll["requested_sources"].insert( {"ra":row[0], "dec":row[1] })
            requested_coords = data[:,:2].tolist()
            
            center, size = gs.find_field( requested_coords )
            try:
                cat = gs.catalog( center, max(size), requested_coords )
            except ValueError:
                return Response( '{ success:false, message:"No suitable sources found in requested field."}',
                                mimetype='application/json')
            
            # handle differences in beast or local operation
            if not app.config['BEAST_MODE']:
                cat.SEDs = cat.SEDs.tolist()
                cat.full_errors = cat.full_errors.tolist()
                cat.coords = cat.coords.tolist()
            else:
                cat = special_dict(cat)
            
            # match requested coords to returned sources
            matches,tmp = gs.identify_matches( requested_coords, cat.coords)
            
            # put into database
            json_list = [ {'query_ID':session['sid']} ]
            i = 0
            for j,match in enumerate(matches):
                if match >= 0:
                    coll['data'].insert( {"index":i, "sed":cat.SEDs[match], "errors":cat.full_errors[match],\
                                            "mode":cat.modes[match], "coords":requested_coords[j], "models":int(cat.models[match])} )
                    json_list.append({ 'ra':requested_coords[j][0], 'dec':requested_coords[j][1], 'mode':cat.modes[match],\
                                       'phot':cat.SEDs[match], 'errors':cat.full_errors[match] })
                    i +=1
            
            # return the catalog in either ascii or JSON
            if response_type == 'json':
                return Response(json.dumps( json_list, indent=2 ), mimetype='application/json')
            elif response_type == 'ascii':
                ascii_out = build_ascii( json_list[1:], "Query ID: {}".format(json_list[0]['query_ID']) )
                response = Response(ascii_out, mimetype='text/plain')
                response.headers['Content-Disposition'] = 'attachment; filename=catalog.txt'
                return response
        
        elif mode == 3:
            # produce matched catalog and zeropoint estimate
            try:
                band = request.args.get('band')
            except:
                return Response( '{ success:false, message:"Request not formatted properly."}',
                                mimetype='application/json')
            # if all's good, create the collection and populate it
            coll = create_collection()
            coll.insert( {"entry":"mode", "mode":mode} )
            coll.insert( {"entry":"passband", "passband":band} )
            for row in data:
                coll["requested_sources"].insert( {"ra":row[0], "dec":row[1], "inst_mag":row[2] })
            
            requested_coords = data[:,:2].tolist()
            inst_mags = data[:,2].tolist()
            
            center, size = gs.find_field( requested_coords )
            try:
                cat = gs.catalog( center, max(size), requested_coords )
            except ValueError:
                return Response( '{ success:false, message:"No suitable sources found in requested field."}',
                                mimetype='application/json')
            
            # handle differences in beast or local operation
            if not app.config['BEAST_MODE']:
                cat.SEDs = cat.SEDs.tolist()
                cat.full_errors = cat.full_errors.tolist()
                cat.coords = cat.coords.tolist()
            else:
                cat = special_dict(cat)
            
            # pull out only the band we care about
            iband = ALL_FILTERS.index(band)
            mod_mags = [ sss[ iband ] for sss in cat.SEDs ]
            
            # calculate zeropoint
            median_zp, mad_zp, matches, all_zps = gs.calc_zeropoint( requested_coords, cat.coords, inst_mags, mod_mags )
            coll.insert( {"entry":"zeropoint_estimate", "zp":median_zp, "mad":mad_zp} )
            for val in all_zps:
                coll['zeropoints'].insert( {"zp":val} )
            
            # put into database and into json or ascii format
            json_list = [ {'query_ID':session['sid'], 'median_zeropoint':median_zp, 'MAD_zeropoint':mad_zp} ]
            i = 0
            for j,match in enumerate(matches):
                if match >= 0:
                    coll['data'].insert( {"index":i, "sed":cat.SEDs[match], "errors":cat.full_errors[match],\
                                            "mode":cat.modes[match], "coords":requested_coords[j], "models":int(cat.models[match])} )
                    json_list.append({ 'ra':requested_coords[j][0], 'dec':requested_coords[j][1], 'mode':cat.modes[match],\
                                       'phot':cat.SEDs[match], 'errors':cat.full_errors[match] })
                    i +=1
            
            # return the catalog in either ascii or JSON
            if response_type == 'json':
                return Response(json.dumps( json_list, indent=2 ), mimetype='application/json')
            elif response_type == 'ascii':
                header = ["Query_ID: {}".format(json_list[0]['query_ID']), "Zeropoint: {}".format(round(json_list[0]['median_zeropoint'],2)),\
                          "M.A.D: {}".format(round(json_list[0]['MAD_zeropoint'],2))]
                ascii_out = build_ascii( json_list[1:], header )
                response = Response(ascii_out, mimetype='text/plain')
                response.headers['Content-Disposition'] = 'attachment; filename=catalog.txt'
                return response                    


def build_ascii( json_list, header_str=None ):
    '''
    Helper function for API.
    json_list: list of dictionaries to make into pretty ascii.
    header_str: either single-line string or list of strings to insert into the
      (commented-out) header of the ascii file.
    '''
    ascii_out = \
    "# Catalog produced by the Photometric Estimate Server\n"+\
    "# http://classy.astro.berkeley.edu/ \n" +\
    "# Generated: {}\n".format(strftime("%H:%M %B %d, %Y")) +\
    '#  Mode = 0: -> B,V,R,I,y modeled from SDSS and 2-MASS\n' +\
    '#  Mode = 1: -> u,g,r,i,z,y,V,I modeled from USNOB-1 and 2-MASS\n'
    if header_str != None:
        if type(header_str) == str:
            ascii_out += '# '+header_str+'\n'
        elif type(header_str) == list:
            for hs in header_str:
                ascii_out += '# '+hs+'\n'
    ascii_out += "# " + "RA".ljust(8) + "DEC".ljust(10) + "".join([f.ljust(8) for f in ALL_FILTERS]) +\
                 "".join([(f+"_err").ljust(8) for f in ALL_FILTERS]) + "Mode\n"
    for obj in json_list:
        ascii_out += "".join([ s.ljust(10) for s in map(lambda x: "%.6f"%x, [obj["ra"], obj["dec"]]) ])
        ascii_out += "".join([ s.ljust(8) for s in map(lambda x: "%.3f"%x, obj["phot"]) ])
        ascii_out += "".join([ s.ljust(8) for s in map(lambda x: "%.4f"%x, obj["errors"]) ])
        ascii_out += str(obj["mode"])+"\n"
    return ascii_out


#######################################################################
def mag2flam( magnitudes, bands ):
    '''
    Converts magnitudes to flam (erg/s/cm^2/A).
    '''
    flam = np.empty_like(magnitudes)
    for i,b in enumerate(bands):
        f0 = FILTER_PARAMS[b][1]
        flam[i] = f0*10**(-.4*magnitudes[i])
    return flam


def flam2mag( flams, bands ):
    '''
    Converts flam back to magnitudes.
    '''
    mags = np.empty_like(flams)
    for i,b in enumerate(bands):
        f0 = FILTER_PARAMS[b][1]
        mags[i] = -2.5*np.log10( flams[i]/f0 )
    return mags



