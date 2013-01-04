########################################################################
from flask import redirect, request, session, url_for, render_template # helper functions
from app import app #the flask object itself, created by __init__.py
import numpy as np
from time import time
import json
from my_code import get_SEDs as gs #my library written to predict SEDs

'''
TO DO:
 - move globals into session, so that apache can thread this, and have a 
   persistant database that carries around the data for each request.
'''

########################################################################
# global variables -- should port these into cookie/session keys later

# ID for this session -- all saved files et cetera will use this identifier.
#  Include a cronjob to delete obsolete/old files -- maybe delete all more than a day old?
CURRENT_ID = str(time()).split('.')[0]
CURRENT_MODE = None
DATA = None #a global container for the data passed along to results page

# The below can stay as global variables, since they don't change across threads
#  Central wavelength and zeropoint for all filters (zp in erg/s/cm^2/A )
FILTER_PARAMS =  {'u': (3551., 8.5864e-9), 'g': (4686., 4.8918e-9),
                  'r': (6165., 2.8473e-9), 'i': (7481., 1.9367e-9),
                  'z': (8931., 1.3564e-9), 'y': (10091., 1.0696e-9),
                  'B': (4400., 6.6000e-9),  'R':(6500., 2.1900e-9),
                  'J':(12350., 3.1353e-10), 'H':(16620., 1.1121e-10),
                  'K':(21590., 4.2909e-11)}
ALL_FILTERS = np.array(['u','g','r','i','z','y','B','R','J','H','K'])

try:
    MODELS = np.load( open('all_models_P.npy','r') )
except:
    raise IOError('cannot find models file')
# convert the MODELS np array into a dictionary of arrays, so we can call by index (faster)
MODELS_DICT = {}
for model in MODELS[1:]:
    MODELS_DICT[model[0]] = model[1:]
del(MODELS) #just to free memory
    

#######################################################################
@app.route('/upload', methods=['GET', 'POST'])
def show_upload():
    '''
    Uses upload.html (and base.html).
    When file uploaded, save it (date/time string), and then pass a pointer
    to the file name along to results! Done!
    '''
    if request.method == 'GET':
        # serve up the upload interface
        return render_template( "upload.html" )
    else:
        # see which mode we're in
        mode = int(request.form['mode'])
        if mode == 1:
            globals()['CURRENT_MODE'] = 1
            data = map(float, [request.form['RA'], request.form['DEC'], request.form['FS']])
        else:
            globals()['CURRENT_MODE'] = mode
            source_file = request.files["source_file"]
            if source_file and allowed_file(source_file.filename):
                # try to parse with numpy
                try:
                    source_file.save( app.config['UPLOAD_FOLDER'] + '{}_source.txt'.format(CURRENT_ID) )
                    data = np.loadtxt( app.config['UPLOAD_FOLDER'] + '{}_source.txt'.format(CURRENT_ID) )
                except:
                    return render_template( "upload.html", feedback="File upload failed! Make sure the file " + \
                                                    "is a .txt file readable with numpy.loadtxt()!")
            else:
                return render_template( "upload.html", feedback="File upload failed! Make sure the file " + \
                                                "is a .txt file readable with numpy.loadtxt()!")
        globals()['DATA'] = data
        return render_template( "upload.html", mode=CURRENT_MODE, data=data[:5] )
    
## show_upload() helper functions
ALLOWED_EXTENSIONS = set(['txt','dat'])
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


#######################################################################
@app.route('/results', methods=['GET'])
def show_results():
    '''
    The main results page, uses results.html (and base.html).
    
    Feed template several lists of source info: 
     ra,dec
     spectrum_id
     multiplier/offset to spectrum
     [ mag estimates ]
     
    '''
    # case out the three different modes
    if CURRENT_MODE == 1:
        ra,dec,fs = DATA
        coords, seds, models = gs.catalog( (ra,dec), (fs, fs), return_models=True ) #square box of size fs
        model_indices, model_offsets = zip(*models)
        globals()['DATA'] = (seds, model_offsets)
        return render_template( "results.html", spec_ids=map(int, model_indices), coords=coords )





#######################################################################
@app.route('/contact', methods=['GET'])
def show_contact():
    '''
    A simple contacts page, with pointers back to my site.
    '''
    return 'contacts!'

@app.route('/servespectrum', methods=['GET'])
def serve_spectrum():
    '''
    Loads spectrum loaded from numpy file (given a spectrum id with url?spec=value)
    and returns it in JSON format, as a set of objects with an x (Angstroms) and y (Flam).
    
    Need to modify such that this also accepts a multiplier/offset to equilibrate
     the spectrum and the mags.
    '''
    spec = int(request.args.get('spec',''))
    sed_index = int(request.args.get('index',''))
    f = '/Users/isaac/Working/code/photo_zp_server/web/app/static/spectra/pickles_uk_{}.npy'.format(spec)
    dat = np.load( f )
    # truncate the data below 2500AA    
    wl = dat[0][ dat[0]>2500 ] #Angstroms
    
    # now match the model spectrum to the SED, using the Y-band
    #  to match (since Y will always be modeled)
    sed_mags = DATA[0][sed_index]
    sed_flam = mag2flam( sed_mags, ALL_FILTERS )
    # pull out the already-calculated model mag for yband
    y_model_flam = mag2flam( [MODELS_DICT[spec][5]], ['y'] ) #need to be array-like
    #D is the multiplier needed to make the model mesh with the fitted model
    D = sed_flam[5]/y_model_flam
    spec = D*dat[1][ dat[0]>2500 ]
    
    # push data into a json-able format: a list of dictionaries
    json_list = [{'x': wl[i], 'y': spec[i]} for i in range(len(wl))]
    return json.dumps( json_list )


@app.route('/servemags', methods=['GET'])
def serve_mags():
    '''
    Loads magnitudes (obs and modeled) from database created by upload.
    Needs database key (given as url?key=value).
    
    Include, in database, both mags and values in FLAM (erg/s/W^2/A) for plotting.
    '''
    sed_index = int(request.args.get('index',''))
    sed_mags = DATA[0][sed_index]
    sed_flam = mag2flam( sed_mags, ALL_FILTERS )
    
    
    # push wavelengths (A) and flam into json-able format
    json_list = [{'x': FILTER_PARAMS[ALL_FILTERS[i]][0], 'y': sed_flam[i], 'name': ALL_FILTERS[i]} for i in range(len(sed_flam))]
    return json.dumps( json_list )






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

