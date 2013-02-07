########################################################################
from flask import redirect, request, session, url_for, render_template, Response # helper functions
from PhotoZPE import app #the flask object itself, created by __init__.py
import numpy as np
from time import time, strftime
import json
import pymongo as pm
from my_code import get_SEDs as gs #my library written to predict SEDs

'''
TO DO:
 - include link to simbad for each source
'''

########################################################################
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
    #MODELS = np.load( '/var/www/photozpe/my_code/all_models_P.npy','r' )
    MODELS = np.load( 'all_models_P.npy', 'r' )
except:
    raise IOError('cannot find models file')
# convert the MODELS np array into a dictionary of arrays, so we can call by index (faster)
MODELS_DICT = {}
for model in MODELS[1:]:
    MODELS_DICT[model[0]] = model[1:]
del(MODELS) #just to free memory

# initialize the database
DB = pm.MongoClient().PZserver

# the maximum number of sources to show on the spectrum page (does not affect catalog download)
max_disp = 300

#######################################################################
@app.route('/upload', methods=['GET', 'POST'])
@app.route('/', methods=['GET', 'POST'])
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
            # produce catalog for single field
            #  first test whether ra,dec, and fs pass muster:
            data = ra, dec, fs = map(float, [request.form['RA'], request.form['DEC'], request.form['FS']])
            if not (0. < ra < 360.) or not (-90. < dec < 90.) or not (0. < fs < 7200.):
                return render_template( "upload.html", feedback="Make sure you've entered valid parameters!")
            # if all's good, create the collection and populate it
            coll = create_collection()
            coll.insert( {"entry":"mode", "mode":mode} )
            coll.insert( {"entry":"search_field", "ra":ra, "dec":dec, "fs":fs} )
        elif mode == 2:
            # produce matched catalog
            source_file = request.files["source_file"]
            print source_file

            if source_file and allowed_file(source_file.filename):
                # try to parse with numpy
                try:
                    source_file.save( app.config['UPLOAD_FOLDER'] + 'tmp_source.txt' )
                    data = np.loadtxt( app.config['UPLOAD_FOLDER'] + 'tmp_source.txt' )[:1000] #only accept first 1000 sources
                except:
                    return render_template( "upload.html", feedback="File upload failed! Make sure the file " + \
                                                    "is a properly-formatted .txt file.")
            else:
                return render_template( "upload.html", feedback="File upload failed! Make sure the file " + \
                                                "is a properly-formatted .txt file.")
            # if all's good, create the collection and populate it
            coll = create_collection()
            coll.insert( {"entry":"mode", "mode":mode} )
            for row in data:
                coll["requested_sources"].insert( {"ra":row[0], "dec":row[1] })
        elif mode == 3:
            # produce matched catalog and zeropoint estimate
            source_file = request.files["source_file"]
            if source_file and allowed_file(source_file.filename):
                # try to parse with numpy
                try:
                    source_file.save( app.config['UPLOAD_FOLDER'] + 'tmp_source.txt' )
                    data = np.loadtxt( app.config['UPLOAD_FOLDER'] + 'tmp_source.txt' )[:1000]  #only accept first 1000 sources
                except:
                    return render_template( "upload.html", feedback="File upload failed! Make sure the file " + \
                                                    "is a properly-formatted .txt file.")
            else:
                return render_template( "upload.html", feedback="File upload failed! Make sure the file " + \
                                                "is a properly-formatted .txt file.")
            # if all's good, create the collection and populate it
            coll = create_collection()
            coll.insert( {"entry":"mode", "mode":mode} )
            coll.insert( {"entry":"passband", "passband":request.form["band"]} )
            for row in data:
                coll["requested_sources"].insert( {"ra":row[0], "dec":row[1], "inst_mag":row[2] })
            
        return render_template( "upload.html", mode=mode, data=data[:5] )
    
## show_upload() helper functions
def create_collection():
    # create a new collection in the database, using the unix time and a random integer
    #  to ensure a datable (yet unique) collection.  Remove these with cronjob!
    sid = str(time()).split('.')[0]+'_'+str(np.random.randint(9999))
    session['sid'] = sid #I use the flask.sessions interface to keep track of data across requests
    return DB[sid]
    
ALLOWED_EXTENSIONS = set(['txt','dat', 'cat'])
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


#######################################################################
@app.route('/results', methods=['GET'])
def show_results():
    '''
    The main results page, uses results.html (and base.html).
    '''
    coll = DB[ session['sid'] ]
    # case out the three different modes
    mode = coll.find_one( {"entry":"mode"} )['mode']

    # first, test to see whether we've already built a database, and simply display it
    if coll['data'].find_one() != None:
        curs = coll['data'].find()
        model_indices, coords = [], []
        for i in range(curs.count()):
            obj = curs.next()
            model_indices.append( obj["models"] )
            coords.append( obj["coords"] )
        if (mode == 1) or (mode == 2):
            return render_template( "results12.html", spec_ids=map(int, model_indices[:max_disp]), coords=coords[:max_disp] )
        elif mode == 3:
            band = coll.find_one( {"entry":"passband"} )["passband"]
            zp_est = coll.find_one( {"entry":"zeropoint_estimate"})["zp"]
            zp_mad = coll.find_one( {"entry":"zeropoint_estimate"})["mad"]
            return render_template( "results3.html", spec_ids=map(int, model_indices[:max_disp]), coords=coords[:max_disp],\
                                        zp=round(zp_est,2), mad=round(zp_mad,2), band=band )
            
    
    if mode == 1:
        search_field = coll.find_one( {"entry":"search_field"})
        ra = search_field['ra']
        dec = search_field['dec']
        fs = search_field['fs']
        coords, seds, models, errors, modes = gs.catalog( (ra,dec), (fs, fs), return_models=True ) #square box of size fs
        for i in range(len(seds)):
            coll['data'].insert( {"index":i, "sed":seds[i].tolist(), "errors":errors[i].tolist(),\
                                    "mode":modes[i], "coords":coords[i].tolist(), "models":int(models[i])} )
        return render_template( "results12.html", spec_ids=map(int, models[:max_disp]), coords=coords[:max_disp] )
        
    elif mode == 2:
        requested_coords = []
        curs = coll['requested_sources'].find()
        for i in range(curs.count()):
            obj = curs.next()
            requested_coords.append( [obj['ra'], obj['dec']] )
        requested_coords = np.array(requested_coords) #put into numpy array for sake of functions below
        
        center, size = gs.find_field( requested_coords )
        coords, seds, models, errors, modes = gs.catalog( center, size, object_coords=requested_coords, return_models=True )
        # match requested coords to returned sources
        matches = gs.identify_matches( requested_coords, coords)
        out_coords, out_model_indices = [],[]
        i = 0
        for match in matches:
            if match != None:
                coll['data'].insert( {"index":i, "sed":seds[match].tolist(), "errors":errors[match].tolist(),\
                                        "mode":modes[match], "coords":coords[match].tolist(), "models":int(models[match])} )
                out_model_indices.append( models[match] )
                out_coords.append( coords[match] )
                i +=1
        return render_template( "results12.html", spec_ids=map(int, out_model_indices[:max_disp]), coords=out_coords[:max_disp] )
    
    elif mode == 3:
        band = coll.find_one( {"entry":"passband"} )["passband"]
        requested_coords, inst_mags = [],[]
        curs = coll['requested_sources'].find()
        for i in range(curs.count()):
            obj = curs.next()
            requested_coords.append( [obj['ra'], obj['dec']] )
            inst_mags.append( obj['inst_mag'] )
        requested_coords = np.array(requested_coords) #put into numpy array for sake of functions below
        inst_mags = np.array(inst_mags)
        
        center, size = gs.find_field( requested_coords )
        coords, seds, models, errors, modes = gs.catalog( center, size, object_coords=requested_coords, return_models=True )
        
        # pull out only the band we care about
        mod_mags = np.array([ sss[ ALL_FILTERS==band ][0] for sss in seds ])
        
        # calculate zeropoint
        median_zp, mad_zp, matches, all_zps = gs.calc_zeropoint( requested_coords, coords, inst_mags, mod_mags, return_zps=True)
        coll.insert( {"entry":"zeropoint_estimate", "zp":median_zp, "mad":mad} )
        for val in all_zps:
            coll['zeropoints'].insert( {"zp":val} )
        # put into db
        out_coords, out_model_indices = [],[]
        i = 0
        for j,match in enumerate(matches):
            if match != None:
                coll['data'].insert( {"index":i, "sed":seds[match].tolist(), "errors":errors[match].tolist(),\
                                        "mode":modes[match], "coords":coords[match].tolist(), "models":int(models[match])} )
                out_model_indices.append( models[match] )
                out_coords.append( coords[match] )
                i +=1

        return render_template( "results3.html", spec_ids=map(int, out_model_indices[:max_disp]), coords=out_coords[:max_disp],\
                                    zp=round(median_zp,2), mad=round(mad_zp,2), band=band )
        


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
    f = app.config['SPECTRA_FOLDER']+'pickles_uk_{}.npy'.format(spec)
    dat = np.load( f )
    # truncate the data below 2500AA    
    wl = dat[0][ dat[0]>2500 ] #Angstroms
    
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
    spec = D*dat[1][ dat[0]>2500 ]
    
    # push data into a json-able format: a list of dictionaries
    json_list = [{'x': wl[i], 'y': spec[i]} for i in range(len(wl))]
    return Response(json.dumps( json_list ), mimetype='application/json')


@app.route('/serveflams', methods=['GET'])
def serve_sed_flams():
    '''
    Loads magnitudes (obs and modeled) from database created by upload, returns FLAM.
    Needs database index (given as url?index=value).
    '''
    sed_index = int(request.args.get('index',''))
    coll = DB[ session['sid'] ]
    curs = coll['data'].find_one( {"index":sed_index} )
    sed_mags = curs["sed"]
    mode = curs["mode"]
    if mode == 1:
        #USNOB+2MASS
        modeled = ['y']*6 + ['n']*5
    else:
        #SDSS+2MASS
        modeled = ['n']*5 + ['y']*3 + ['n']*3
    sed_flam = mag2flam( sed_mags, ALL_FILTERS )
    
    # push everything into json-able format
    json_list = [{'x': FILTER_PARAMS[ALL_FILTERS[i]][0], 'y': sed_flam[i], 'name': ALL_FILTERS[i],\
                    'modeled':modeled[i]} for i in range(len(sed_flam))]
    return Response(json.dumps( json_list ), mimetype='application/json')


@app.route('/servemags', methods=['GET'])
def serve_sed_mags():
    '''
    Loads magnitudes (obs and modeled) from database created by results page, returns FLAM.
    Needs database key (given as url?key=value).
    '''
    sed_index = int(request.args.get('index',''))    
    coll = DB[ session['sid'] ]
    curs = coll['data'].find_one( {"index":sed_index} )
    sed_mags = curs["sed"]
    sed_errs = curs["errors"]
    mode = curs["mode"]
    
    if mode == 1:
        #USNOB+2MASS
        modeled = ['y']*6 + ['n']*5
    else:
        #SDSS+2MASS
        modeled = ['n']*5 + ['y']*3 + ['n']*3
    # push everything into json-able format
    json_list = [{'x': FILTER_PARAMS[ALL_FILTERS[i]][0], 'y': sed_mags[i], 'err': sed_errs[i], \
                 'modeled': modeled[i], 'name': ALL_FILTERS[i]} for i in range(len(sed_mags))]
    return Response(json.dumps( json_list ), mimetype='application/json')


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
    "#\n#  Mode is the set of observations used to fit the model\n" +\
    "#   0=SDSS+2MASS, 1=USNOB+2MASS\n"+\
    "# " + "{}      {}       ".format("RA","DEC") + (' '*6).join(ALL_FILTERS) + " "*6 + '  '.join([val+"_err" for val in ALL_FILTERS]) + "  Mode\n"
    
    coll = DB[ session['sid'] ]
    curs = coll['data'].find()
    for i in range(curs.count()):
        obj = curs.next()
        catalog_txt += " ".join(map(lambda x: "%.6f"%x, obj["coords"]))+" "
        catalog_txt += " ".join(map(lambda x: "%.3f"%x, obj["sed"]))+" "
        catalog_txt += " ".join(map(lambda x: "%.4f"%x, obj["errors"]))+" "
        catalog_txt += str(obj["mode"])+"\n"
        
    response = Response(catalog_txt, mimetype='text/plain')
    response.headers['Content-Disposition'] = 'attachment; filename=catalog.txt'
    return response


@app.route('/servezp', methods=['GET'])
def serve_zeropoints():
    '''
    Loads and returns zeropoint estimates as inserted into the database by results page.
    Needs database key (given as url?key=value).
    '''
    coll = DB[ session['sid'] ]
    curs = coll['zeropoints'].find()
    json_list = []
    for i in range(curs.count()):
        obj = curs.next()
        json_list.append( {'zp': obj['zp']} )
    return Response(json.dumps( json_list ), mimetype='application/json')
    
    
    
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

