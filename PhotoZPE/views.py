########################################################################
from flask import redirect, request, session, url_for, render_template, Response, send_from_directory # helper functions
from PhotoZPE import app #the flask object itself, created by __init__.py
import numpy as np
from time import time, strftime
import json
import pymongo as pm
from my_code import get_SEDs as gs #my library written to predict SEDs

'''
TO DO:
 - put sources on image of sky intead of list
 - have nice error pages
'''

########################################################################
# The below can stay as global variables, since they don't change across threads
#  Central wavelength and zeropoint for all filters (zp in erg/s/cm^2/A )
FILTER_PARAMS =  gs.FILTER_PARAMS
ALL_FILTERS = np.array(gs.ALL_FILTERS)

try:
    MODELS = np.load( app.root_path+'/static/spectra/all_models.npy' )
except:
    raise IOError('cannot find models file')
# convert the MODELS np array into a dictionary of arrays, so we can call by index (faster)
MODELS_DICT = {}
for model in MODELS[1:]:
    MODELS_DICT[model[0]] = model[1:]
del(MODELS) #just to free memory

try:
    SPEC_TYPES = np.loadtxt( app.root_path+'/static/spectra/pickles_types.txt', dtype=str )
except:
    raise IOError('cannot find spectral types file')
    

# initialize the database
DB = pm.MongoClient().PZserver

# the maximum number of sources to show on the spectrum page (does not affect catalog download)
max_disp = 300

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


@app.route('/', methods=['GET', 'POST'])
def home():
    # homepage simply points to upload
    return redirect(url_for('show_upload'))


@app.route('/favicon.ico')
def favicon():
    # quick redirect to show favicon
    return send_from_directory(app.root_path+'/static/img','favicon.ico')


#######################################################################
@app.route('/results', methods=['GET'])
def show_results():
    '''
    The main results page, uses results.html (and base.html).
    '''
    # two-tiered try-except clause to figure out whether the user submitted
    #  their own SID and then whether the known SID has a database entry tied to it
    try:
        # see whether a session ID was passed as an argument
        sid = request.args.get('sid')
        coll = DB[ sid ]
        session['sid'] = sid
        # The gymnastics above are neccessary because of strange errors
        #  when accessing the database wit a redefined session['sid'].
        #  I don't fully understand what was wrong before, and don't fully
        #  understand why this works.
    except:
        try:
            coll = DB[ session['sid'] ]
        except:
            return redirect(url_for('show_upload'))
    try:
        # if database entry not found, shunt them back to upload
        mode = coll.find_one( {"entry":"mode"} )['mode']
    except:
        return redirect(url_for('show_upload'))
    
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
        cat = gs.catalog( (ra,dec), fs ) #square box of size fs
        for i in range(len(cat.SEDs)):
            coll['data'].insert( {"index":i, "sed":cat.SEDs[i].tolist(), "errors":cat.full_errors[i].tolist(),\
                                    "mode":cat.modes[i], "coords":cat.coords[i].tolist(), "models":int(cat.models[i])} )
        return render_template( "results12.html", spec_ids=map(int, cat.models[:max_disp]), coords=cat.coords[:max_disp] )
    
    elif mode == 2:
        requested_coords = []
        curs = coll['requested_sources'].find()
        for i in range(curs.count()):
            obj = curs.next()
            requested_coords.append( [obj['ra'], obj['dec']] )
        requested_coords = np.array(requested_coords) #put into numpy array for sake of functions below
        
        center, size = gs.find_field( requested_coords )
        cat = gs.catalog( center, max(size) )
        
        # match requested coords to returned sources
        matches = gs.identify_matches( requested_coords, cat.coords)
        out_coords, out_model_indices = [],[]
        i = 0
        for j,match in enumerate(matches):
            if match != None:
                coll['data'].insert( {"index":i, "sed":cat.SEDs[match].tolist(), "errors":cat.full_errors[match].tolist(),\
                                        "mode":cat.modes[match], "coords":requested_coords[j].tolist(), "models":int(cat.models[match])} )
                out_model_indices.append( cat.models[match] )
                out_coords.append( requested_coords[j] )
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
        cat = gs.catalog( center, max(size) )
        
        # pull out only the band we care about
        mod_mags = np.array([ sss[ ALL_FILTERS==band ] for sss in cat.SEDs ])
        
        # calculate zeropoint
        median_zp, mad_zp, matches, all_zps = gs.calc_zeropoint( requested_coords, cat.coords, inst_mags, mod_mags, return_zps=True)
        coll.insert( {"entry":"zeropoint_estimate", "zp":median_zp, "mad":mad_zp} )
        for val in all_zps:
            coll['zeropoints'].insert( {"zp":val} )
        # put into db
        out_coords, out_model_indices = [],[]
        i = 0
        for j,match in enumerate(matches):
            if match != None:
                coll['data'].insert( {"index":i, "sed":cat.SEDs[match].tolist(), "errors":cat.full_errors[match].tolist(),\
                                        "mode":cat.modes[match], "coords":requested_coords[j].tolist(), "models":int(cat.models[match])} )
                out_model_indices.append( cat.models[match] )
                out_coords.append( requested_coords[j] )
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
    Loads magnitudes (obs and modeled) from database created by results page, returns FLAM
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
    "#\n#  Mode is the set of observations used to fit the model\n" +\
    "#   0=SDSS+2MASS, 1=USNOB+2MASS\n"+\
    "# " + "{}      {}       ".format("RA","DEC") + (' '*6).join(ALL_FILTERS) +\
    " "*6 + '  '.join([val+"_err" for val in ALL_FILTERS]) + "  Mode\n"
    
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
        Produce & return a catalog in a field.
        '''
        try:
            ra = float( request.args.get('ra') )
            dec = float( request.args.get('dec') )
            size = float( request.args.get('size') )
        except:
            return Response( '{ success:false, message:"Request not formatted properly."}',
                            mimetype='application/json')
        if not (0. < ra < 360.) or not (-90. < dec < 90.) or not (0. < size < 7200.):
            return Response( '{ success:false, message:"Requested field either too large or entered incorrectly."}',
                            mimetype='application/json')
    
        # create a database entry, so we can view these results through the GUI too
        coll = create_collection()
        coll.insert( {"entry":"mode", "mode":1} )
        coll.insert( {"entry":"search_field", "ra":ra, "dec":dec, "fs":size} )
        cat = gs.catalog( (ra,dec), size )
    
        # put the catalog entries both into the database and into a response
        json_list = [ {'query_ID':session['sid']} ]
        for i in range(len(cat.SEDs)):
            coll['data'].insert( {"index":i, "sed":cat.SEDs[i].tolist(), "errors":cat.full_errors[i].tolist(),\
                                    "mode":cat.modes[i], "coords":cat.coords[i].tolist(), "models":int(cat.models[i])} )
            json_list.append({ 'ra':cat.coords[i][0], 'dec':cat.coords[i][1], 'mode':cat.modes[i], 'phot':cat.SEDs[i].tolist(),\
                               'errors':cat.full_errors[i].tolist() })
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
        Produce and return a cross-matched catalog, perhaps with zeropoint estimate.
        '''
        try:
            mode = int( request.args.get('mode') )
        except:
            return Response( '{ success:false, message:"Request not formatted properly."}',
                            mimetype='application/json')
        if mode == 2:
            # produce matched catalog only
            source_file = request.files["source_file"]
            if source_file and allowed_file(source_file.filename):
                # try to parse with numpy
                try:
                    source_file.save( app.config['UPLOAD_FOLDER'] + 'tmp_source.txt' )
                    data = np.loadtxt( app.config['UPLOAD_FOLDER'] + 'tmp_source.txt' )[:1000] #only accept first 1000 sources
                except:
                    return Response( '{ success:false, message:"Uploaded file incorrectly formatted."}',
                                    mimetype='application/json')
            else:
                return Response( '{ success:false, message:"Uploaded file incorrectly formatted."}',
                                mimetype='application/json')
            # if all's good, create the collection and populate it
            coll = create_collection()
            coll.insert( {"entry":"mode", "mode":mode} )
            for row in data:
                coll["requested_sources"].insert( {"ra":row[0], "dec":row[1] })
            
            requested_coords = data[:,:2] #put ra,dec into numpy array for sake of functions below
            
            center, size = gs.find_field( requested_coords )
            cat = gs.catalog( center, max(size) )
            # match requested coords to returned sources
            matches = gs.identify_matches( requested_coords, cat.coords)
            json_list = [ {'query_ID':session['sid']} ]
            i = 0
            for j,match in enumerate(matches):
                if match != None:
                    coll['data'].insert( {"index":i, "sed":cat.SEDs[match].tolist(), "errors":cat.full_errors[match].tolist(),\
                                            "mode":cat.modes[match], "coords":requested_coords[j].tolist(), "models":int(cat.models[match])} )
                    json_list.append({ 'ra':requested_coords[j][0], 'dec':requested_coords[j][1], 'mode':cat.modes[match],\
                                       'phot':cat.SEDs[match].tolist(), 'errors':cat.full_errors[match].tolist() })
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
            try:
                band = request.args.get('band')
            except:
                return Response( '{ success:false, message:"Request not formatted properly."}',
                                mimetype='application/json')
            # produce matched catalog and zeropoint estimate
            source_file = request.files["source_file"]
            if source_file and allowed_file(source_file.filename):
                # try to parse with numpy
                try:
                    source_file.save( app.config['UPLOAD_FOLDER'] + 'tmp_source.txt' )
                    data = np.loadtxt( app.config['UPLOAD_FOLDER'] + 'tmp_source.txt' )[:1000]  #only accept first 1000 sources
                except:
                    return Response( '{ success:false, message:"Uploaded file incorrectly formatted."}',
                                    mimetype='application/json')
            else:
                return Response( '{ success:false, message:"Uploaded file incorrectly formatted."}',
                                mimetype='application/json')
            # if all's good, create the collection and populate it
            coll = create_collection()
            coll.insert( {"entry":"mode", "mode":mode} )
            coll.insert( {"entry":"passband", "passband":band} )
            for row in data:
                coll["requested_sources"].insert( {"ra":row[0], "dec":row[1], "inst_mag":row[2] })
            
            requested_coords = data[:,:2] #put ra,dec into numpy array for sake of functions below
            inst_mags = data[:,2]
            
            center, size = gs.find_field( requested_coords )
            cat = gs.catalog( center, max(size) )
            
            # pull out only the band we care about
            mod_mags = np.array([ sss[ ALL_FILTERS==band ] for sss in cat.SEDs ])
            
            # calculate zeropoint
            median_zp, mad_zp, matches, all_zps = gs.calc_zeropoint( requested_coords, cat.coords, inst_mags, mod_mags, return_zps=True)
            coll.insert( {"entry":"zeropoint_estimate", "zp":median_zp, "mad":mad_zp} )
            for val in all_zps:
                coll['zeropoints'].insert( {"zp":val} )
            # put into database and into json or ascii format
            json_list = [ {'query_ID':session['sid'], 'median_zeropoint':median_zp, 'MAD_zeropoint':mad_zp} ]
            i = 0
            for j,match in enumerate(matches):
                if match != None:
                    coll['data'].insert( {"index":i, "sed":cat.SEDs[match].tolist(), "errors":cat.full_errors[match].tolist(),\
                                            "mode":cat.modes[match], "coords":requested_coords[j].tolist(), "models":int(cat.models[match])} )
                    json_list.append({ 'ra':requested_coords[j][0], 'dec':requested_coords[j][1], 'mode':cat.modes[match],\
                                       'phot':cat.SEDs[match].tolist(), 'errors':cat.full_errors[match].tolist() })
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
    "#\n#  Mode is the set of observations used to fit the model\n" +\
    "#   0=SDSS+2MASS, 1=USNOB+2MASS\n"
    if header_str != None:
        if type(header_str) == str:
            ascii_out += '# '+header_str+'\n'
        elif type(header_str) == list:
            for hs in header_str:
                ascii_out += '# '+hs+'\n'
    ascii_out += "# " + "{}       {}       ".format("RA","DEC") + (' '*6).join(ALL_FILTERS) + " "*6 +\
    '  '.join([val+"_err" for val in ALL_FILTERS]) + "  Mode\n"
    for obj in json_list:
        ascii_out += " ".join(map(lambda x: "%.6f"%x, [obj["ra"], obj["dec"]]))+" "
        ascii_out += " ".join(map(lambda x: "%.3f"%x, obj["phot"]))+" "
        ascii_out += " ".join(map(lambda x: "%.4f"%x, obj["errors"]))+" "
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



