########################################################################
from flask import redirect, request, url_for, render_template # helper functions
from app import app #the flask object itself, created by __init__.py
import numpy as np
import json

########################################################################
# the website itself

# INDEX
@app.route('/upload', methods=['GET', 'POST'])
@app.route('/', methods=['GET'])
def show_upload():
    '''
    Uses upload.html (and base.html).
    '''
    if request.method == 'GET':
        # serve up the upload interface
        return render_template( "upload.html" )
    else:
        # just uploaded a file
        #need to manage uploads safely and properly, look into it here
        return render_template( "upload.html", feedback='successful upload!')
        

@app.route('/results', methods=['GET'])
def show_results():
    '''
    The main results page, uses results.html (and base.html).
    
    Feed template several lists of source info: 
     ra,dec
     spectrum_id
     [ mag estimates ]
    '''
    return render_template( "results.html", spec_ids=range(10,40) )


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
    '''
    f = '/Users/isaac/Working/code/photo_zp_server/web/app/static/spectra/pickles_uk_{}.npy'.format(request.args.get('spec',''))
    dat = np.load( f )
    # push dat into a json-able format: a list of dictionaries
    json_list = [{'x': dat[0,i], 'y': dat[1,i]} for i in range(dat.shape[1])]
    return json.dumps( json_list )

@app.route('/servemags', methods=['GET'])
def serve_mags():
    return 'serving mags'