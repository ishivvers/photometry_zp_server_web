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
    '''
    
    return render_template( "results.html" )


@app.route('/contact', methods=['GET'])
def show_contact():
    '''
    A simple contacts page, with pointers back to my site.
    '''
    return 'contacts!'

@app.route('/served', methods=['GET'])
def serve_data():
    '''
    Test data-server to get D3 going.
    '''
    dat = np.load( '/Users/isaac/Working/code/photo_zp_server/web/app/static/pickles_uk_55.npy' )
    # push dat into a json-able format: a list of dictionaries
    json_list = [{'x': dat[0,i], 'y': dat[1,i]} for i in range(dat.shape[1])]
    return json.dumps( json_list )

