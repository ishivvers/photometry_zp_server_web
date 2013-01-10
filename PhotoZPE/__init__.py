
from flask import Flask
app = Flask(__name__) # the main flask object

from PhotoZPE import views # the site views, defined in views.py


#app.config['UPLOAD_FOLDER'] = '/Users/isaac/Working/code/photo_zp_server/web/tmp/'
app.config['UPLOAD_FOLDER'] = '/home/isaac/Working/code/photo_zp_server/tmp/'
app.config['SECRET_KEY'] = 'random!modnar'
#app.config['SPECTRA_FOLDER'] = '/Users/isaac/Working/code/photo_zp_server/web/app/static/spectra/'
app.config['SPECTRA_FOLDER'] = '/home/isaac/Working/code/photo_zp_server/PhotoZPE/static/spectra/'
