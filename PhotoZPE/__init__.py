
from flask import Flask
from sys import platform

app = Flask(__name__) # the main flask object
from PhotoZPE import views # the site views, defined in views.py

if platform == 'Linux':
    app.config['UPLOAD_FOLDER'] = '/var/www/photozpe/tmp/'
    app.config['SPECTRA_FOLDER'] = '/var/www/photozpe/PhotoZPE/static/spectra/'
else:
    app.config['UPLOAD_FOLDER'] = '/Users/isaac/Working/code/photo_zp_server/tmp/'
    app.config['SPECTRA_FOLDER'] = '/Users/isaac/Working/code/photo_zp_server/PhotoZPE/static/spectra/'
app.config['SECRET_KEY'] = 'random!modnar'