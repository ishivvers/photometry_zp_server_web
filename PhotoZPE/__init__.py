
from flask import Flask
from sys import platform

app = Flask(__name__) # the main flask object
from PhotoZPE import views # the site views, defined in views.py

if 'linux' in platform.lower():
    app.config['UPLOAD_FOLDER'] = '/home/isaac/Working/code/photo_zp_server/tmp/'
    app.config['BEAST_MODE'] = True
else:
    app.config['UPLOAD_FOLDER'] = '/Users/isaac/Working/code/photo_zp_server/tmp/'
    app.config['BEAST_MODE'] = False
app.config['SPECTRA_FOLDER'] = app.root_path+'/static/spectra/'
app.config['SECRET_KEY'] = 'random!modnar'
