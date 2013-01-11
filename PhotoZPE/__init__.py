
from flask import Flask

app = Flask(__name__) # the main flask object
from PhotoZPE import views # the site views, defined in views.py

app.config['UPLOAD_FOLDER'] = '/var/www/photozpe/tmp/'
app.config['SECRET_KEY'] = 'random!modnar'
app.config['SPECTRA_FOLDER'] = '/var/www/photozpe/PhotoZPE/static/spectra/'
