from flask import Flask
from sys import platform

app = Flask(__name__) # the main flask object
from PhotoZPE import views # the site views, defined in views.py

app.config['SPECTRA_FOLDER'] = app.root_path+'/static/spectra/'
app.config['SECRET_KEY'] = 'random!modnar'
