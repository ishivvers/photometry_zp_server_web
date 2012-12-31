
from flask import Flask
app = Flask(__name__) # the main flask object
from app import views # the site views, defined in views.py
