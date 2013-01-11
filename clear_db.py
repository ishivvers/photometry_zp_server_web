#!/usr/bin/env python
'''
A python script to clear stored DB entries more than N days old, and
delete old uploaded files.
Part of the PhotoZP Server project, this script should be run regularly
as a cronjob.
'''

import pymongo as pm
from time import time
from os import listdir, system

age_cut = 1 #days
DB = pm.MongoClient().PZserver
tmp_folder = '/var/www/photozpe/tmp/'

now = int(time())
time_cut = now - age_cut*60*60*24

# clear db
collections = DB.collection_names()
for coll in collections:
    try:
        root_time = int( coll.split('.')[0].split('_')[0] )
        if root_time < time_cut:
            DB.drop_collection( coll )
    except:
        pass

# clear tmp folder
txt_files = listdir(tmp_folder)
for f in txt_files:
    try:
        root_time = int( f.split('.')[0].split('_')[0] )
        if root_time < time_cut:
            system( 'rm {}{}'.format(tmp_folder, f) )
    except:
        pass