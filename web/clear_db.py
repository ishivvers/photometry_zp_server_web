#!/usr/bin/env python
'''
A python script to clear stored DB entries more than N days old.
Part of the PhotoZP Server project, this script should be run regularly
as a cronjob.
'''

import pymongo as pm
from time import time

age_cut = 1 #days
DB = pm.MongoClient().PZserver

now = int(time())
time_cut = now - age_cut*60*60*24

collections = DB.collection_names()
for coll in collections:
    try:
        root_time = int( coll.split('.')[0].split('_')[0] )
        if root_time < time_cut:
            DB.drop_collection( coll )
    except:
        pass
