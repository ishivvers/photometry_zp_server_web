#!/indirect/o/ishivvers/my_python/bin/python
'''
A python script to rotate DB logfiles, using the built-in
MongoDB logRotate command.
Part of the PhotoZP Server project, this script should be run regularly
as a cronjob.
'''

import pymongo as pm

DB = pm.MongoClient().admin
DB.command('logRotate')

