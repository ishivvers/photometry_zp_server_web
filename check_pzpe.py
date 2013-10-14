#!/usr/bin/env python
'''
A quick script that checks for the photozpe page to be up and working,
 and emails me if it is not.
Should be regularly run by a cronjob.
'''

import smtplib          # used for email
from email.mime.text import MIMEText
import urllib2
from subprocess import Popen,PIPE

toAddress="ishivvers@berkeley.edu"
fromAddress="pzpe_checker@astro.berkeley.edu"
testPage='http://classy.astro.berkeley.edu/'

def send_msg(toAddress, fromAddress, subject, message):
    '''Send a message'''
    # Set up message text
    msg=MIMEText(message)
    msg["From"]=fromAddress
    msg["To"]=toAddress
    msg["Subject"]=subject

    # Send email
    server = smtplib.SMTP('chandra.berkeley.edu')
    server.set_debuglevel(0)
    server.sendmail(fromAddress,toAddress,msg.as_string())
    server.quit()

try:
    # first test that the page is up
    page = urllib2.urlopen(testPage, timeout=30).read()
    start = page.find('<title>')
    assert (start>0)
    end = page.find('</title>')
    assert (end>0)
    assert (page[start:end].strip('<title>')=='PhotoZPE')

    # also test that the whole chain works, by making a quick API call
    res = Popen("hostname",shell=True,stdout=PIPE,stderr=PIPE)
    hostname,err = res.communicate()
    if hostname.strip() == 'UCBerk':
        targetfile = '/home/isaac/Downloads/tmp.txt'
    elif hostname.strip() == 'lupus.berkeley.edu':
        targetfile = '/o/ishivvers/scratch/tmp.txt'
    else:
        raise Exception('What machine are we on, here?')
    res = Popen("wget -O "+targetfile+" 'http://classy.astro.berkeley.edu/api?method=1&ra=200.&dec=20.&size=150.&response=ascii'",
                shell=True,stdout=PIPE,stderr=PIPE)
    out,err = res.communicate() #wait for the download to finish
    lines = open(targetfile,'r').readlines()
    Popen("rm "+targetfile, shell=True)
    assert (lines[0]=='# Catalog produced by the Photometric Estimate Server\n')

except Exception as E:
    send_msg(toAddress, fromAddress, '[PhotZPE] Page Down Error',
        'I report that the PhotoZPE page appears to be down!\n\n'+str(E)+\
        '\n\nIn solidarity,\nthe classiest robot')
