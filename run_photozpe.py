#!/usr/local/bin/Enthought/epd-7.3-1-rh5-x86_64/bin/python
'''
Run our server!
To have this visible from outside, make sure you map port 80 to 5555,
 using nginx, for example.
Guide: http://sebschmidt.blogspot.com/2011/07/nginx-how-to-setup-tornado-and-apache.html
'''

from PhotoZPE import app
from sys import platform
if __name__ == '__main__':
    #if platform == 'darwin':
    app.run(debug = True)  # run with flask server
    '''
    else:
        from tornado.wsgi import WSGIContainer
        from tornado.httpserver import HTTPServer
        from tornado.ioloop import IOLoop
        
        http_server = HTTPServer(WSGIContainer(app))
        http_server.listen(5555)
        IOLoop.instance().start()
    '''
