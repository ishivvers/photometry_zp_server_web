'''
Run our server!
To have this visible from outside, make sure you map port 80 to 5555,
 using nginx, for example.
Guide: http://sebschmidt.blogspot.com/2011/07/nginx-how-to-setup-tornado-and-apache.html
'''
from PhotoZPE import app
from sys import platform

# the flask run command provides the best mechanism for debugging
#app.run(debug=True)

if __name__ == '__main__':
    from tornado.wsgi import WSGIContainer
    from tornado.httpserver import HTTPServer
    from tornado.ioloop import IOLoop

    
    ## main logging enabled for the nginx proxy, not the tornado server ##
    ## uncomment the below for debug-level logging in photo_zpe.log ##
    '''
    from tornado import options
    from time import strftime
    # log to text files
    import logging        
    log = logging.getLogger()
    log.setLevel(logging.DEBUG)
    log.addHandler( logging.handlers.TimedRotatingFileHandler('photo_zpe.log', when='W0') )
    logging.info("#"*10+strftime(" %Y-%m-%d %H:%M:%S ")+"#"*10+"\n###### Starting tornado web server ######")
    '''

    # start the server and let her run
    http_server = HTTPServer(WSGIContainer(app))
    http_server.listen(5555)
    IOLoop.instance().start()
