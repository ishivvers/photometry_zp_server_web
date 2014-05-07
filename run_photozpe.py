'''
Run our server!
To have this visible from outside, make sure you map port 80 to the listening app ports,
 using nginx, for example.
Guide: http://sebschmidt.blogspot.com/2011/07/nginx-how-to-setup-tornado-and-apache.html
'''
from PhotoZPE import app


if __name__ == '__main__':
    from sys import argv
    # will increment port by any integer passed as first cl argument
    port = 5000
    try:
        port += int(argv[1])
    except:
        pass
    
    # working tornado example, though it seems like flask works just fine
    '''
    from tornado.wsgi import WSGIContainer
    from tornado.httpserver import HTTPServer
    from tornado.ioloop import IOLoop
    from tornado import options
    from time import strftime
    # log to text files
    import logging        
    log = logging.getLogger()
    log.setLevel(logging.DEBUG)
    log.addHandler( logging.handlers.TimedRotatingFileHandler('photo_zpe.log', when='W0') )
    logging.info("#"*10+strftime(" %Y-%m-%d %H:%M:%S ")+"#"*10+"\n###### Starting tornado web server ######")
    
    http_server = HTTPServer(WSGIContainer(app))
    http_server.listen( port )
    IOLoop.instance().start()
    '''
    
    # start the simple flask server and let her see the larger world
    app.run(host='0.0.0.0', port=port ) #set debug=True for debugging interface
