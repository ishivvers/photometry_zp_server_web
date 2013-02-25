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
    if platform == 'darwin':
        app.run(debug = True)  # run with flask server
    else:
        from tornado.wsgi import WSGIContainer
        from tornado.httpserver import HTTPServer
        from tornado.ioloop import IOLoop
        from tornado import options
        
        # text-only logging
        #options.options['log_file_prefix'].set('/var/www/logs/photo_zpe.log')
        #options.parse_command_line()
        
        # fancy mongoDB logging
        #  see http://www.joet3ch.com/blog/2011/09/08/alternative-tornado-logging/
        import logging
        from mongolog.handlers import MongoHandler
        def log_request(self, handler):
            
            log = logging.getLogger('')
            log.setLevel(logging.DEBUG)
            log.addHandler(MongoHandler.to(db='PZserver', collection='log'))
            
            if handler.get_status() < 400:
                log_method = log.info
            elif handler.get_status() < 500:
                log_method = log.warn
            else:
                log_method = log.error
            
            request_time = 1000.0 * handler.request.request_time()
            log_message = '%d %s %.2fms' % (handler.get_status(), handler._request_summary(), request_time)
            log_method(log_message)
        
        # start the server and let her run
        http_server = HTTPServer(WSGIContainer(app))
        http_server.listen(5554)
        IOLoop.instance().start()
