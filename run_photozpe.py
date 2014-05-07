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
    
    # start the simple flask server and let her see the larger world
    app.run(host='0.0.0.0', port=port, debug=True ) #set debug=True for debugging interface
