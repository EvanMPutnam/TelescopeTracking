"""
Description.  Tracking software with the help of libraries.
"""
#TODO: we need to debug the convert to CartesianToRaDec function.  Need to test.
from pyorbital import tlefile
from pyorbital.orbital import Orbital
from datetime import datetime
import time
import numpy
from astropy.coordinates import SkyCoord


def convertCartesianToRaDec(x, y, z):
    '''
    Returns RA/DECin degrees in decimal notation
    :param x:
    :param y:
    :param z:
    :return:
    '''
    a = None
    if(y > 0 and x > 0):
        a = numpy.degrees(numpy.arctan(y/x))
    elif (x < 0):
        a = 180+numpy.degrees(numpy.arctan(y/x))
    elif (x >0 and y < 0):
        a = 360 + numpy.degrees(numpy.arctan(y/x))
    r = numpy.sqrt(x**2+y**2+z**2)
    d = numpy.degrees(numpy.arcsin(z/r))
    #print("RA:"+str(a) +" Dec:"+ str(d))
    c = SkyCoord(a, d, frame='icrs', unit='deg')
    print()
    print(c.ra.hms)
    print(c.dec.dms)
    print(str(a) + ' '+str(d))
    return a, d


def track():
    #tle = tlefile.read('ISS', 'orbit.txt')
    orbit = Orbital('ISS', 'orbit.txt')
    while True:
        now = datetime.utcnow()
        #print(now)
        coords = (orbit.get_position(now)[0])
        #print(orbit.get_lonlatalt(now))
        #print(coords)
        convertCartesianToRaDec(coords[0], coords[1], coords[2])


        time.sleep(1)



track()




