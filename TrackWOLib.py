'''
Not accurate.  Have 3 possibilities.
Math is wrong, Convert to RA/Dec is wrong, OR my code is wrong...
'''


#TODO 1 test with this! http://www.satflare.com/track.asp?q=25544

#TODO 4 add telescope controls LAST and test with telescope

#TODO 5 rewrite telescope software.  Number 4 dependent.

#TODO 6 measure exact time elapsed per iteration as this is important...

#TODO 7 to improve performance can remove the list objects and just have them be variables


'''
Commands
`.' ACK Used to determine the telescope axis alignment, ie. Polar or Alt-Az
:U# Used to toggle between low and high precision
coordinates.
:GD# Used to get Declination Coordinates
:GR# Used to get Right Ascension Coordinates
:Srxx:xx:xx# Used to Set right ascension coordinates
:Sdxx:xx:xx# Used to Set declination coordinates
:Q# Used to Stop all motors
:Qn#
:Qs#
:Qw#
:Qe#
:RS# Used to set the motor rate to Slew rate (Fastest).
:RM# Used to set the motor rate to Find rate (Second Fastest).
:RC# Used to set the motor rate to Centering rate (Second Slowest).
:RG# Used to set the motor rate to Guide rate (Slowest).
:FQ# Used to stop the focus motor
:FF# Used to set focus motor rate to Fast.
:FS# Used to set focus motor rate to Slow.
:F+# Used to move focus motor in.
:F-# Used to move focus motor out.
'''


'''
Keplerian elements can be found here!
Found here: http://heavens-above.com/orbit.aspx?satid=25544
'''

'''
TODO GUI:
    Rewrite the telescope stuff.  Maybe in wx or tkinter???
    Lets translate all of the available commands and build to an executable
    Make some sort of installer so that we can directly track stuff.
'''


import calendar
import datetime
import numpy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import csv
import time

from astropy import units as u
from astropy.coordinates import SkyCoord



def testLatLong():
    #https://www.mathworks.com/matlabcentral/fileexchange/7941-convert-cartesian--ecef--coordinates-to-lat--lon--alt?focused=5062924&tab=function
    #Based on above equation
    pass


def convertCartesianToRaDec(x, y, z):
    '''
    Returns in degrees in decimal notation
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
    print("RA:"+str(a) +" Dec:"+ str(d))

    c = SkyCoord(a, d, frame='icrs', unit='deg')

    #print("RA: ",c.ra.hms)
    #print("DEC: ", c.dec.degree)


    return a, d


def updatePosition(dt ,inc, RA, ecc, w, Mic, a, mu, M, v, Eout, r, i):

    M.append(Mic + dt * numpy.sqrt((mu / (a ** 3))))
    #print("DT" + str(dt))

    eN = Eout[i] - (
        (Eout[i] - ecc * numpy.degrees(numpy.sin(Eout[i])) - M[i]) / (1 - ecc * numpy.degrees(numpy.sin(Eout[i]))))
    Eout.append(eN)

    vN = 2 * numpy.degrees(numpy.arctan((((1 + ecc) / (1 - ecc)) ** 0.5) * numpy.degrees(numpy.tan(Eout[i] / 2))))
    v.append(vN)

    rN = (a * (1 - ecc ** 2)) / (1 + ecc * numpy.degrees(numpy.cos(v[i])))
    r.append(rN)

    xN = r[i] * (numpy.degrees(numpy.cos(RA)) * numpy.degrees(numpy.cos(w + v[i])) -
                 numpy.degrees(numpy.sin(RA)) *
                 numpy.degrees(numpy.sin(w + v[i])) *
                 numpy.degrees(numpy.cos(inc)))

    yN = r[i] * (numpy.degrees(numpy.sin(RA)) * numpy.degrees(numpy.cos(w + v[i])) -
                 numpy.degrees(numpy.cos(RA)) *
                 numpy.degrees(numpy.sin(w + v[i])) *
                 numpy.degrees(numpy.cos(inc)))

    zN = r[i] * (numpy.degrees(numpy.sin(inc)) * numpy.degrees(numpy.sin(w + v[i])))


    #print("Coords: ",xN, yN, zN)
    #convertCartesianToRaDec(xN, yN, zN)
    return (xN, yN, zN)





def trackData(show):

    #INclination
    inc = 51.6430#51.6407

    #Right ascension of ascending node
    RA = 344.3160#55.1621

    #Eccentricity
    ecc = 0.0003105#0.0003577

    #Argument of perigee
    w = 49.4097#19.6636

    #Mean anom in degrees
    Mic = 310.7326#340.4654

    #Earth constant
    mu = 3.986e5  # Gravitational parameter for earth

    #Revs per day
    n = 15.54044848

    #Semi major axis
    a = (mu ** (1 / 3)) / (((2 * n * numpy.pi) / (86400)) ** (2 / 3))
    print("semi " + str(a))


    Eout = []
    Eout.append(Mic)
    M = []
    v = []


    r = []


    #22 1 2018 13:23:06
    date = "30 1 2018 15:49:24 UTC" #How to format dates
    dateFormater = "%d %m %Y %H:%M:%S %Z"
    epochTime = calendar.timegm(time.strptime(date, dateFormater))

    print(epochTime)



    startTime = datetime.datetime.utcnow()
    startTime = (startTime-datetime.datetime(1970, 1, 1)).total_seconds()
    print(startTime)

    #Where to start in iteration
    dt = (int)(startTime-epochTime)
    print(dt)

    x = []
    y = []
    z = []


    #Stuff to show live plotting!
    showcase = show

    fig = None
    if showcase == True:
        plt.ion()
        fig = plt.figure()

    ax=plt.axes(projection='3d')


    i = 0
    while True:
        if (i > 5400):
            break

        start = time.time()
        data = updatePosition(dt, inc, RA, ecc, w, Mic, a, mu, M, v, Eout, r, i)
        convertCartesianToRaDec(data[0], data[1], data[2])


        #Graph the active orbit.  Might not be time critical
        if showcase == True and i > 3:
            x.append(data[0])
            y.append(data[1])
            z.append(data[2])
            ax.scatter3D(x, y, z)
            plt.show()


            # Handles time critical stuff.  Not exact.
            #end = time.time()
            if end-start > 1:
                #plt.pause(0)
                pass
            else:
                #plt.pause(1-(end-start))
                pass
            plt.pause(0.1)


            dt += 1
            i += 1
        elif showcase == False and i > 3:
            #Do telescope tracking stuff here!


            #Handles time critical stuff.  Not exact.
            end = time.time()
            if end-start < 1:
                #time.sleep(1-(end-start))
                pass

            dt += 1
            i += 1
        else:
            # Handles time critical stuff.  Not exact.
            end = time.time()
            if end-start < 1:
                #time.sleep(1-(end-start))
                pass

            dt += 1
            i += 1









if __name__ == '__main__':
    trackData(True)
    #plotStuff()
    #vega = csvReader("star.csv")
    #convertCartesianToRaDec(vega.x, vega.y, vega.z)









