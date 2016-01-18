#!/usr/bin/env python

import argparse

import matplotlib.pyplot as plt
import numpy as np
import sys,os
import math

from astropy import units as u
from astropy.coordinates import SkyCoord,EarthLocation,AltAz
from astropy.time import Time

# my version of conversion of ha,dec,lat into alt,az
# all in degrees
# gives same results as idlutils/goddard/pro/astro/hadec2altaz.pro
# used in plate design
def hadec2altaz(ha, dec, lat) : 
    d2r = math.pi/180.
    sh = np.sin(ha*d2r)
    ch = np.cos(ha*d2r)
    sd = np.sin(dec*d2r)
    cd = np.cos(dec*d2r)
    sl = np.sin(lat*d2r)
    cl = np.cos(lat*d2r)
    x = np.array([cd*ch,cd*sh,sd])
    r = np.array([[sl,0,-cl],[0,1,0],[cl,0,sl]])
    x = r.dot(x)    
    r = np.sqrt(x[0]**2+x[1]**2)
    az = np.arctan2(-x[1],-x[0]) /d2r
    alt = np.arctan2(x[2],r) / d2r
    return alt,az

# helper function for shifting angles to desired range
def normalize_angle(angle):
    while angle <= -180:
        angle += 360
    while angle > 180:
        angle -= 360
    return angle

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ra', type=float, default=25.0,
        help='Right Ascension (degrees)')
    parser.add_argument('--dec', type=float, default=12.0,
        help='Declination (degrees)')
    parser.add_argument('--mjd', type=float, default=55000.0,
        help='Modified Julien Date (days)')
    parser.add_argument('--scan-offset', type=float, default=15.0,
        help='Scan offset (hours)')
    parser.add_argument('--scan-width', type=float, default=4.0,
        help='Scan width (hours)')
    parser.add_argument('--scan-steps', type=int, default=100,
        help='Number of sampling points')
    parser.add_argument('--astropy-apo', action='store_true',
        help='User apo observatory from astropy')
    args = parser.parse_args()

    if args.astropy_apo:
        sdss = EarthLocation.of_site('apo')
    else:
        sdss = EarthLocation(lat=32.7797556*u.deg, lon=-(105+49./60.+13/3600.)*u.deg, height=2797*u.m)

    coord = SkyCoord(ra=args.ra*u.degree, dec=args.dec*u.degree, frame='icrs')

    # scan of time
    hours = args.scan_offset + np.linspace(-0.5*args.scan_width, 0.5*args.scan_width, args.scan_steps)
    my_alt = np.zeros((hours.size))
    py_alt = np.zeros((hours.size))
    py_ha = np.zeros((hours.size))

    for i in range(hours.size):
        mjd_value = args.mjd*u.day + hours[i]*u.hour
        time = Time(val=mjd_value, scale='tai', format='mjd', location=sdss)
        # altitude from astropy
        py_alt[i] = coord.transform_to(AltAz(obstime=time, location=sdss)).alt.to(u.deg).value
        # this is supposed to be the hour angle from astropy
        py_ha[i] = time.sidereal_time('apparent').to(u.deg).value - args.ra 
        # simple rotation to get alt,az based on ha
        my_alt[i], az = hadec2altaz(py_ha[i], args.dec, sdss.latitude.to(u.deg).value) 
        print hours[i], py_ha[i], py_alt[i], my_alt[i]

    py_ha = np.array(map(normalize_angle, py_ha.tolist()))
    ii = np.argsort(py_ha)
    py_ha=py_ha[ii]
    py_alt=py_alt[ii]
    my_alt=my_alt[ii]

    fig = plt.figure(figsize=(8,6))
    plt.plot(py_ha, py_alt - my_alt, 'o', c='b')
    plt.title('Compare hadec2altaz')
    # plt.title('(ra,dec) = (%.2f,%.2f)' % (args.ra, args.dec))
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('astropy_alt - rotation_alt [deg]')
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main()



