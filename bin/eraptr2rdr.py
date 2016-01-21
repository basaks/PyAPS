#!/opt/local/bin/python

#########Troposphere phase contributions in InSAR########
# Computes a Radar Coordinates map of Atmospheric Phase #
# Screen from ECMWF ERA-I data                          #
# Original Fortran code, R. Jolivet                     #
# Original Python code and routines:                    #
# 							#
# This file is created by R. Jolivet, the 03-02-2012    #
#########################################################                   

#--------------------------------------------------------
#Import Routines from PyAPS
import PyAPS
import pygrib
import sys
from sys import stdout
import os.path
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
from scipy.interpolate import interp1d,LinearNDInterpolator
from scipy.integrate import cumtrapz
from pylab import *
from numpy import random

#--------------------------------------------------------
#Arguments
if len(sys.argv) < 6:
    print 'Usage: eraptr2rdr.py inputfile dem.hgt outname wvl inc'
    print 'inputfile         : Formatted text file for all ERA-I grid points available '
    print 'dem.dem           : ROI-PAC style DEM simulation file in radar geometry with .rsc file included'
    print 'outname           : Output float-32 file with phase values in radians'
    print 'wvl               : Wavelength in cm'
    print 'inc               : Incidence angle in degrees'
    print '  '
    print ' The input file is formatted as follows:'
    print '---------------------------------------------------------------------------------'
    print 'Geopotential Height (m) | Pressure (Pa) | Temperature (K) | Relative Humidity (%)'
    print '---------------------------------------------------------------------------------'
    print 'Lon :  longitude'
    print 'Lat :  latitude'
    print '---------------------------------------------------------------------------------'
    print 'z1 p1 t1 r1'
    print 'z2 p2 t2 r2'
    print '...........'
    print 'zn pn tn rn'
    print '---------------------------------------------------------------------------------'
    print 'Geopotential Height (m) | Pressure (Pa) | Temperature (K) | Relative Humidity (%)'
    print '---------------------------------------------------------------------------------'
    print 'Lon :  longitude'
    print 'Lat :  latitude'
    print '---------------------------------------------------------------------------------'
    print 'z1 p1 t1 r1'
    print 'z2 p2 t2 r2'
    print '...........'
    print 'zn pn tn rn'
    print 'And So On....'

    sys.exit(1)

#--------------------------------------------------------
# Parsing Input Parameters
#--------------------------------------------------------

print 'PROGRESS: PARSING INPUT PARAMETERS'

#ECMWF grib input
fname = sys.argv[1]
if(os.path.isfile(fname) == False):
    print 'Input File not found: ', fname
    sys.exit(1)

#Dem input
dname = sys.argv[2]
if(os.path.isfile(dname) == False):                                                                                                  
    print 'DEM File not found: ', dname
    sys.exit(1)

#outfile
oname = sys.argv[3]

#Wavelength, conversion cm to m
wvl = float(sys.argv[4])/100.0

#Incidence angle, conversion deg to rad
inc = float(sys.argv[5])*pi/180.0

#--------------------------------------------------------
# Initialize Constants 
#--------------------------------------------------------
# Reading atmo constants dictionary
cdic=PyAPS.initconst()
cdic['wvl']=wvl
cdic['inc']=inc

#Values for interpolation
minAlt = cdic['minAlt']
maxAlt = cdic['maxAlt']
nhgt = cdic['nhgt']
minAltp = cdic['minAltP']
hgt = linspace(minAlt, maxAlt, nhgt)

# Scaling for interpolation
# For rdr geom grid is about 200*200 pixels
# nght gives the spacing
#hgtscale = 0.5 for nhgt=151
hgtscale=((maxAlt-minAlt)/nhgt)/200

#--------------------------------------------------------
# Reading Dem File
#--------------------------------------------------------
[lon, lat, nx, ny, bufspc] = PyAPS.rd_rsc(dname)

#--------------------------------------------------------
#Get the min lat lon and max lat lon, roughly
#--------------------------------------------------------
minlon=lon.min()
maxlon=lon.max()
minlat=lat.min()
maxlat=lat.max()

#--------------------------------------------------------
#Read data in ERA-I text file
#--------------------------------------------------------
[lvls, latlist, lonlist, gph, tmp, vpr] = PyAPS.read_eratxt(fname, cdic)
#--------------------------------------------------------
#Compute ERA-I grid points coordinates into radar coordinates
#--------------------------------------------------------
[xi, yi] = PyAPS.glob2rdr(nx, ny, lat, lon, latlist, lonlist)

#--------------------------------------------------------
#Interpolate data in height
#--------------------------------------------------------
[Pi, Ti, Vi] = PyAPS.intP2H(lvls, hgt, gph, tmp, vpr, cdic)

#--------------------------------------------------------
#Computing Delay functions
#--------------------------------------------------------
[DDry ,DWet] = PyAPS.PTV2del(Pi, Ti, Vi, hgt, cdic)
Delfn=DDry+DWet

#--------------------------------------------------------
#Building the interpolation function
#--------------------------------------------------------
fnc = PyAPS.make3dintp(Delfn, xi, yi, hgt, hgtscale)

#--------------------------------------------------------
#Writing to File
#--------------------------------------------------------

print 'PROGRESS: WRITING TO FILE'

yarr = arange(1, ny+1)
xarr = arange(1, nx+1)
fin = open(dname, 'rb')
fout = open(oname, 'wb')

for m in range(0, 2*ny, 2):
    dem = fromfile(file=fin, dtype=float32, count=nx)
    dem = fromfile(file=fin, dtype=float32, count=nx)
    dem[dem<minAltp] = minAltp
    demy = dem.astype(float64)
    y=ones((nx,))*yarr[int(m/2)]
    d=demy/hgtscale
    llh=hstack([xarr[:, newaxis], y[:, newaxis], d[:, newaxis]])
    res = fnc(llh)
    resy = res.astype(float32)
    resy.tofile(fout)

fin.close()
fout.close()
print 'PROGRESS: COMPLETE'

