#!/opt/local/bin/python

#########Troposphere phase contributions in InSAR########
# Computes a Radar Coordinates map of Atmospheric Phase #
# Screen from ECMWF ERA-I data                          #
# Original Fortran code, R. Jolivet                     #
# Original Python code and routines: P. Agram           #
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
from numpy import random
import math
import numpy as np

#--------------------------------------------------------
#Arguments
if len(sys.argv) < 6:
    print 'Usage: era2rdr.py ei.oper.an.pl.regn128sc.yyyymmddhh dem.hgt outname wvl inc'
    print 'ei.oper.an.pl.regn128sc.yyyymmddhh : GRB file with weather model data. Can be downloaded from http://dss.ucar.edu/datasets/ds627.0/ '
    print 'dem.dem                         : ROI-PAC style DEM simulation file in radar geometry with .rsc file included'
    print 'outname                         : Output float-32 file with phase values in radians'
    print 'wvl                             : Wavelength in cm'
    print 'inc                             : Incidence angle in degrees'
    sys.exit(1)

#--------------------------------------------------------
# Parsing Input Parameters
#--------------------------------------------------------

print 'PROGRESS: PARSING INPUT PARAMETERS'

#ECMWF grib input
fname = sys.argv[1]
if not os.path.isfile(fname):
    print 'ERA File not found: ', fname
    sys.exit(1)

#Dem input
dname = sys.argv[2]
if not os.path.isfile(dname):
    print 'DEM File not found: ', dname
    sys.exit(1)

#outfile
oname = sys.argv[3]

#Wavelength, conversion cm to m
wvl = float(sys.argv[4])/100.0

#Incidence angle, conversion deg to rad
inc = float(sys.argv[5])*math.pi/180.0

#--------------------------------------------------------
# Initialize Constants 
#--------------------------------------------------------
# Reading atmo constants dictionary
cdic = PyAPS.initconst()
cdic['wvl'] = wvl
cdic['inc'] = inc

#Values for interpolation
minAlt = cdic['minAlt']
maxAlt = cdic['maxAlt']
nhgt = cdic['nhgt']
minAltp = cdic['minAltP']
hgt = np.linspace(minAlt, maxAlt, nhgt)

# Scaling for interpolation
# For rdr geom grid is about 200*200 pixels
# nght gives the spacing
#hgtscale = 0.5 for nhgt=151
hgtscale = ((maxAlt-minAlt)/nhgt)/200

#--------------------------------------------------------
# Reading Dem File
#--------------------------------------------------------
[lon, lat, nx, ny, bufspc] = PyAPS.rd_rsc(dname)

#--------------------------------------------------------
#Get the min lat lon and max lat lon, roughly
#--------------------------------------------------------
minlon = lon.min()
maxlon = lon.max()
minlat = lat.min()
maxlat = lat.max()

#--------------------------------------------------------
#Read data in ERA-I GRIB file
#--------------------------------------------------------
[lvls, latlist, lonlist, gph, tmp, vpr] = PyAPS.get_era(
    fname, minlat - bufspc, maxlat + bufspc, minlon - bufspc,
    maxlon + bufspc, cdic)

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
[DDry, DWet] = PyAPS.PTV2del(Pi, Ti, Vi, hgt, cdic)
Delfn = DDry+DWet

#--------------------------------------------------------
#Building the interpolation function
#--------------------------------------------------------
fnc = PyAPS.make3dintp(Delfn, xi, yi, hgt, hgtscale)

#--------------------------------------------------------
#Writing to File
#--------------------------------------------------------

print 'PROGRESS: WRITING TO FILE'

yarr = np.arange(1, ny+1)
xarr = np.arange(1, nx+1)
fin = open(dname, 'rb')
fout = open(oname, 'wb')
for m in range(0, 2*ny, 2):
    dem = np.fromfile(file=fin, dtype=np.float32, count=nx)
    dem[dem < minAltp] = minAltp
    demy = dem.astype(np.float64)
    y = np.ones((nx,))*yarr[int(m/2)]
    d = demy/hgtscale
    llh = np.hstack([xarr[:, np.newaxis], y[:, np.newaxis], d[:, np.newaxis]])
    res = fnc(llh)
    resy = res.astype(np.float32)
    resy.tofile(fout)

fin.close()
fout.close()
print 'PROGRESS: COMPLETE'

