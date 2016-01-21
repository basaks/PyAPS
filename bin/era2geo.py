#!/opt/local/bin/python

######################Troposphere phase contributions in InSAR########
# Adapted from single station version provided by Romain Jolivet     #
# Adapted to directly use EI.OPER.AN.PL GRIB files and a single DEM File    #
# Python version written by Piyush Agram                             #
# Date: Jan 12, 2012                                                 #
######################################################################

import sys
import os.path
import math
import numpy as np
import matplotlib.pyplot as plt
import PyAPS

if len(sys.argv) < 6:
    print 'Usage: era2ph.py ei.oper.an.pl.regn128sc.yyyymmddhh dem.dem outname wvl inc'
    print 'ei.oper.an.pl.regn128sc.yyyymmddhh : GRB file with weather model data. Can be downloaded from http://dss.ucar.edu/datasets/ds627.0/ '
    print 'dem.dem                         : ROI-PAC style DEM file with .rsc file included'
    print 'outname                         : Output float-32 file with phase values in radians'
    print 'wvl                             : Wavelength in cm'
    print 'inc                             : Incidence angle in degrees'
    sys.exit(1)

###############Parsing input paramters##############
print 'PROGRESS: PARSING INPUT PARAMETERS'
fname = sys.argv[1]
if not os.path.isfile(fname):
    print 'ERA File not found: ', fname
    sys.exit(1)

dname = sys.argv[2]

if not os.path.isfile(dname):
    print 'DEM File not found: ', dname
    sys.exit(1)

if not os.path.isfile(dname+'.rsc'):
    print 'DEM RSC File not found: ', dname
    sys.exit(1)

oname = sys.argv[3]
wvl = float(sys.argv[4])/100.0          #Conversion from cm to meters
inc = float(sys.argv[5])*math.pi/180.0	#Conversion to radians
####################Completed parsing inputs


#####Reading DEM.rsc file ############################
[lon, lat, nx, ny, bufspc] = PyAPS.geo_rsc(dname)

##############Completed reading DEM.rsc file##############
cdict=PyAPS.initconst()
cdict['wvl'] = wvl
cdict['inc'] = inc

plotflag = 'n'
hgt = np.linspace(cdict['minAlt'], cdict['maxAlt'], cdict['nhgt']) #Heights for interpolation

# Scaling for interpolation
# For geo geom grid is about 0.703*0.703 
# nght gives the spacing
#hgtscale = 142.25 for nhgt=151
hgtscale=((cdict['maxAlt']-cdict['minAlt'])/cdict['nhgt'])/0.703


#################Reading in Weather Data from NARR file#############
[lvls,latlist,lonlist,gph,tmp,vpr] = PyAPS.get_era(fname,lat[0]-bufspc,lat[1]+bufspc,lon[0]-bufspc,lon[1]+bufspc,cdict)

##########Interpolating to heights from Pressure levels###########
[Presi,Tempi,Vpri] = PyAPS.intP2H(lvls,hgt,gph,tmp,vpr,cdict)
###########Computing the delay function ###############################
[DDry,DWet] = PyAPS.PTV2del(Presi,Tempi,Vpri,hgt,cdict)
Delfn = DDry+DWet
del DDry
del DWet

#########Reading in DEM and writing output #################
#######################make_geomap.f90#######################
fnc = PyAPS.make3dintp(Delfn,lonlist,latlist,hgt,hgtscale)

print 'PROGRESS: INTERPOLATION FUNCTION READY'
minAltp = cdict['minAltP']
laty = np.linspace(lat[1], lat[0], ny)
lonx = np.linspace(lon[0], lon[1], nx)

fout = open(oname, 'wb')

for m in range(ny):
    dem = np.fromfile(dname, dtype=np.int16, count=nx)
    dem[dem<minAltp] = minAltp
    demy = dem.astype(np.float64)
    llh = np.zeros((nx, 3))
    llh[:, 0] = lonx
    llh[:, 1] = laty[m]
    llh[:, 2] = demy/hgtscale
    res = fnc(llh)
    resy = res.astype(np.float32)
    resy.tofile(fout)


fout.close()
print 'PROGRESS: COMPLETED'
if plotflag in ('y', 'Y'):
    plt.show()

