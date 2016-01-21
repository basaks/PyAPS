######################Troposphere phase contributions in InSAR########
# Adapted from single station version provided by Romain Jolivet     #
# Python version written by Piyush Agram                             #
# Date: Jan 12, 2012                                                 #
# Modified to handle radar geometry by Romain Jolivet                #
# Date: Feb 27, 2012         				             #
######################################################################

import pygrib
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
from scipy.interpolate import interp1d,LinearNDInterpolator
from scipy.integrate import cumtrapz
import numpy as np
import math


#######Initialization of constants######################
def initconst():
    #Various constants used in the delay computation.
    constdict = {}
    constdict['k1'] = 0.776   #(K/Pa)
    constdict['k2'] = 0.716   #(K/Pa)
    constdict['k3'] = 3750    #(K^2.Pa)
    constdict['g'] = 9.81     #(m/s^2)
    constdict['Rd'] = 287.05  #(J/Kg/K)
    constdict['Rv'] = 461.495 #(J/Kg/K)
    constdict['mma'] = 29.97  #(g/mol)
    constdict['mmH'] = 2.0158 #(g/mol)
    constdict['mmO'] = 16.0   #(g/mol)

    constdict['a1w'] = 611.21 # hPa
    constdict['a3w'] = 17.502 #
    constdict['a4w'] = 32.19  # K
    constdict['a1i'] = 611.21 # hPa
    constdict['a3i'] = 22.587 #
    constdict['a4i'] = -0.7   # K
    constdict['T3'] = 273.16  # K
    constdict['Ti'] = 250.16  # K
    constdict['nhgt'] = 151   # Number of levels for interpolation
    constdict['minAlt'] = 0.0
    constdict['maxAlt'] = 15000.0
    constdict['minAltP'] = 200.0
    return constdict
###############Completed the list of constants################


###############Reading input RSC file for radar###############
def rd_rsc(inname):
    print("PROGRESS: READING %s RSC FILE" %inname)

    rpacdict = {}
    infile = open(inname+'.rsc', 'r')
    line = infile.readline()
    while line:
        llist = line.split()
        if len(llist) > 0:
            rpacdict[llist[0]] = llist[1]
        line = infile.readline()
    infile.close()

    nx = int(rpacdict['WIDTH'])
    ny = int(rpacdict['FILE_LENGTH'])
    lat = np.zeros((4, 1))
    lon = np.zeros((4, 1))
    lat[0] = float(rpacdict['LAT_REF1'])
    lon[0] = float(rpacdict['LON_REF1'])
    lat[1] = float(rpacdict['LAT_REF2'])
    lon[1] = float(rpacdict['LON_REF2'])
    lat[2] = float(rpacdict['LAT_REF3'])
    lon[2] = float(rpacdict['LON_REF3'])
    lat[3] = float(rpacdict['LAT_REF4'])
    lon[3] = float(rpacdict['LON_REF4'])
    bufspc = 1.2

    del rpacdict
    return lon, lat, nx, ny, bufspc

#######################Finished rd_rsc###########################


###############Reading input RSC file for geo###############
def geo_rsc(input_name):
    print("PROGRESS: READING %s RSC FILE" % input_name)

    rpacdict = {}
    infile = open(input_name+'.rsc', 'r')
    line = infile.readline()
    while line:
        llist = line.split()
        if len(llist) > 0:
                rpacdict[llist[0]] = llist[1]
        line = infile.readline()
    infile.close()

    nx = int(rpacdict['WIDTH'])
    ny = int(rpacdict['FILE_LENGTH'])
    lat = np.zeros((2, 1))
    lon = np.zeros((2, 1))
    lat[1] = float(rpacdict['Y_FIRST'])
    lon[0] = float(rpacdict['X_FIRST'])
    if lon[0] < 0:
        lon[0] += 360.0

    dx = float(rpacdict['X_STEP'])
    dy = float(rpacdict['Y_STEP'])

    lat[0] = lat[1] + dy*ny
    lon[1] = lon[0] + dx*nx
    if rpacdict['PROJECTION'] in ('LL', 'LATLON'):
        bufspc = 1.2
    elif rpacdict['PROJECTION'] in 'UTM':
        bufspc = 100000.0

    del rpacdict
    return lon, lat, nx, ny, bufspc

#######################Finished geo_rsc###########################


##########Interpolating to heights from Pressure levels###########
def intP2H(lvls, hgt, gph, tmp, vpr, cdic):
    """
    Interpolates the pressure level data to altitude.
    :param gph: Geopotential height
    :param tmp: Temperature
    :param vpr: Vapor pressure
    :param hgt: regular grid of height values.
    :param cdic: ?
    gph,tmp,vpr are of size(nstn,nlvls)
    """

    min_alt = cdic['minAlt']     # Hardcoded parameter.
    max_alt = cdic['maxAlt']     # Hardcoded parameter.

    print 'PROGRESS: INTERPOLATING FROM PRESSURE TO HEIGHT LEVELS'
    nstn = gph.shape[1]           # Number of stations
    nhgt = len(hgt)               # Number of height points
    Presi = np.zeros((nstn, nhgt))
    Tempi = np.zeros((nstn, nhgt))
    Vpri = np.zeros((nstn, nhgt))

    for i in range(nstn):
        temp = gph[:, i]       # Obtaining height values
        hx = temp.copy()
        s_flag = False
        e_flag = False
        if hx.min() > min_alt:          # Add point at start
            s_flag = True
            hx = np.concatenate((hx, [min_alt-1]))

        if hx.max() < max_alt:		#Add point at end
            e_flag = True
            hx = np.concatenate(([max_alt+1], hx))

        hx = -hx             #Splines needs monotonically increasing.

        hy = lvls.copy()     #Interpolating pressure values
        if s_flag:
            val = hy[-1] + (hx[-1] - hx[-2])*(hy[-1] - hy[-2])/(hx[-2]-hx[-3])
            hy = np.concatenate((hy, [val]))
        if e_flag:
            val = hy[0] - (hx[0] - hx[1])*(hy[0] - hy[1])/(hx[1]-hx[2])
            hy = np.concatenate(([val], hy))

        tck = interp1d(hx, hy, kind='cubic')

        temp = tck(-hgt)      #Again negative for consistency with hx
        if s_flag & e_flag:
            Presi[i, :] = temp[1:nhgt+1].copy()
        elif s_flag:
            Presi[i, :] = temp[0:nhgt].copy()
        elif e_flag:
            Presi[i, :] = temp[1:nhgt+1].copy()
        else:
            Presi[i, :] = temp.copy()

        del temp

        temp = tmp[:, i]		#Interpolating temperature
        hy = temp.copy()
        if s_flag:
            val = hy[-1] + (hx[-1] - hx[-2])*(hy[-1] - hy[-2])/(hx[-2]-hx[-3])
            hy = np.concatenate((hy, [val]))
        if e_flag:
            val = hy[0] - (hx[0] - hx[1]) * (hy[0] - hy[1])/(hx[1]-hx[2])
            hy = np.concatenate(([val], hy))

        tck = interp1d(hx, hy, kind='cubic')
        temp = tck(-hgt)
        if s_flag & e_flag:
            Tempi[i, :] = temp[1:nhgt+1].copy()
        elif s_flag:
            Tempi[i, :] = temp[0:nhgt].copy()
        elif e_flag:
            Tempi[i, :] = temp[1:nhgt+1].copy()
        else:
            Tempi[i, :] = temp.copy()

        del temp

        temp = vpr[:, i]          #Interpolating vapor pressure
        hy = temp.copy()
        if s_flag:
            val = hy[-1] + (hx[-1] - hx[-2])*(hy[-1] - hy[-2])/(hx[-2]-hx[-3])
            hy = np.concatenate((hy, [val]))
        if e_flag:
            val = hy[0] - (hx[0] - hx[1]) * (hy[0] - hy[1])/(hx[1]-hx[2])
            hy = np.concatenate(([val], hy))

        tck = interp1d(hx, hy, kind='cubic')
        temp = tck(-hgt)
        if s_flag & e_flag:
            Vpri[i, :] = temp[1:nhgt+1].copy()
        elif s_flag:
            Vpri[i, :] = temp[0:nhgt].copy()
        elif e_flag:
            Vpri[i, :] = temp[1:nhgt+1].copy()
        else:
            Vpri[i, :] = temp.copy()

        del temp


    return Presi,Tempi,Vpri
###########Completed interpolation to height levels #####################


###########Computing the delay function ###############################
def PTV2del(Presi,Tempi,Vpri,hgt,cdict):
    #Computes the delay function given P,T and Vpr.
    #Computes refractive index at each altitude and
    #integrates the delay using cumtrapz.
    #hgt refers to the grid on which P,T, Vpr were interpolated.
    #cdict is a dictionary of constants.
    print 'PROGRESS: COMPUTING DELAY FUNCTIONS'
    nhgt = len(hgt)			#Number of height points
    nstn = Presi.shape[0]		#Number of stations
    WonT = Vpri/Tempi
    WonT2 = WonT/Tempi

    k1 = cdict['k1']
    Rd = cdict['Rd']
    Rv = cdict['Rv']
    k2 = cdict['k2']
    k3 = cdict['k3']
    inc = cdict['inc']
    wvl = cdict['wvl']
    g = cdict['g']

    #Dry delay
    DDry2 = np.zeros((nstn, nhgt))
    for i in range(nstn):
        DDry2[i, :] = k1*Rd*(Presi[i, :] - Presi[i, -1])*1.0e-6/(math.cos(inc)*g)

    DDry2 = 4.0*math.pi*DDry2/wvl

    #Wet delay
    S1 = cumtrapz(WonT, x=hgt, axis=-1)
    val = 2*S1[:, -1]-S1[:, -2]
    val = val[:, None]
    S1 = np.concatenate((S1, val), axis=-1)
    del WonT

    S2 = cumtrapz(WonT2, x=hgt, axis=-1)
    val = 2*S2[:, -1]-S2[:, -2]
    val = val[:, None]
    S2 = np.concatenate((S2, val), axis=-1)
    DWet2 = 1.0e-6*((k2-k1*Rd/Rv)*S1+k3*S2)/math.cos(inc)
    DWet2 = -4*math.pi*DWet2/wvl
    for k in range(nstn):
        DWet2[k, :] = DWet2[k, :] - DWet2[k, -1]

    return DDry2, DWet2
####################Completed delay function#################


##########Conversion from Geo coordinate to Radar coordinates######
def glob2rdr(nx, ny, lat, lon, latl, lonl):
    #Transfert these latlon coordinates into radar geometry (xi,yi)
    #with a simple linear transformation given the first pixel
    #and the pixel spacing of the simulation

    # Mapping function is :
    # range = a1*lat+b1*lon+c1
    # azimu = a2*lat+b2*lon+c2
    # First  point is (1,1)   i.e. Near Range, First Lane <==> Lat[1],Lon[1]
    # Second point is (nx,1)  i.e. Far Range, First Lane <==> Lat[2],Lon[2]
    # Third  point is (1,ny)  i.e. Near Range, Last Lane <==> Lat[3],Lon[3]
    # Fourth point is (nx,ny) i.e. Far Range, Last Lane <==> Lat[4],Lon[4]

    #	A=array([(lat[0], lon[0], 1.),(lat[1], lon[1], 1.),(lat[2], lon[2], 1.),(lat[3], lon[3], 1.)])
    A = np.hstack([lat, lon, np.ones((4, 1))])
    b = np.array([[1, 1], [nx, 1], [1, ny], [nx, ny]])
    mfcn = np.linalg.lstsq(A, b)[0]

    #Get grid points xi yi coordinates from this mapping function
    nstn = len(latl)
    A = np.array([latl, lonl, np.ones(nstn)]).conj().T
    xi = np.dot(A, mfcn[:, 0])
    yi = np.dot(A, mfcn[:, 1])

    plotflag = 'n'
    if plotflag in ('y', 'Y'):
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(lonl, latl, s=8, c='k')
        xline=[lon[0], lon[1], lon[3], lon[2], lon[0]]
        yline=[lat[0], lat[1], lat[3], lat[2], lat[0]]
        plt.plot(xline, yline, '-r')
        plt.title('Area of interest and %d stations used' % nstn)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        plt.subplot(212)
        plt.scatter(xi, yi, s=8, c='k')
        p = ptch.Rectangle((1, 1), nx, ny, edgecolor="Red", fill=False)
        plt.gca().add_patch(p)
        plt.title('Area of interest in Radar Geometry')
        plt.xlabel('Range')
        plt.ylabel('Azimuth')
        plt.show()

    return xi, yi
#################Completed transforming geo 2 radar##################

####Setting up 3D interpolation function in geo/xy coordinates#######
def make3dintp(Delfn,lonlist,latlist,hgt,hgtscale):

    ##Delfn   = Ddry + Dwet. Delay function.
    ##lonlist = list of lons for stations. / x
    ##latlist = list of lats for stations. / y
    nstn = Delfn.shape[0]
    nhgt = Delfn.shape[1]
    xyz = np.zeros((nstn*nhgt, 3))
    Delfn = np.reshape(Delfn, (nstn*nhgt, 1))
    count = 0
    for m in range(nstn):
        for n in range(nhgt):
            xyz[count,0] = lonlist[m]
            xyz[count,1] = latlist[m]
            xyz[count,2] = hgt[n]/hgtscale     #For same grid spacing as lat/lon
            count += 1

    xyz[:, 2] = xyz[:, 2] + 0.001*np.random.rand((nstn*nhgt))/hgtscale  # For unique Delaunay
    del latlist
    del lonlist
    del hgt
    print 'PROGRESS: BUILDING INTERPOLATION FUNCTION'
    fnc = LinearNDInterpolator(xyz,Delfn)

    return fnc
###########Completed 3D interpolation in geo coordinates############

#############Clausis-Clayperon for ECMWF###########################
def cc_era(tmp, cdic):
    #Input are :
    # tmp: a vector containing the temperature profiles at each station (size is (n-levels,n-stations))
    # cdic: a dictionnary containing the constants used in the Clausius-Clapeyron law

    #Output is :
    # esat: Water vapor saturation partial pressure (size is (n-levels,n-stations))
    a1w = cdic['a1w']
    a3w = cdic['a3w']
    a4w = cdic['a4w']
    a1i = cdic['a1i']
    a3i = cdic['a3i']
    a4i = cdic['a4i']
    T3 = cdic['T3']
    Ti = cdic['Ti']

    esatw = a1w*np.exp(a3w*(tmp-T3)/(tmp-a4w))
    esati = a1i*np.exp(a3i*(tmp-T3)/(tmp-a4i))
    esat = esati.copy()
    for k in range(len(tmp)):
        if tmp[k] >= T3:
            esat[k] = esatw[k]
        elif tmp[k] <= Ti:
            esat[k] = esati[k]
        else:
            wgt = (tmp[k]-Ti)/(T3-Ti)
            esat[k] = esati[k] + (esatw[k]-esati[k])*wgt*wgt

    return esat
###############Completed CC_ERA#########################################


#############Clausius-Clapeyron for ECMWF as used in Jolivet et al 2011#
def cc_eraorig(tmp, cdic):
    #This routine takes temperature profiles and returns Saturation water vapor partial pressure
    #Using the Clausius-Clapeyron law applied in Jolivet et al. 2011,GRL, doi:10.1029/2011GL048757
    #It can be used in case you are using Relative Humidity from ECMWF models

    #Input:
    # tmp: temperature profiles at a station (size(n-stations))
    # cdic: Dictionnary containing parameter values

    #Output:
    # esat: saturation water vapor partial pressure at each station (size(n-stations))

    a1w = cdic['a1w']
    T3 = cdic['T3']
    Rv = cdic['Rv']

    esat = a1w*np.exp((2.5e6/Rv) * ((1/T3) - (1/tmp)))

    return esat
###############Completed CC_ERAORIG#####################################

#############Clausius-Clapeyron for NARR ###########
def cc_narr(tmp,cdic):
    #This routine takes temperature profiles and returns Saturation water vapor partial pressure
    #Using the Clausius-Clapeyron law applied in Jolivet et al. 2011,GRL, doi:10.1029/2011GL048757
    #It can be used in case you are using Relative Humidity from ECMWF models

    #Input:
    # tmp: temperature profiles at a station (size(n-stations))
    # cdic: Dictionnary containing parameter values

    #Output:
    # esat: saturation water vapor partial pressure at each station (size(n-stations))

    a1w = cdic['a1w']
    a3w = cdic['a3w']
    a4w = cdic['a4w']
    T3 = cdic['T3']
    Rv = cdic['Rv']
    esat = a1w*np.exp(a3w*(tmp-T3)/(tmp-a4w))
    return esat
###############Completed CC_NARR#####################################


########Read in ERA data from a given ERA Interim file##################
def get_era(fname,minlat,maxlat,minlon,maxlon,cdic):
    ###Read data from ERA grib file.
    ###Note that Lon values should be between [0-360].
    ###GRB file with weather model data can be downloaded from http://dss.ucar.edu/datasets/ds627.0/
    print 'PROGRESS: READING GRIB FILE'
    lvls = np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000])
    nlvls = len(lvls)

    alpha = cdic['Rv']/cdic['Rd']
    gphind = np.arange(nlvls)*12+1

    grbs = pygrib.open(fname)
    grbs.seek(gphind[0])
    grb = grbs.read(1)[0]
    lats, lons = grb.latlons()
    g = cdic['g']
    mask = (lats > minlat) & (lats < maxlat) \
        & (lons > minlon) & (lons < maxlon)
    [ii, jj] = np.where(mask)
    del mask
    latlist = lats[ii, jj]
    lonlist = lons[ii, jj]
    nstn = np.size(ii)

    ####Create arrays for 3D storage
    gph = np.zeros((nlvls, nstn))   #Potential height
    tmp = gph.copy()                  #Temperature
    vpr = gph.copy()                  #Vapor pressure
    print 'Number of stations:', nstn
    lvls = np.multiply(lvls, 100.0)              #Conversion to absolute pressure
    for i in range(nlvls):
        grbs.seek(gphind[i])   #Reading potential height.
        grb = grbs.read(3)
        val = grb[0].values
        gph[i, :] = val[ii, jj]/g

        val = grb[1].values   #Reading temperature
        temp = val[ii, jj]
        tmp[i, :] = temp

        #For use with relative humidity
        esat = cc_era(temp, cdic)
        grbs.seek(gphind[i]+6)
        grb = grbs.read(1)
        val = grb[0].values
        temp = val[ii, jj]/100.0
        vpr[i, :] = temp*esat

    #For use with Specific Humidity
    #		val = grb[2].values  #Specific humidity
    #	        temp = val[ii,jj]
    #		vpr[i,:] = temp*lvls[i]*alpha/(1+(alpha - 1)*temp)

    #Original version by Jolivet et al 2011 with relative humidity
    #		esat = cc_eraorig(temp,cdic)
    #		grbs.seek(gphind[i]+6)
    #		grb = grbs.read(1)
    #		val = grb[0].values
    #		temp = val[ii,jj]/100.0
    #		vpr[i,:] = temp*esat

    return lvls, latlist, lonlist, gph, tmp, vpr

###############Completed GET_ERA########################################


##############Read Input text file #####################################
def read_eratxt(fname, cdic):

    lvls = []
    latlist = []
    lonlist = []
    gpht = []
    tmpt = []
    reht = []

    g = cdic['g']

    f = open(fname,'r')
    tmp = f.readlines()
    i = 0
    nstn = 0
    maxloop = int(np.size(tmp))
    while i < maxloop:
        if tmp[i][0] == '-':
            nstn += 1
            lonlat = tmp[i+3].rsplit(' ')
            lonlist.append(float(lonlat[3]))
            latlist.append(float(lonlat[9]))
            i += 5
            new = 'y'
        else:
            if new in ('y'):
                n = 1
                val = tmp[i].rsplit(' ')
                gpht.append(float(val[0]))
                lvls.append(float(val[1]))
                tmpt.append(float(val[2]))
                reht.append(float(val[3]))
                i += 1
                new = 'n'
            else:
                n += 1
                val = tmp[i].rsplit(' ')
                gpht.append(float(val[0]))
                lvls.append(float(val[1]))
                tmpt.append(float(val[2]))
                reht.append(float(val[3]))
                i += 1

    gpht = np.array(gpht)/g
    gph = np.flipud(gpht.reshape((n, nstn), order='F'))
    del gpht

    tmpt = np.array(tmpt)
    esat = cc_eraorig(tmpt, cdic)
    tmp = np.flipud(tmpt.reshape((n, nstn), order='F'))
    del tmpt

    vprt = (np.array(reht)/100.)*esat
    vpr = np.flipud(vprt.reshape((n, nstn), order='F'))
    del vprt
    del esat

    lvls = np.flipud(np.array(lvls))
    lvls = lvls[0:n]

    lonlist = np.array(lonlist)
    latlist = np.array(latlist)

    return lvls, latlist, lonlist, gph, tmp, vpr
