import numpy as np
import netCDF4
import scipy.ndimage
try:
    import pandas as pd
    use_pandas = True
except ImportError:
    use_pandas = False
    pass

import utils
import poppy.grid

def get_amoc(ncfiles, latlim=(30,60), zlim=(500,9999)):
    """Retrieve AMOC time series from a set of CESM/POP model output files

    Parameters
    ----------
    ncfiles : list of str
        paths to input files
    latlim : tuple (float,float)
        Latitude limits between which to find the maximum AMOC
    zlim : tuple (float,float)
        Depth limits between which to find the maximum AMOC

    Returns
    -------
    Atlantic meridional overturning circulation (AMOC) maximum between 
    30N and 60N below 500 m water depth (following Shields et al. 2012).

    Note
    ----
    Only works with POP data that has the diagnostic variable 'MOC' included.

    """
    n = len(ncfiles)
    print('Processing {} files ...'.format(n))

    maxn = utils.get_ulimitn()

    with netCDF4.Dataset(ncfiles[0]) as ds:
        dsvar = ds.variables
        zax = dsvar['moc_z'][:]/100.
        kza = np.argmin(np.abs(zax-zlim[0]))
        kzo = np.argmin(np.abs(zax-zlim[1]))
        nz = kzo-kza+1

        latax = dsvar['lat_aux_grid'][:]
        ja = np.argmin(np.abs(latax-latlim[0]))
        jo = np.argmin(np.abs(latax-latlim[1]))        
        nlat = jo-ja+1
        
    if n <= maxn:
        with netCDF4.MFDataset(ncfiles) as ds:
            dsvar = ds.variables
            timeax = netCDF4.num2date(dsvar['time_bound'][:,0],
                                      dsvar['time_bound'].units,'proleptic_gregorian')
            amoc = dsvar['MOC'][:,1,0,kza:kzo+1,ja:jo+1]
    else:
        timeax = np.zeros(n,'object')
        amoc = np.zeros((n,nz,nlat))
        for i,fname in enumerate(sorted(ncfiles)):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = netCDF4.num2date(dsvar['time_bound'][0,0],
                                             dsvar['time_bound'].units,'proleptic_gregorian')
                amoc[i,...] = dsvar['MOC'][0,1,0,kza:kzo+1,ja:jo+1]
                
    window_size = 12
    if window_size > 1:
        maxmeanamoc = np.max(np.max(scipy.ndimage.convolve1d(
            amoc,weights=np.ones(int(window_size))/float(window_size),
            axis=0,mode='constant',cval=np.nan),axis=-1),axis=-1)
        maxmeanamoc[:window_size+1] = np.nan
        maxmeanamoc[-window_size:] = np.nan
    else:
        maxmeanamoc = np.max(np.max(amoc,axis=-1),axis=-1)

    return maxmeanamoc, timeax


componentnames = {
    0 : 'Total',
    1 : 'Eulerian-Mean Advection',             
    2 : 'Eddy-Induced Advection (bolus) + Diffusion',
    3 : 'Eddy-Induced (bolus) Advection',
    4 : 'Submeso Advection',
    }


def get_mht(ncfiles, latlim=(30,60), component=0):
    """Get MHT time series from CESM/POP data
    
    Parameters
    ----------
    ncfiles : list of str
        paths to input files
    latlim : tup
        latitude limits for maximum
    component : int
        see metrics.componentnames
    """
    n = len(ncfiles)
    print('Processing {} files ...'.format(n))
    maxn = utils.get_ulimitn()

    with netCDF4.Dataset(ncfiles[0]) as ds:
        latax = ds.variables['lat_aux_grid'][:]
        ja = np.argmin(np.abs(latax-latlim[0]))
        jo = np.argmin(np.abs(latax-latlim[1]))
        nlat = jo-ja+1
        
    if n <= maxn:
        with netCDF4.MFDataset(ncfiles) as ds:
            dsvar = ds.variables
            timeax = netCDF4.num2date(dsvar['time_bound'][:,0],
                                      dsvar['time_bound'].units,'proleptic_gregorian')
            nheat = dsvar['N_HEAT'][:,0,component,ja:jo+1]
    else:
        timeax = np.zeros(n,'object')
        nheat = np.zeros((n,nlat))
        for i,fname in enumerate(sorted(ncfiles)):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = netCDF4.num2date(dsvar['time_bound'][0,0],
                                             dsvar['time_bound'].units,'proleptic_gregorian')
                nheat[i,:] = dsvar['N_HEAT'][0,0,component,ja:jo+1]
                
    window_size = 12
    maxmeannheat = np.max(scipy.ndimage.convolve1d(
        nheat,weights=np.ones(int(window_size))/float(window_size),
        axis=0,mode='constant',cval=np.nan),axis=-1)

    return maxmeannheat, timeax


def get_mst(ncfiles, lat0=55, component=0):
    """Get MST time series from CESM/POP data
    
    Parameters
    ----------
    ncfiles : list of str
        paths to input files
    lat0 : float
        latitude to take the mean at
    component : int
        see metrics.componentnames
    """
    n = len(ncfiles)
    print('Processing {} files ...'.format(n))
    maxn = utils.get_ulimitn()

    with netCDF4.Dataset(ncfiles[0]) as ds:
        dsvar = ds.variables
        latax = dsvar['lat_aux_grid'][:]
        j0 = np.argmin(np.abs(latax-lat0))
        
    if n <= maxn:
        with netCDF4.MFDataset(ncfiles) as ds:
            dsvar = ds.variables
            timeax = netCDF4.num2date(dsvar['time_bound'][:,0],
                                      dsvar['time_bound'].units,'proleptic_gregorian')
            nsalt = dsvar['N_SALT'][:,0,component,j0]
    else:
        timeax = np.zeros(n,'object')
        nsalt = np.zeros(n)
        for i,fname in enumerate(sorted(ncfiles)):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = netCDF4.num2date(dsvar['time_bound'][0,0],
                                             dsvar['time_bound'].units,'proleptic_gregorian')
                nsalt[i] = dsvar['N_SALT'][0,0,component,j0]
                
    window_size=12
    window = np.ones(int(window_size))/float(window_size)
    meannsalt = np.convolve(nsalt,window,'same')
    meannsalt[:window_size+1] = np.nan
    meannsalt[-window_size:] = np.nan

    return meannsalt, timeax


def get_timeseries(ncfiles, varn, grid='T', reducefunc=np.mean, latlim=(), lonlim=()):
    """Get time series of any 2D POP field reduced by a numpy function
    
    Parameters
    ----------
    ncfiles : list of str
        paths to input files
    varn : str
        variable name
    grid : str ('T' or 'U')
        which grid the variable is on
    reducefunc : function
        function to reduce the selected region
        must have an 'axis' parameter to select reduction axis
    latlim : tup
        latitude limits for maximum
    lonlim : tup
        longitude limits for maximum
    """
    n = len(ncfiles)
    print('Processing {} files ...'.format(n))
    maxn = utils.get_ulimitn()

    with netCDF4.Dataset(ncfiles[0]) as ds:
        mask = poppy.grid.get_mask_lonlat(ds,lonlim=lonlim,latlim=latlim,grid=grid)
        mask &= ds.variables['KM'+grid][:]>0
        jj,ii = np.where(mask)
        units = ds.variables[varn].units
        
    if n <= maxn:
        with netCDF4.MFDataset(ncfiles) as ds:
            dsvar = ds.variables
            timeax = netCDF4.num2date(dsvar['time_bound'][:,0],
                                      dsvar['time_bound'].units,'proleptic_gregorian')
            tseries = reducefunc(dsvar[varn][:,jj,ii],axis=-1)
    else:
        timeax = np.zeros(n,'object')
        tseries = np.zeros((n))
        for i,fname in enumerate(sorted(ncfiles)):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = netCDF4.num2date(dsvar['time_bound'][0,0],
                                             dsvar['time_bound'].units,'proleptic_gregorian')
                tseries[i] = reducefunc(dsvar[varn][0,0,jj,ii])
                
    if use_pandas:
        index = pd.Index(utils.datetime_to_decimal_year(timeax), name='Model year')
        ts = pd.Series(tseries, index=index, name='{} ({})'.format(varn, units))
        ts.latlim = latlim
        ts.lonlim = lonlim
        ts.varn = varn
        ts.reducefunc = str(reducefunc)
        return ts
    else:
        return tseries, timeax

