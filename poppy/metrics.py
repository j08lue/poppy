import numpy as np
import netCDF4
import scipy.ndimage
import subprocess
import traceback
import calendar
try:
    import pandas as pd
    use_pandas = True
except ImportError:
    use_pandas = False
    print 'Pandas could not be imported. Functions will return data as tuple '\
            '(tseries, timeaxis).'
    pass

from . import grid as poppygrid

### HELP FUNCTIONS

def _get_timeax(ds):
    dsvar = ds.variables
    time = dsvar['time']
    timedata = time[:]
    timeunits = time.units
    if timeunits.startswith('days since 0000'):
        timeunits = timeunits.replace('days since 0000', 'days since 0001')
        timedata -= (365+16)
    return netCDF4.num2date(timedata, units=timeunits, calendar=time.calendar)

def _get_decimal_year(ds):
    dsvar = ds.variables
    time = dsvar['time']
    if time.units != 'days since 0000-01-01 00:00:00':
        raise NotImplementedError()
    return time[:]/365.

def get_ulimitn(default=1e3):
    """Try to get the maximum number of open files on the system. Works with Unix."""
    return 100 # There are performance issues with netCDF4.MFDataset
    try:
        maxn = int(subprocess.check_output(['ulimit -n'],shell=True))
    except:
        maxn = default
        traceback.print_exc()
    return maxn

def datetime_to_decimal_year(dd, ndays=None):
    """Compute decimal year from datetime instances 

    Parameters
    ----------
    dd : datetime instance or iterable of such
        dates to convert
    ndays : int, optional 
        number of days in year
        if not given, it will be determined from the date
        (assuming proleptic gregorian)
    """
    def _convert(d, ndays):
        if ndays is None:
            ndays = [365,364][calendar.isleap(d.year)]
        nsec = float(ndays*24*3600)
        doy = int(d.strftime('%j')) # ugly but currently ncdftime-safe
        return d.year + ((doy-1)*24*3600 + d.hour*3600 + d.minute*60 + d.second) / nsec
    _convert_vec = np.vectorize(_convert)
    return np.squeeze(_convert_vec(dd, ndays))[()]

def _nfiles_diag(n):
    if n == 0:
        raise ValueError('No files found. Check your glob pattern.')
    else:
        print 'Processing {} files ...'.format(n)

def _pandas_add_meta_data(ts, meta):
    return ts # NOT WORKING PROPERLY ANYWAYS!
    for k,v in meta.iteritems():
        ts.k = v
        ts._metadata.append(k)
    return ts

def _pandas_copy_meta_data(source, target, addmeta={}):
    return target # NOT WORKING PROPERLY ANYWAYS!
    for k in source._metadata:
        try:
            setattr(target, k, getattr(source, k))
            target._metadata.append(k)
        except AttributeError:
            print 'Warning: Series/dataframe has no attribute \'{}\''.format(k)
    _pandas_add_meta_data(target, addmeta)
    return target


### METRICS FUNCTIONS

def get_amoc(ncfiles, latlim=(30,60), zlim=(500,9999), window_size=12):
    """Retrieve AMOC time series from a set of CESM/POP model output files

    Parameters
    ----------
    ncfiles : list of str
        paths to input files
    latlim : tuple (float,float)
        Latitude limits between which to find the maximum AMOC
    zlim : tuple (float,float)
        Depth limits between which to find the maximum AMOC
    window_size : int
        Smoothing window with to apply before taking maximum

    Returns
    -------
    Atlantic meridional overturning circulation (AMOC) maximum between 
    30N and 60N below 500 m water depth (following Shields et al. 2012).

    Note
    ----
    Only works with POP data that has the diagnostic variable 'MOC' included.

    """
    n = len(ncfiles)
    _nfiles_diag(n)
    
    maxn = get_ulimitn()

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
            timeax = _get_timeax(ds)
            amoc = dsvar['MOC'][:,1,0,kza:kzo+1,ja:jo+1]
    else:
        timeax = np.zeros(n,'object')
        amoc = np.zeros((n,nz,nlat))
        for i,fname in enumerate(ncfiles):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = _get_timeax(ds)[0]
                amoc[i] = dsvar['MOC'][0,1,0,kza:kzo+1,ja:jo+1]
                
    if window_size > 1:
        maxmeanamoc = np.max(np.max(scipy.ndimage.convolve1d(
            amoc,weights=np.ones(int(window_size))/float(window_size),
            axis=0,mode='constant',cval=np.nan),axis=-1),axis=-1)
        maxmeanamoc[:window_size+1] = np.nan
        maxmeanamoc[-window_size:] = np.nan
    else:
        maxmeanamoc = np.max(np.max(amoc,axis=-1),axis=-1)

    if use_pandas:
        index = pd.Index(datetime_to_decimal_year(timeax), name='ModelYear')
        ts = pd.Series(maxmeanamoc, index=index, name='AMOC')
        _pandas_add_meta_data(ts, meta=dict(
           latlim = latlim,
           zlim = zlim,
            ))
        return ts
    else:
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
    _nfiles_diag(n)
    maxn = get_ulimitn()

    with netCDF4.Dataset(ncfiles[0]) as ds:
        latax = ds.variables['lat_aux_grid'][:]
        ja = np.argmin(np.abs(latax-latlim[0]))
        jo = np.argmin(np.abs(latax-latlim[1]))
        nlat = jo-ja+1
        
    if n <= maxn:
        with netCDF4.MFDataset(ncfiles) as ds:
            dsvar = ds.variables
            timeax = _get_timeax(ds)
            nheat = dsvar['N_HEAT'][:,0,component,ja:jo+1]
    else:
        timeax = np.zeros(n,'object')
        nheat = np.zeros((n,nlat))
        for i,fname in enumerate(ncfiles):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = _get_timeax(ds)[0]
                nheat[i,:] = dsvar['N_HEAT'][0,0,component,ja:jo+1]
                
    window_size = 12
    maxmeannheat = np.max(scipy.ndimage.convolve1d(
        nheat,weights=np.ones(int(window_size))/float(window_size),
        axis=0,mode='constant',cval=np.nan),axis=-1)

    if use_pandas:
        index = pd.Index(datetime_to_decimal_year(timeax), name='ModelYear')
        ts = pd.Series(maxmeannheat, index=index, name='MHT')
        _pandas_add_meta_data(ts, meta=dict(
            latlim = latlim,
            component = component,
            ))
        return ts
    else:
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
    _nfiles_diag(n)
    maxn = get_ulimitn()

    with netCDF4.Dataset(ncfiles[0]) as ds:
        dsvar = ds.variables
        latax = dsvar['lat_aux_grid'][:]
        j0 = np.argmin(np.abs(latax-lat0))
        
    if n <= maxn:
        with netCDF4.MFDataset(ncfiles) as ds:
            dsvar = ds.variables
            timeax = _get_timeax(ds)
            nsalt = dsvar['N_SALT'][:,0,component,j0]
    else:
        timeax = np.zeros(n,'object')
        nsalt = np.zeros(n)
        for i,fname in enumerate(sorted(ncfiles)):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = _get_timeax(ds)[0]
                nsalt[i] = dsvar['N_SALT'][0,0,component,j0]
                
    window_size=12
    window = np.ones(int(window_size))/float(window_size)
    meannsalt = np.convolve(nsalt,window,'same')
    meannsalt[:window_size+1] = np.nan
    meannsalt[-window_size:] = np.nan

    if use_pandas:
        index = pd.Index(datetime_to_decimal_year(timeax), name='ModelYear')
        ts = pd.Series(meannsalt, index=index, name='MST')
        ts = _pandas_add_meta_data(ts, meta=dict(
            lat0 = lat0,
            component = component,
            ))
        return ts
    else:
        return meannsalt, timeax


def get_timeseries(ncfiles, varn, grid, 
        reducefunc=np.mean, 
        latlim=None, lonlim=None, k=0):
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
    k : int
        layer
    """
    n = len(ncfiles)
    _nfiles_diag(n)
    maxn = get_ulimitn()

    with netCDF4.Dataset(ncfiles[0]) as ds:
        dshape = ds.variables[varn].shape
        ndims = len(dshape)
        if ndims not in (3,4):
            raise IndexError('Fields with {} dimensions not supported.'.format(ndims))
        if latlim is None and lonlim is None:
            mask = None
        else:
            mask = poppygrid.get_mask_lonlat(ds, lonlim=lonlim, latlim=latlim, grid=grid)
            mask &= ds.variables['KM'+grid][:]>0
            jj,ii = np.where(mask)
        
    if n <= maxn:
        with netCDF4.MFDataset(ncfiles) as ds:
            dsvar = ds.variables
            timeax = _get_decimal_year(ds)
            if ndims == 4:
                if mask is None:
                    tseries = reducefunc(dsvar[varn][:,k,:,:], axis=(-1,-2))
                else:
                    tseries = reducefunc(dsvar[varn][:,k,jj,ii], axis=-1)
            else:
                if mask is None:
                    tseries = reducefunc(dsvar[varn][:,:,:], axis=(-1,-2))
                else:
                    tseries = reducefunc(dsvar[varn][:,jj,ii], axis=-1)
    else:
        print '... file by file ...'
        timeax = np.zeros(n, 'object')
        tseries = np.zeros((n))
        for i,fname in enumerate(ncfiles):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = _get_decimal_year(ds)[0]
                if ndims == 4:
                    if mask is None:
                        tseries[i] = reducefunc(dsvar[varn][0,k,:,:])
                    else:
                        tseries[i] = reducefunc(dsvar[varn][0,k,jj,ii])
                else:
                    if mask is None:
                        tseries[i] = reducefunc(dsvar[varn][0,:,:])
                    else:
                        tseries[i] = reducefunc(dsvar[varn][0,jj,ii])
            if np.mod(i,100) == 0:
                print '{}/{}'.format(i,n)
                
    if use_pandas:
        index = pd.Index(timeax, name='ModelYear')
        ts = pd.Series(tseries, index=index, name=varn)
        _pandas_add_meta_data(ts, meta=dict(
            latlim = latlim,
            lonlim = lonlim,
            varn = varn,
            reducefunc = str(reducefunc),
            k = k,
            grid = grid,
            ))
        return ts
    else:
        return tseries, timeax

