"""Get AMOC evolution as a metric for model ocean state"""
import numpy as np
import netCDF4
import glob

import utils


def get_amoc(pattern,latlim=(30,60),zlim=(500,9999)):
    """Retrieve AMOC evolution from a set of CESM/POP model output files

    Parameters
    ----------
    pattern : str
        File search pattern used with glob.glob() to find the netCDF files to read
    latlim : tuple (float,float)
        Latitude limits between which to find the maximum AMOC
    zlim : tuple (float,float)
        Depth limits between which to find the maximum AMOC

    Returns
    -------
    Atlantic meridional overturning circulation (AMOC) maximum between 30N and 60N below 500 m water depth (following Shields et al. 2012).

    Note
    ----
    Only works with POP data that has the diagnostic variable 'MOC' included.

    """
    ncfiles = glob.glob(pattern)
    n = len(ncfiles)
    if n == 0: raise ValueError('No files found for pattern {}'.format(pattern))
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
        with netCDF4.MFDataset(pattern) as ds:
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
                
    return amoc,timeax

