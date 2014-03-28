import numpy as np

import poppy.grid

def _fill0(a):
    return np.ma.filled(a,0.)


def get_vertical_stream_function(ds,region='Global',t=0):
    """Get vertical stream function for a given region
    
    Parameters
    ----------
    ds : netCDF4.Dataset
        open netCDF dataset
    region : str
        region ID to be used with ``poppy.grid.get_regmasks``
    t : int
        time level (default: 0)

    Returns
    -------
    psi,lon,lat
    """
    dsvar = ds.variables

    lat = dsvar['ULAT'][:]
    dz = dsvar['dz'][:] * 1e-2
    dx = dsvar['DXU'][:,:] * 1e-2

    ny,nx = lat.shape
    nz = len(dz)

    if region is None or region == 'Global':
        bmask = np.ones((ny,nx),bool)
    else:
        bmask = poppy.grid.get_regmasks(dsvar['REGION_MASK'][:])[region]
    imask = np.asarray(bmask,dtype='i4')

    latax = lat[:,np.int(np.round(np.mean(np.where(bmask)[-1])))]
    latlim = (np.min(lat[bmask]),np.max(lat[bmask]))
    zax = dsvar['z_t'][:]*1e-2

    Vdz = np.zeros((nz,ny))
    for k in xrange(nz):
        Vdz[k,:] = np.sum((
                _fill0(dsvar['VVEL'][t,k,:,:]) * 1e-2
                * imask
                * dx
                ),axis=-1) * dz[k]

    # compute streamfunction in Sv
    psi = np.cumsum(Vdz,axis=0)
    psi *= 1e-6

    return zax,latax,latlim,psi



def get_barotropic_stream_function(ds,region,lon0=None,t=0):
    """Get barotropic stream function for a given region
    
    Parameters
    ----------
    ds : netCDF4.Dataset
        open netCDF dataset
    region : str
        region ID to be used with ``poppy.grid.get_regmasks``
    lon0 : float
        longitude at which to start the integration
    t : int
        time level (default: 0)

    Returns
    -------
    psi,lon,lat
    """
    dsvar = ds.variables

    lon = dsvar['ULONG'][:]
    lat = dsvar['ULAT'][:]

    regmask = poppy.grid.get_regmasks(dsvar['REGION_MASK'][:],int)[region]

    U = np.zeros(lat.shape)
    for k in xrange(len(dsvar['dz'])):
        U += (_fill0(dsvar['UVEL'][t,k]) * 1e-2
                * dsvar['dz'][k] * 1e-2)
    U *= dsvar['DYU'][:,:] * 1e-2
    U *= regmask

    if lon0 is not None:
        ix = poppy.grid.ix_rewrap_lon(lon,lon0,lat,latlim=(-60,60))
        lon = lon[ix]
        lat = lat[ix]
        U = U[ix]
        regmask = regmask[ix]

    # compute stream function
    psi = np.cumsum(U[::-1],axis=0)[::-1]
    psi *= 1e-6 # convert to Sv
    psimasked = np.ma.masked_where(regmask==0,psi)

    return psimasked,lon,lat
