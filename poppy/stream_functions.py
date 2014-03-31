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



def get_barotropic_stream_function(ds,region=None,lon0=None,t=0):
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

    if region is None:
        regmask = np.ones(lon.shape)
    else:
        regmask = poppy.grid.get_regmasks(dsvar['REGION_MASK'][:],int)[region]

    V = np.zeros(lat.shape)
    for k in xrange(len(dsvar['dz'])):
        V += (_fill0(dsvar['VVEL'][t,k]) * 1e-2
                * dsvar['dz'][k] * 1e-2)
    V *= dsvar['DXU'][:,:] * 1e-2
    V *= regmask

    if lon0 is not None:
        ix = poppy.grid.ix_rewrap_lon(lon,lon0,lat,latlim=(-60,60))
        lon = lon[ix]
        lat = lat[ix]
        V = V[ix]
        regmask = regmask[ix]

    # compute stream function:
    # cumulative zonal integral of meridional velocity
    psi = np.cumsum(V[:,::-1],axis=1)[:,::-1]
    psi *= -1.
    psi *= 1e-6 # convert to Sv
    psimasked = np.ma.masked_where(regmask==0,psi)

    return psimasked,lon,lat


def get_barotropic_stream_function_section(ds,ii,jj,region=None,t=0):
    """Get barotropic stream function section
    
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
    t = t or 0
    dsvar = ds.variables

    lon = dsvar['ULONG'][jj,ii]
    lat = dsvar['ULAT'][jj,ii]

    if region is None:
        regmask = np.ones(lon.shape)
    else:
        regmask = poppy.grid.get_regmasks(dsvar['REGION_MASK'][jj,ii],int)[region]

    if isinstance(jj,int):
        V = np.sum((_fill0(dsvar['VVEL'][t,:,jj,ii]) * 1e-2
            * dsvar['dz'][:][:,np.newaxis] * 1e-2),axis=0)
        V *= dsvar['DXU'][jj,ii] * 1e-2
        V *= regmask
        psi = np.cumsum(V[::-1])[::-1]
        psi *= -1.

    elif isinstance(ii,int):
        U = np.sum((_fill0(dsvar['UVEL'][t,:,jj,ii]) * 1e-2
            * dsvar['dz'][:][:,np.newaxis] * 1e-2),axis=0)
        U *= dsvar['DYU'][jj,ii] * 1e-2
        U *= regmask
        psi = np.cumsum(U[::-1])[::-1]

    else:
        raise ValueError('Either ii or jj must be a single integer for integration along a grid line.')

    # compute stream function
    psi *= 1e-6 # convert to Sv
    psimasked = np.ma.masked_where(regmask==0,psi)

    return psimasked,lon,lat
