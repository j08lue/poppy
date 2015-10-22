"""POP grid handling"""
import numpy as np
import bisect
from scipy.interpolate import NearestNDInterpolator
import netCDF4


def find_k_depth(ds,depth):
    """In an open netCDF *ds*, find index of z-layer closest to the given depth"""
    z = ds.variables['z_w'][:]/100.
    kz = bisect.bisect_left(z,depth)
    return kz,z[kz]


def find_j_lat(ds,targetlat,varname='ULAT',region={},returnabs=True):
    """In an open netCDF *ds*, find index of latitude band closest to given *lat0*"""
    if region:
        i1,i2,j1,j2 = [region[key] for key in ['i1','i2','j1','j2']]
        meanlat = np.mean(
            ds.variables[varname][j1:j2,i1:i2],
            axis=1)
    else:
        meanlat = np.mean(
            ds.variables[varname][:],
            axis=1)

    jlat = np.argmin(np.abs(meanlat-targetlat))
    
    if returnabs and region:
        return jlat+j1,meanlat[jlat]
    else:
        return jlat,meanlat[jlat]


def get_regmasks(region_mask, fmt=bool):
    """Get boolean region masks based on the POP variable `REGION_MASK`
    
    Parameters
    ----------
    region_mask : ndarray
        values from dataset variable 'REGION_MASK'
    fmt : (int or bool)
        output type

    Note
    ----
    The actual values must be passed in region_mask, not the netCDF variable.
    """
    # initialize
    regions = ['Global', 'Atlantic', 'Pacific', 'Indo-Pacific']
    regmask = np.ones(region_mask.shape, dict(names=regions, formats=[fmt]*len(regions)))
    # regions
    regmask['Pacific'] = (region_mask == 2)
    regmask['Atlantic'] = (region_mask >= 6) & (region_mask != 7)
    regmask['Indo-Pacific'] = (region_mask >= 2) & (region_mask <= 3)
    return regmask


def get_mask_lonlat(fname,lonlim=None,latlim=None,grid='T'):
    """Mask the region confined by `lonlim` and `latlim` in the grid from `fname`
    
    Parameters
    ----------
    fname : str or open netCDF4.Dataset
        dataset to get the grid from 
    lonlim,latlim : tup
        limits for region
    grid : str, [T,U]
        which grid to use (T or U)
    """
    def _get_mask(ds,lonlim,latlim,grid):
        ny,nx = ds.variables[grid+'LONG'].shape
        mask = np.ones((ny,nx),bool)
        if lonlim is not None:
            lon = np.mod(ds.variables[grid+'LONG'][:],360)
            lonlim = np.mod(lonlim,360)
            if np.diff(lonlim) > 0:
                mask &= (lon >= lonlim[0]) & (lon <= lonlim[1])
            else:
                mask &= (lon >= lonlim[0]) | (lon <= lonlim[1])
        if latlim is not None:
            lat = ds.variables[grid+'LAT'][:]
            mask &= (lat >= latlim[0]) & (lat <= latlim[1])
        if not mask.any():
            raise ValueError('All masked.')
        return mask

    try:
        return _get_mask(fname,lonlim,latlim,grid)
    except AttributeError:
        with netCDF4.Dataset(fname) as ds:
            return _get_mask(ds,lonlim,latlim,grid)


def get_zintstr(zint,zaxlim):
    """Generate a descriptive string for *zint*
    depending on the depth axis ranges *zaxlim*"""
    zintstr = ''
    if zint and len(zint) == 2:
        if zint[0] > 0:
            if zint[1] > zaxlim[-1]:
                zintstr = 'below {:.0f}m'.format(zint[0])
            else:
                zintstr = '{0[0]:.0f}m-{0[1]:.0f}m'.format(zint)
        else:
            if zint[1] <= zaxlim[-1]:
                zintstr = 'above {:.0f}m'.format(zint[1])    
    return zintstr


def ix_rewrap_lon(lon,lon0,lat=None,latlim=()):
    """Locate the point in the 2D array `lon` where to rewrap it zonally 
    such that it starts with `lon0` in the west
    Assumes that `lon` has the shape (y,x). Transpose beforehand if not.

    Parameters
    ----------
    lon : 1d or 2d array (y,x)
        longitude grid
    lon0 : float
        longitude where to rewrap
    lat : 2d array, optional
        latitunde grid (y,x)
        required for latlim
    latlim : (float,float), optional
        latitude boundaries where to use values from lon
        useful when grid is highly irregular at higher latitudes
        in that case, set e.g. latlim = (-60,60)
        
    Returns
    -------
    indices for x and y axes of variables to rewrap the data
    """
    lon360 = np.mod(lon,360)
    lon0 = np.mod(lon0,360)

    if np.ndim(lon) == 2:
        if lat is not None and latlim:
            meanlon = np.mean(np.ma.masked_where(((lat>latlim[1]) | (lat<latlim[0])),lon360),axis=0)
        else:
            meanlon = np.mean(lon360,axis=0)
        i0 = np.argmin(np.abs(meanlon-lon0))
        ny,nx = lon.shape
        ii = np.concatenate([np.arange(i0,nx),np.arange(0,i0)])
        jj = np.arange(0,ny)
        return np.ix_(jj,ii)

    elif np.ndim(lon) == 1:
        i0 = np.argmin(np.abs(lon-lon0))
        nx = len(lon)
        return np.concatenate([np.arange(i0,nx),np.arange(0,i0)])

    else:
        raise ValueError('lon must be a 1D or 2D array!')



def find_region(longrid=None,latgrid=None,lonlim=(),latlim=(),latlim_include=()):
    """Find region in a near-regular lon/lat grid confined by lonlim and/or latlim

    Parameters
    ----------
    longrid, latgrid : 2D ndarrays
        Grids, must be 2D and with axes (lat,lon)
    lonlim, latlim : tuple
        Limits confining the region
    latlim_include : (float,float), optional
        latitude boundaries where to use values from lon
        useful when grid is highly irregular at higher latitudes
        in that case, set e.g. latlim_include = (-60,60)
    """
    if longrid is not None:
        ny, nx = longrid.shape
    elif latgrid is not None:
        ny, nx = latgrid.shape
    else:
        raise ValueError('At least one of longrid and latgrid must be provided.')
    
    if latgrid is not None and latlim_include:
        included_in_mean = ((latgrid>=latlim_include[0]) & (latgrid<=latlim_include[1]))
    else:
        included_in_mean = np.ones((ny,nx), dtype=bool)
    
    if longrid is not None and lonlim:
        nx = longrid.shape[-1]
        longrid = np.mod(longrid,360)
        lonlim = np.mod(lonlim,360)
        meanlons = np.mean(np.ma.masked_where(~included_in_mean,longrid),axis=0)
        ia = np.argmin(np.abs(meanlons-lonlim[0]))
        io = np.argmin(np.abs(meanlons-lonlim[1]))
        if ia > io:
            ii = np.concatenate([np.arange(ia,nx),np.arange(0,io+1)])
        else:
            ii = np.arange(ia,io+1)
    else:
        ii = np.arange(nx)

    if latgrid is not None and latlim:
        ny = latgrid.shape[0]
        meanlats = np.mean(latgrid[:,ii],axis=1)
        ja = np.argmin(np.abs(meanlats-latlim[0]))
        jo = np.argmin(np.abs(meanlats-latlim[1]))
        jj = np.arange(ja,jo+1)
    else:
        jj = np.arange(ny)
    
    return np.ix_(jj,ii)



def get_ndx(datafileref,datafile):
    """Get ratio of grid spacings between two grids"""

    with netCDF4.Dataset(datafileref) as dsref:

        lonref = dsref.variables['ULONG'][:]
        latref = dsref.variables['ULAT'][:]
        dxref = dsref.variables['DXU'][:]

        pointsref = np.vstack([latref.flat,lonref.flat]).T

        with netCDF4.Dataset(datafile) as ds:

            lon = ds.variables['ULONG'][:]
            lat = ds.variables['ULAT'][:]
            dx = ds.variables['DXU'][:]

            points = np.vstack([lat.flat,lon.flat]).T

            interpolator = NearestNDInterpolator(pointsref,dxref.flatten())

            dxref_int = np.reshape(interpolator(points),lon.shape)

            ndx = dxref_int/dx

    return ndx


def get_ugrid_from_file(fname,lon0_rewrap=None,latlim=(-60,60)):
    """Get grid from POP ocean model data netCDF file, 
    rewrap at *lon0_rewrap* using *latlim* (optional)"""
    with netCDF4.Dataset(fname) as ds:
        lon = ds.variables['ULONG'][:]
        lat = ds.variables['ULAT'][:]

    if lon0_rewrap is not None:
        ix = ix_rewrap_lon(lon,lon0_rewrap,lat,latlim=latlim)
        lon = lon[ix]
        lat = lat[ix]

    return lon,lat

