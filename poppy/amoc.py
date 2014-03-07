"""Get AMOC evolution as a metric for model ocean state"""
import numpy as np

def get_max_amoc(ds,zmin=500,latlim=(30,60),t=0):
    """Retrieve AMOC from a CESM/POP model output file

    Parameters
    ----------
    ds : open netCDF4.Dataset
        Dataset
    zmin : float
        minimum depth below which to find the maximum
    latlim : (float,float)
        latitude bounds between which to find the maximum
    t : int
        time level (default 0)

    Returns
    -------
    Atlantic meridional overturning circulation (AMOC) maximum between 
    30N and 60N below 500 m water depth (default, following Shields et al. 2012).

    Note
    ----
    Only works with POP data that has the diagnostic variable 'MOC' included.

    """
    zax = ds.variables['moc_z'][:]/100.
    kmin = np.argmin(np.abs(zax-500))

    latax = ds.variables['lat_aux_grid'][:]
    ja = np.argmin(np.abs(latax-latlim[0]))
    jo = np.argmin(np.abs(latax-latlim[1]))        
    
    return np.max(ds.variables['MOC'][t,1,0,kmin:,ja:jo+1])
