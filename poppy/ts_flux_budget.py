import numpy as np
from oceanpy.fluxbudget import budget_over_region_2D
from oceanpy.stats import central_differences

from poppy.surface_salt_fluxes import salt_flux_to_fw_flux

def _fill0(a):
    return np.ma.filled(a,0.)


def fluxbudget_VVEL(ds,mask,varn,kza=0,kzo=None,S0=34.8,t=0):
    """Integrate horizontal flux using VVEL*SCALAR"""
    dsvar = ds.variables
    dxu = dsvar['DXU'][:] * 1e-2
    dyu = dsvar['DYU'][:] * 1e-2
    dz = dsvar['dz'][:] * 1e-2
    if kzo is None: kzo = len(dz)
    fluxbudget = 0.
    for k in xrange(kza,kzo):
        uflux = _fill0(dsvar['UVEL'][t,k]) * 1e-2
        uflux *= dyu
        uflux *= dz[k]
        vflux = _fill0(dsvar['VVEL'][t,k]) * 1e-2
        vflux *= dxu
        vflux *= dz[k]
        if not varn:
            scalar = None
        elif varn == 'heat':
            scalar = _fill0(dsvar['TEMP'][t,k])
        elif varn == 'salt':
            scalar = _fill0(dsvar['SALT'][t,k])
        elif varn == 'freshwater':
            scalar = (S0 - _fill0(dsvar['SALT'][t,k])) / S0
        fluxbudget += budget_over_region_2D(uflux,vflux,scalar=scalar,mask=mask,grid='ArakawaB')
    if varn == 'heat':
        fluxbudget *= (1e3 * 4e3 * 1e-15) # PW
    elif varn == 'salt':
        fluxbudget *= salt_flux_to_fw_flux(ds)
    return fluxbudget


def fluxbudget_UESVNS(ds,mask,kza=0,kzo=None,t=0):
    """Integrate horizontal flux using UES and VNS variables"""
    dsvar = ds.variables
    dz = dsvar['dz'][:] * 1e-2
    tarea = dsvar['UAREA'][:] * 1e-4
    if kzo is None: kzo = len(dz)
    fluxbudget = 0.
    for k in xrange(kza,kzo):
        uflux = _fill0(dsvar['UES'][t,k])
        uflux *= tarea
        uflux *= dz[k]
        vflux = _fill0(dsvar['VNS'][t,k])
        vflux *= tarea
        vflux *= dz[k]
        fluxbudget += budget_over_region_2D(uflux,vflux,scalar=None,mask=mask)
    fluxbudget *= salt_flux_to_fw_flux(ds)
    return fluxbudget


def fluxbudget_bolus_visop(ds,mask,varn,kza=0,kzo=None,S0=34.8,t=0):
    """Compute flux of `varn` into region `mask` due to eddy (bolus) velocity"""
    dsvar = ds.variables
    dxt = dsvar['DXT'][:] * 1e-2
    dyt = dsvar['DYT'][:] * 1e-2
    dz = dsvar['dz'][:] * 1e-2
    if kzo is None: kzo = len(dz)
    fluxbudget = 0.
    for k in xrange(kza,kzo):
        # get bolus velocity
        uflux = _fill0(dsvar['UISOP'][t,k]) * 1e-2 # m s-1
        vflux = _fill0(dsvar['VISOP'][t,k]) * 1e-2 # m s-1
        # get scalar data
        if varn == 'heat':
            scalar = _fill0(dsvar['TEMP'][t,k])
        elif varn == 'salt':
            scalar = _fill0(dsvar['SALT'][t,k])
        elif varn == 'freshwater':
            scalar = (S0 - _fill0(dsvar['SALT'][t,k])) / S0
        # multiply flux by scalar
        uflux *= scalar
        vflux *= scalar
        # multiply by horizontal grid spacing
        uflux *= dyt
        vflux *= dxt
        # multiply by vertical grid spacing
        uflux *= dz[k]
        vflux *= dz[k]
        # compute budget
        fluxbudget += budget_over_region_2D(uflux,vflux,scalar=None,mask=mask)
    if varn == 'heat':
        fluxbudget *= (1e3 * 4e3 * 1e-15) # PW
    elif varn == 'salt':
        fluxbudget *= salt_flux_to_fw_flux(ds)
    return fluxbudget

fluxbudget_bolus = fluxbudget_bolus_visop


def fluxbudget_diffusion(ds,mask,varn,kza=0,kzo=None,S0=34.8,t=0):
    """Compute flux of `varn` into region `mask` due to diffusion"""
    dsvar = ds.variables
    dxt = dsvar['DXT'][:] * 1e-2
    dyt = dsvar['DYT'][:] * 1e-2
    dz = dsvar['dz'][:] * 1e-2
    if kzo is None: kzo = len(dz)
    fluxbudget = 0.
    for k in xrange(kza,kzo):
        # get scalar data
        if varn == 'heat':
            scalar = _fill0(dsvar['TEMP'][t,k])
        elif varn == 'salt':
            scalar = _fill0(dsvar['SALT'][t,k])
        elif varn == 'freshwater':
            scalar = (S0 - _fill0(dsvar['SALT'][t,k])) / S0
        # get gradient
        uflux = central_differences(scalar,dxt,axis=1) # [scalar] m-1
        vflux = central_differences(scalar,dyt,axis=0) # [scalar] m-1
        # multiply gradient by diffusion coefficient
        kappa = _fill0(dsvar['KAPPA_ISOP'][t,k] * 1e-4) # m2 s-1
        uflux *= kappa
        vflux *= kappa
        # multiply by horizontal grid spacing
        uflux *= dyt
        vflux *= dxt
        # multiply by vertical grid spacing
        uflux *= dz[k]
        vflux *= dz[k]
        # compute budget
        fluxbudget += budget_over_region_2D(uflux,vflux,scalar=None,mask=mask)
    # convert to right units
    if varn == 'heat':
        fluxbudget *= (1e3 * 4e3 * 1e-15) # PW
    elif varn == 'salt':
        fluxbudget *= salt_flux_to_fw_flux(ds)
    return fluxbudget


def fluxbudget_bolus_advection_tendency(ds,mask,varn,t=0):
    dsvar = ds.variables
    if varn == 'heat':
        integrand = _fill0(dsvar['ADVT_ISOP'][t][mask]) * 1e-2
    elif varn == 'salt':
        integrand = _fill0(dsvar['ADVS_ISOP'][t][mask]) * 1e-2
    else:
        raise ValueError('This function only works for heat and salt transport.')
    integrand *= dsvar['TAREA'][:][mask] * 1e-4
    integral = np.sum(integrand)
    if varn == 'salt':
        integral *= salt_flux_to_fw_flux(ds)
    return integral


def transport_divergence(ds,mask,kza=0,kzo=None,t=0):
    dsvar = ds.variables
    dxu = dsvar['DXU'][:] * 1e-2
    dyu = dsvar['DYU'][:] * 1e-2
    tarea = dsvar['TAREA'][:] * 1e-4
    dz = dsvar['dz'][:] * 1e-2    
    if kzo is None: kzo = len(dz)
    transport_divergence = 0.
    for k in xrange(kza,kzo):
        uflux = _fill0(dsvar['UES'][t,k])
        uflux *= dyu
        uflux *= dz[k]
        uflux *= mask
        vflux = _fill0(dsvar['VNS'][t,k])
        vflux *= dxu
        vflux *= dz[k]
        vflux *= mask
        divergence = central_differences(uflux,dxu,axis=1) + central_differences(vflux,dyu,axis=0)
        divergence *= mask
        transport_divergence += np.sum(divergence*tarea)
    transport_divergence *= salt_flux_to_fw_flux(ds)
    return transport_divergence


def transport_divergence_from_vertical(ds,mask,kza=0,kzo=None,t=0):
    dsvar = ds.variables
    dz = dsvar['dz'][:] * 1e-2
    if kzo is None: kzo = len(dz)
    transport_divergence = 0.
    for k in xrange(kza,kzo):
        wflux = _fill0(dsvar['WTS'][t,k][mask])
        wflux *= dz[k]
        wflux *= dsvar['TAREA'][:][mask] * 1e-4
        transport_divergence += np.sum(wflux)
    transport_divergence *= salt_flux_to_fw_flux(ds)
    return transport_divergence
