import numpy as np
from oceanpy.stats import central_differences

def _fill0(a):
    return np.ma.filled(a,0.)

def mean_velocity_component(ds,varn,regmask=1,kza=0,kzo=None,S0=34.8):
    """Mean velocity component using VNT or VNS

    From https://bb.cgd.ucar.edu/node/1000983 :
    UET is obtained using
    
        UET(i,j) = 0.5*( u(i,j)*dyu(i,j) + u(i,j-1)*dyu(i,j-1) ) * 0.5*(T(i,j) + T(i+1,j)) # units m2 s-1 degC
        UET(i,j) = UET(i,j)/TAREA(i,j) # result has units degC s-1
    
    This means that UET has to be multiplied by DXU (and VNT by DYU) to get units [degC m s-1]
    """
    dsvar = ds.variables
    dxu = dsvar['DXU'][:]/100.
    dyu = dsvar['DYU'][:]/100.
    dz = dsvar['dz'][:]/100.
    if kzo is None: kzo = len(dz)
    meanvel = np.zeros(regmask.shape[0])
    for k in xrange(kza,kzo):
        if varn == 'heat':
            layer = dsvar['VNT'][0,k] # degC s-1
            layer *= dyu
        elif varn == 'salt':
            layer = dsvar['VNS'][0,k] # PPT s-1
            layer *= dyu
        elif varn == 'freshwater':
            layer = (S0 - dsvar['SALT'][0,k]) / S0
            layer *= dsvar['VVEL'][0,k]/100.
        layer *= dz[k]
        layer *= dxu
        layer *= regmask
        meanvel += np.sum(layer,axis=-1)
    if varn == 'heat':
        meanvel *= (1e3 * 4e3 * 1e-15)
    elif varn == 'salt':
        meanvel *= 1e-6 # Sv PPT
    return meanvel


def diffusion_component(ds,varn,regmask=1,kza=0,kzo=None,S0=34.8):
    """Temperature/Salt diffusion"""
    dsvar = ds.variables
    dxt = dsvar['DXT'][:]/100.
    dyt = dsvar['DYT'][:]/100.
    dz = dsvar['dz'][:]/100.
    if kzo is None: kzo = len(dz)
    diffusion = np.zeros(regmask.shape[0])
    for k in xrange(kza,kzo):
        layer = _fill0(dsvar['KAPPA_ISOP'][0,k]/1e4) # m2 s-1
        if varn == 'heat':
            scalar = _fill0(dsvar['TEMP'][0,k])
        elif varn == 'salt':
            scalar = _fill0(dsvar['SALT'][0,k])
        elif varn == 'freshwater':
            scalar = (S0 - _fill0(dsvar['SALT'][0,k])) / S0
        gradient = central_differences(scalar,dyt,axis=0) # [scalar] m s-1
        #gradient = np.zeros(scalar.shape)
        #gradient[1:,:] = np.diff(scalar,axis=0)
        #gradient[0] = gradient[1]
        #gradient /= dyt
        layer *= gradient
        layer *= dz[k]
        layer *= dxt
        layer *= regmask
        diffusion += np.sum(layer,axis=-1)
        diffusion *= -1.
    if varn == 'heat':
        diffusion *= (1e3 * 4e3 * 1e-15) # PW
    elif varn == 'salt':
        diffusion *= 1e-6 # Sv PPT
    return diffusion


def bolus_velocity_component_vnt_isop(ds,varn,regmask=1,kza=0,kzo=None,S0=0):
    """Eddy-induced velocity / bolus velocity using VNT_ISOP"""
    dsvar = ds.variables
    dxu = dsvar['DXU'][:]/100.
    dyu = dsvar['DYU'][:]/100.
    dz = dsvar['dz'][:]/100.
    if kzo is None: kzo = len(dz)
    bolus = np.zeros(regmask.shape[0])
    for k in xrange(kza,kzo):
        if varn == 'heat':
            layer = dsvar['VNT_ISOP'][0,k] # degC s-1
        elif varn == 'salt':
            layer = dsvar['VNS_ISOP'][0,k] # PPT s-1
        elif varn == 'freshwater':
            raise NotImplementedError('Salinity normalization does not work with this function.\n \
                    Use `_bolus_velocity_component_visop` instead.')
        layer *= dyu
        layer *= dz[k]
        layer *= dxu
        layer *= regmask
        bolus += np.sum(layer,axis=-1)
    if varn == 'heat':
        bolus *= (1e3 * 4e3 * 1e-15) # convert [degC m3 s-1] to [PW]
    elif varn == 'varn':
        bolus *= 1e-6 # Sv PPT
    return bolus


def bolus_velocity_component_visop(ds,varn,regmask=1,kza=0,kzo=None,S0=0):
    """Eddy-induced velocity / bolus velocity using VISOP variable"""
    dsvar = ds.variables
    dxu = dsvar['DXU'][:]/100.
    dz = dsvar['dz'][:]/100.
    if kzo is None: kzo = len(dz)
    bolus = np.zeros(regmask.shape[0])
    for k in xrange(kza,kzo):
        layer = dsvar['VISOP'][0,k]/100.
        if varn == 'heat':
            layer *= dsvar['TEMP'][0,k]
        elif varn == 'salt':
            layer *= dsvar['SALT'][0,k]
        elif varn == 'freshwater':
            layer = (S0 - _fill0(dsvar['SALT'][0,k])) / S0
        layer *= dxu
        layer *= dz[k]
        layer *= regmask
        bolus += np.sum(layer,axis=-1)
    if varn == 'heat':
        bolus *= (1e3 * 4e3 * 1e-15) # convert [degC m3 s-1] to [PW]
    elif varn == 'salt':
        bolus *= 1e-6 # Sv PPT
    return bolus

bolus_velocity_component = bolus_velocity_component_visop

