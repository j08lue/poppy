import numpy as np
from oceanpy.stats import central_differences


def mean_velocity_component(ds,varn,regmask=1,kza=0,kzo=None,S0=34.8):
    """Mean velocity component using VNT or VNS

    From https://bb.cgd.ucar.edu/node/1000983 :
    UET is obtained using
    
        UET(i,j) = 0.5*( u(i,j)*dyu(i,j) + u(i,j-1)*dyu(i,j-1) ) * 0.5*(T(i,j) + T(i+1,j)) # units m2 s-1 degC
        UET(i,j) = UET(i,j)/TAREA(i,j) # result has units degC s-1
    
    This means that UET has to be multiplied by DXU (and VNT by DYU) to get units [degC m s-1]
    """
    dxu = ds.variables['DXU'][:]/100.
    dyu = ds.variables['DYU'][:]/100.
    dz = ds.variables['dz'][:]/100.
    if kzo is None: kzo = len(dz)
    meanvel = np.zeros(regmask.shape[0])
    for k in xrange(kza,kzo):
        if varn == 'heat':
            layer = ds.variables['VNT'][0,k] # degC s-1
            layer *= dyu
        elif varn == 'salt':
            layer = ds.variables['VNS'][0,k] # PPT s-1
            layer *= dyu
        elif varn == 'freshwater':
            layer = (S0 - ds.variables['SALT'][0,k]) / S0
            layer *= ds.variables['VVEL'][0,k]/100.
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
    dxt = ds.variables['DXT'][:]/100.
    dyt = ds.variables['DYT'][:]/100.
    dz = ds.variables['dz'][:]/100.
    if kzo is None: kzo = len(dz)
    diffusion = np.zeros(regmask.shape[0])
    for k in xrange(kza,kzo):
        layer = np.ma.filled(ds.variables['KAPPA_ISOP'][0,k]/1e4,0) # m2 s-1
        if varn == 'heat':
            scalar = np.ma.filled(ds.variables['TEMP'][0,k],0)
        elif varn == 'salt':
            scalar = np.ma.filled(ds.variables['SALT'][0,k],0)
        elif varn == 'freshwater':
            scalar = (S0 - np.ma.filled(ds.variables['SALT'][0,k],0)) / S0
        layer *= central_differences(scalar,dyt,axis=0) # [scalar] m s-1
        layer *= dz[k]
        layer *= dxt
        layer *= regmask
        diffusion += np.sum(layer,axis=-1)
    if varn == 'heat':
        diffusion *= (1e3 * 4e3 * 1e-15) # PW
    elif varn == 'salt':
        diffusion *= 1e-6 # Sv PPT
    return diffusion


def bolus_velocity_component_vnt_isop(ds,varn,regmask=1,kza=0,kzo=None,S0=0):
    """Eddy-induced velocity / bolus velocity using VNT_ISOP"""
    dxu = ds.variables['DXU'][:]/100.
    dyu = ds.variables['DYU'][:]/100.
    dz = ds.variables['dz'][:]/100.
    if kzo is None: kzo = len(dz)
    bolus = np.zeros(regmask.shape[0])
    for k in xrange(kza,kzo):
        if varn == 'heat':
            layer = ds.variables['VNT_ISOP'][0,k] # degC s-1
        elif varn == 'salt':
            layer = ds.variables['VNS_ISOP'][0,k] # PPT s-1
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
    dxu = ds.variables['DXU'][:]/100.
    dz = ds.variables['dz'][:]/100.
    if kzo is None: kzo = len(dz)
    bolus = np.zeros(regmask.shape[0])
    for k in xrange(kza,kzo):
        layer = ds.variables['VISOP'][0,k]/100.
        if varn == 'heat':
            layer *= ds.variables['TEMP'][0,k]
        elif varn == 'salt':
            layer *= ds.variables['SALT'][0,k]
        elif varn == 'freshwater':
            layer = (S0 - np.ma.filled(ds.variables['SALT'][0,k],0)) / S0
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

