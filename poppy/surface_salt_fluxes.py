"""
The net surface forcing of POP model salinity is
SFWF - QFLUX/latent_heat_fusion/1.e4            [kg FW /m^2/s]

where bold-faced indicates variables found in POP model output netcdf, latent heat of fusion can also
be found there, and the 1.e4 factor is needed for correct units. The QFLUX term isn't included in SFWF
because it's not considered a flux exchanged between model components--it is the FW flux associated with
frazil ice formation which occurs within the ocean model!


This net surface flux can be broken into components as follows:
SFWF - QFLUX/latent_heat_fusion/1.e4   =  
    PREC_F + EVAP_F + ROFF_F + IOFF_F + MELT_F
    + SALT_F*(sflux_factor/salinity_factor)
    - QFLUX/latent_heat_fusion/1.e4
    + (any weak or strong salinity restoring)


MELT_F represents any FW flux computed by the ice model and sent to the ocean (associated with ice melt/growth), while SALT_F represents any SALT flux computed by the ice model and sent to the ocean (associated with the fact that ice melt/growth requires that salt be added/removed from the ocean because ice has a constant salinity of ~4 psu!). This salt flux must of course be converted to appropriate units; in the above equation, it is converted from (kg SALT/m^2/s) to (kg FW/m^2/s). 

Note that we don't generally save salinity restoring fluxes, so any difference between the RHS and LHS of the above equation computed from POP output fields might be attributable to salinity restoring.

SFWF - (PREC_F + EVAP+F + ROFF_F + IOFF_F + MELT_F) - SALT_F*(sflux_factor/salinity_factor)
= (any weak or strong salinity restoring)


The equation above can be compared to the POP code in forcing_coupled.F90; specifically the line
        STF(:,:,2,iblock) = RCALCT(:,:,iblock)*(  &
            (PREC_F(:,:,iblock) + EVAP_F(:,:,iblock) +  &
                MELT_F(:,:,iblock) + ROFF_F(:,:,iblock) + IOFF_F(:,:,iblock))*salinity_factor &
            + SALT_F(:,:,iblock)*sflux_factor)


For guidance on conversion of FW flux to virtual salt flux, etc, refer to POP_ConstantsMod.F90
"""
import netCDF4


def net_salinity_forcing(ncfile):
    """
    The net surface forcing of POP model salinity is
        SFWF - QFLUX/latent_heat_fusion/1.e4  [kg FW /m^2/s]
    """
    def _get_data(ds):
        dsv = ds.variables
        return dsv['SFWF'][0] - dsv['QFLUX'][0]/dsv['latent_heat_fusion'][:]/1.e4 # units [kg FW /m^2/s]
    try:
        return _get_data(ncfile)
    except AttributeError:
        with netCDF4.Dataset(ncfile) as ds:
            return _get_data(ds)


def salinity_restoring(ncfile):
    """
    Note that we don't generally save salinity restoring fluxes, 
    so any difference between the RHS and LHS of the above equation 
    computed from POP output fields might be attributable to salinity restoring.

        SFWF - (PREC_F + EVAP+F + ROFF_F + IOFF_F + MELT_F) - SALT_F*(sflux_factor/salinity_factor)
            = (any weak or strong salinity restoring)
    """
    def _get_data(ds):
        dsv = ds.variables
        data = dsv['SFWF'][0]
        data -= (dsv['PREC_F'][0] + dsv['EVAP_F'][0] + dsv['ROFF_F'][0] + dsv['IOFF_F'][0] + dsv['MELT_F'][0]) 
        data -= dsv['SALT_F'][0]*(dsv['sflux_factor'][:]/dsv['salinity_factor'][:])
        return data
    try:
        return _get_data(ncfile)
    except AttributeError:
        with netCDF4.Dataset(ncfile) as ds:
            return _get_data(ds)


def salt_flux_to_fw_flux(ds):
    """Returns the factor with which to multiply salt flux with to get fresh water flux
    i.e. (kg SALT/m^2/s) to (kg FW/m^2/s)
    """
    dsvar = ds.variables
    return dsvar['sflux_factor'][0] / dsvar['salinity_factor'][0]
