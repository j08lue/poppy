import netCDF4
import numpy as np
import calendar


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


def get_time_datetime(ds):
    """Get datetime instances from 'time' variable
    using netCDF4.num2date but accounting for reference year 0000
    
    Parameters
    ----------
    ds : open netCDF4.Dataset
        POP input data
    """
    dsvar = ds.variables
    time = dsvar['time']
    timedata = time[:]
    timeunits = time.units
    if timeunits.startswith('days since 0000'):
        timeunits = timeunits.replace('days since 0000', 'days since 0001')
        timedata -= 364
    return netCDF4.num2date(timedata, units=timeunits, calendar=time.calendar)


def get_time_decimal_year(ds, **kwarg):
    """Get time from 'time' variable in dataset and convert to decimal year
    
    Parameters
    ----------
    ds : open netCDF4.Dataset
        POP input data
    """
    return datetime_to_decimal_year(get_time_datetime(ds))
    


