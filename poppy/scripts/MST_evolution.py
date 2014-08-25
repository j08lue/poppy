#!/usr/bin/env python
"""Plot MHT from multi-file DS"""
import numpy as np
import netCDF4
import matplotlib.pyplot as plt ; plt.close('all')
import glob
import cPickle as pickle

import utils

def plot_mst_evolution(pattern,lat0=55,component=0,savefig=False,figname='',
        savetofile=None):
    """Plot MST evolution from CESM/POP data"""
    componentnames = {
        0 : 'Total',
        1 : 'Eulerian-Mean Advection',             
        2 : 'Eddy-Induced Advection (bolus) + Diffusion',
        3 : 'Eddy-Induced (bolus) Advection',
        4 : 'Submeso Advection',
        }

    ncfiles = glob.glob(pattern)
    n = len(ncfiles)
    print('Processing {} files ...'.format(n))

    maxn = utils.get_ulimitn()

    with netCDF4.Dataset(ncfiles[0]) as ds:
        dsvar = ds.variables
        latax = dsvar['lat_aux_grid'][:]
        j0 = np.argmin(np.abs(latax-lat0))
        
    if n <= maxn:
        with netCDF4.MFDataset(pattern) as ds:
            dsvar = ds.variables
            timeax = netCDF4.num2date(dsvar['time_bound'][:,0],
                                      dsvar['time_bound'].units,'proleptic_gregorian')
            nsalt = dsvar['N_SALT'][:,0,component,j0]
    else:
        timeax = np.zeros(n,'object')
        nsalt = np.zeros(n)
        for i,fname in enumerate(sorted(ncfiles)):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = netCDF4.num2date(dsvar['time_bound'][0,0],
                                             dsvar['time_bound'].units,'proleptic_gregorian')
                nsalt[i] = dsvar['N_SALT'][0,0,component,j0]
                
    window_size=12
    window = np.ones(int(window_size))/float(window_size)
    meannsalt = np.convolve(nsalt,window,'same')
    meannsalt[:window_size+1] = np.nan
    meannsalt[-window_size:] = np.nan

    if savetofile:
        dataout = dict(
                component=component,componentnames=componentnames,lat0=lat0,timeax=timeax,nsalt=nsalt,
                meannsalt=meannsalt,ncfiles=ncfiles)
        with open(savetofile,'wb') as fout:
            pickle.dump(dataout,fout)
    else:
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(timeax,meannsalt,'r')
        ax.grid(True)
        ax.set_ylabel('Meridional salt transport (Sv PPT)')
        ax.set_xlabel('integration year')
        ax.set_title('{} meridional salt transport at {}N'.format(componentnames[component],lat0))
        ax.xaxis_date()
        #fig.autofmt_xdate()

        if savefig or figname:
            figname = figname or 'MST_evolution.png'
            fig.savefig(figname,dpi=300)
        else:
            fig.show()


if __name__ == "__main__":

    defaults = dict(
        pattern = '*.pop.h.????-??.nc',
        component = 0,
        lat0 = 55,
        savefig = False,
        )

    import argparse
    import parser_functions as pf
    parser = argparse.ArgumentParser(description="Plot N_SALT evolution from CESM/POP data")
    parser.add_argument('-p','--pattern',type=pf.parse_pattern,
            help='file search pattern. Provide inside "" from command line!')
    parser.add_argument('-c','--component',type=int,choices=[0,1,2,3,4],help='transport component')
    parser.add_argument('-l','--lat0',type=pf.parse_coords,help='latitude where to extract salt transport')
    parser.add_argument('-s','--savefig',action='store_true',help='set this to save the figure with default filename')
    parser.add_argument('-f','--figname',type=str,help='save the figure to the given filename (implies -s)')
    parser.add_argument('--savetofile',type=str,help='save the data to this file (pickle)')
    
    parser.set_defaults(**defaults)
    args = parser.parse_args()
 
    plot_mst_evolution(**vars(args))
