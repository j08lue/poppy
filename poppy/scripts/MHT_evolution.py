#!/usr/bin/env python
"""Plot MHT from multi-file DS"""
import numpy as np
import netCDF4
import matplotlib.pyplot as plt ; plt.close('all')
import scipy.ndimage
import glob
import cPickle as pickle

import utils


def plot_mht_evolution(pattern,latlim=(30,60),component=0,savefig=False,figname='',
        savetofile=None):
    """Plot MHT evolution from CESM/POP data"""
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
        latax = ds.variables['lat_aux_grid'][:]
        ja = np.argmin(np.abs(latax-latlim[0]))
        jo = np.argmin(np.abs(latax-latlim[1]))
        nlat = jo-ja+1
        
    if n <= maxn:
        with netCDF4.MFDataset(pattern) as ds:
            dsvar = ds.variables
            timeax = netCDF4.num2date(dsvar['time_bound'][:,0],
                                      dsvar['time_bound'].units,'proleptic_gregorian')
            nheat = dsvar['N_HEAT'][:,0,component,ja:jo+1]
    else:
        timeax = np.zeros(n,'object')
        nheat = np.zeros((n,nlat))
        for i,fname in enumerate(sorted(ncfiles)):
            with netCDF4.Dataset(fname) as ds:
                dsvar = ds.variables
                timeax[i] = netCDF4.num2date(dsvar['time_bound'][0,0],
                                             dsvar['time_bound'].units,'proleptic_gregorian')
                nheat[i,:] = dsvar['N_HEAT'][0,0,component,ja:jo+1]
                
    window_size = 12
    maxmeannheat = np.max(scipy.ndimage.convolve1d(
        nheat,weights=np.ones(int(window_size))/float(window_size),
        axis=0,mode='constant',cval=np.nan),axis=-1)

    if savetofile:
        dataout = dict(
                component=component,componentnames=componentnames,latlim=latlim,timeax=timeax,nheat=nheat,
                maxmeannheat=maxmeannheat,ncfiles=ncfiles)
        with open(savetofile,'wb') as fout:
            pickle.dump(dataout,fout)
    else:
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(timeax,maxmeannheat,'g')
        ax.grid(True)
        ax.set_ylabel('Meridional heat transport (PW)')
        ax.set_xlabel('integration year')
        ax.set_title('{0} MHT maximum between {1[0]}N and {1[1]}N'.format(componentnames[component],latlim))
        ax.xaxis_date()
        #fig.autofmt_xdate()
        
        if savefig or figname:
            figname = figname or 'MHT_evolution.png'
            fig.savefig(figname,dpi=300)
        else:
            fig.show()


if __name__ == "__main__":

    defaults = dict(
        pattern = '*.pop.h.????-??.nc',
        component = 0,
        latlim=(30,60),
        savefig = False,
        )

    import argparse
    import parser_functions as pf
    parser = argparse.ArgumentParser(description="Plot N_HEAT evolution from CESM/POP data")
    parser.add_argument('-p','--pattern',type=pf.parse_pattern,
            help='file search pattern. Provide inside "" from command line!')
    parser.add_argument('-c','--component',type=int,choices=[0,1,2,3,4],help='transport component')
    parser.add_argument('-s','--savefig',action='store_true',help='set this to save the figure with default filename')
    parser.add_argument('-l','--latlim',type=pf.parse_coords,help='latitude limits where to search for max mht')
    parser.add_argument('--figname',type=str,help='save the figure to the given filename (implies -s)')
    parser.add_argument('--savetofile',type=str,help='save the data to this file (pickle)')
    
    parser.set_defaults(**defaults)
    args = parser.parse_args()
 
    plot_mht_evolution(**vars(args))
