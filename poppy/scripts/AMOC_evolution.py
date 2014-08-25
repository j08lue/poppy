#!/usr/bin/env python
"""Plot MOC from multi-file DS

MOC dimension key:
(time,transport_reg,moc_components,moc_z,lat_aux_grid)

transport_regions : (Global,Atlantic)
moc_components : (Eulerian Mean,Eddy-induced,Submeso)
moc_z : depth axis

"""
import numpy as np
import matplotlib.pyplot as plt ; plt.close('all')
import scipy.ndimage
import cPickle as pickle

from amoc import get_amoc

def plot_amoc_evolution(pattern,latlim=(30,60),zlim=(500,9999),
        window_size=12,savefig=False,figname='',
        savetofile=None):
    """Plot AMOC evolution from CESM/POP data"""

    amoc,timeax = get_amoc(pattern,latlim=latlim,zlim=zlim)

    if window_size > 1:
        maxmeanamoc = np.max(np.max(scipy.ndimage.convolve1d(
            amoc,weights=np.ones(int(window_size))/float(window_size),
            axis=0,mode='constant',cval=np.nan),axis=-1),axis=-1)
        maxmeanamoc[:window_size+1] = np.nan
        maxmeanamoc[-window_size:] = np.nan
    else:
        maxmeanamoc = np.max(np.max(amoc,axis=-1),axis=-1)

    if savetofile:
        dataout = dict(
                amoc=amoc,timeax=timeax,
                maxmeanamoc=maxmeanamoc,
                latlim=latlim,zlim=zlim,
                window_size=window_size)
        with open(savetofile,'wb') as fout:
            pickle.dump(dataout,fout)
    else:
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(timeax,maxmeanamoc,'b')
        ax.set_xlim(timeax.min(),timeax.max())
        ax.grid(True)
        ax.set_ylabel('AMOC (Sv)')
        ax.set_xlabel('integration year')
        ax.xaxis_date()

        if savefig or figname:
            figname = figname or 'AMOC_evolution.png'
            fig.savefig(figname,dpi=300)
        else:
            fig.show()


if __name__ == "__main__":

    defaults = dict(
        pattern = '*.pop.h.????-??.nc',
        latlim = (30,60),
        zlim = (500,9999),
        window_size = 12,
        savefig = False,
        )

    import argparse
    import parser_functions as pf
    parser = argparse.ArgumentParser(description="Plot AMOC evolution from CESM/POP data")
    parser.add_argument('-p','--pattern',type=pf.parse_pattern,
            help='file search pattern. Provide inside "" from command line!')
    parser.add_argument('-l','--latlim',type=pf.parse_coords,
            help='Latitude limits for AMOC region, e.g. 30,60. Takes negative values as e.g. m30 for -30.')
    parser.add_argument('-z','--zlim',type=lambda s: map(float,s.split(',')),
            help='Depth limits for AMOC region, e.g. 500,9999 for below 500 metres.')
    parser.add_argument('--window_size',type=int,
            help='window size for running mean (set to 0 to disable smoothing)')
    parser.add_argument('-s',dest='savefig',action='store_true',
            help='set this to save the figure with default filename')
    parser.add_argument('--figname',type=str,
            help='save the figure to the given filename (implies -s)')
    parser.add_argument('--savetofile',type=str,
            help='save the data to this file (pickle)')
    
    parser.set_defaults(**defaults)
    args = parser.parse_args()
 
    plot_amoc_evolution(**vars(args))
