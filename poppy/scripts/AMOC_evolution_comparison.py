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
import os
import itertools

from amoc import get_amoc


def plot_amoc_evolution(patterns,labels=[],defaultpattern='*.pop.h.????-??.nc',window_size=12,savefig=False,figname=''):
    """Plot AMOC evolution from CESM/POP data"""

    if defaultpattern:
        for i in xrange(len(patterns)):
            if os.path.isdir(patterns[i]):
               patterns[i] = os.path.join(patterns[i],defaultpattern)
    
    fig = plt.figure()
    ax = fig.gca()
    lscycler = itertools.cycle(["-","--","-.",":"])
    
    for i,pattern in enumerate(patterns):
        amoc,timeax = get_amoc(pattern)

        maxmeanamoc = np.max(np.max(scipy.ndimage.convolve1d(
            amoc,weights=np.ones(int(window_size))/float(window_size),
            axis=0,mode='constant',cval=np.nan),axis=-1),axis=-1)
        maxmeanamoc[:window_size+1] = np.nan
        maxmeanamoc[-window_size:] = np.nan

        if len(labels) == len(patterns):
            label = labels[i]
        else:
            label = 'pattern{:01d}'.format(i)
        ax.plot(timeax,maxmeanamoc,ls=next(lscycler),label=label)

    if len(patterns) > 1:
        ax.legend(loc=0)
    ax.autoscale(axis='x',tight=True)
    ax.grid(True)
    ax.set_ylabel('AMOC (Sv)')
    ax.set_xlabel('integration year')
    ax.xaxis_date()
    #fig.autofmt_xdate()

    if savefig or figname:
        figname = figname or 'AMOC_evolution_comparison.png'
        fig.savefig(figname,dpi=300)
    else:
        fig.show()



if __name__ == "__main__":

    defaults = dict(
        patterns = ['.'],
        labels = [],
        defaultpattern = '*.pop.h.????-??.nc',
        window_size = 12,
        savefig = False,
        )

    import argparse
    parser = argparse.ArgumentParser(description="Plot AMOC evolution from CESM/POP data")
    parser.add_argument('--patterns',type=(lambda s: s.split(',')),help='file search patterns, separated by comma. Provide inside "" from command line')
    parser.add_argument('--labels',type=(lambda s: s.split(',')),help='labels for lines corresponding to patterns')
    parser.add_argument('--defaultpattern',type=str,help='extend directories with this pattern')
    parser.add_argument('--window_size',type=int,help='window size for running mean (set to 0 to disable)')
    parser.add_argument('-s',dest='savefig',action='store_true',help='set this to save the figure with default filename')
    parser.add_argument('--figname',type=str,help='save the figure to the given filename (implies -s)')
    
    parser.set_defaults(**defaults)
    args = parser.parse_args()
 
    plot_amoc_evolution(**vars(args))
