#!/usr/bin/env python

import matplotlib.pyplot as plt ; plt.close('all')
plt.style.use('ggplot')
import glob

import poppy.metrics 


def plot_amoc_time_series(pattern, latlim=(30,60), zlim=(500,9999), savefig=False, figname=''):
    """Plot AMOC time series from CESM/POP data"""
    ncfiles = sorted(glob.glob(pattern))
    ts = poppy.metrics.get_amoc(ncfiles, latlim=latlim, zlim=zlim)
    ts.plot()
    plt.ylabel('AMOC (Sv)')
    plt.xlabel('integration year')

    if savefig or figname:
        figname = figname or 'AMOC_time_series.png'
        plt.gca().figure.savefig(figname,dpi=300)
    else:
        plt.show()


if __name__ == "__main__":

    def parse_coords(s):
        result = map(float,(s.replace('m','-')).split(','))
        return result[0] if len(result) == 0 else result

    import argparse
    parser = argparse.ArgumentParser(description="Plot AMOC time series from CESM/POP data")
    parser.add_argument('pattern', nargs='+',
            help='file search pattern. Provide inside "" from command line!',
            default='*.pop.h.????-??.nc')
    parser.add_argument('-l','--latlim',type=parse_coords,
            help='Latitude limits for AMOC region, e.g. 30,60. Takes negative values as e.g. m30 for -30.', 
            default=(30,60))
    parser.add_argument('-z','--zlim',type=lambda s: map(float,s.split(',')),
            help='Depth limits for AMOC region, e.g. 500,9999 for below 500 metres.',
            default=(500,9999))
    parser.add_argument('-s',dest='savefig',action='store_true',
            help='set this to save the figure with default filename')
    parser.add_argument('--figname',type=str,
            help='save the figure to the given filename (implies -s)')
    
    args = parser.parse_args()

    if len(args.files) == 1:
        args.files = sorted(glob.glob(args.files[0]))
 
    plot_amoc_time_series(vars(args))
 
