#!/usr/bin/env python

import matplotlib.pyplot as plt ; plt.close('all')
import cPickle as pickle
import glob

import poppy.metrics 
poppy.metrics.use_pandas = False


def save_amoc_time_series(pattern, outfname, latlim=(30,60), zlim=(500,9999)):
    """Plot AMOC time series from CESM/POP data"""
    ncfiles = sorted(glob.glob(pattern))
    maxmeanamoc,timeax = poppy.metrics.get_amoc(ncfiles, latlim=latlim, zlim=zlim)
    dataout = dict(
            timeax=timeax,
            maxmeanamoc=maxmeanamoc,
            latlim=latlim,zlim=zlim)
    with open(outfname,'wb') as fout:
        pickle.dump(dataout,fout)


def plot_amoc_time_series(pattern, latlim=(30,60), zlim=(500,9999), savefig=False, figname=''):
    """Plot AMOC time series from CESM/POP data"""
    ncfiles = sorted(glob.glob(pattern))
    maxmeanamoc,timeax = poppy.metrics.get_amoc(ncfiles, latlim=latlim, zlim=zlim)
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(timeax,maxmeanamoc,'b')
    ax.set_xlim(timeax.min(),timeax.max())
    ax.grid(True)
    ax.set_ylabel('AMOC (Sv)')
    ax.set_xlabel('integration year')
    ax.xaxis_date()

    if savefig or figname:
        figname = figname or 'AMOC_time_series.png'
        fig.savefig(figname,dpi=300)
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
    parser.add_argument('-l','--latlim', type=parse_coords,
            help='Latitude limits for AMOC region, e.g. 30,60. Accepts negative values as e.g. m30 for -30.', 
            default=(30,60))
    parser.add_argument('-z', '--zlim', type=lambda s: map(float,s.split(',')),
            help='Depth limits for AMOC region, e.g. 500,9999 for below 500 metres.',
            default=(500,9999))
    parser.add_argument('-s', dest='savefig', action='store_true',
            help='set this to save the figure with default filename')
    parser.add_argument('--figname', help='save the figure to the given filename (implies -s)')
    parser.add_argument('--savetofile', dest='outfname', 
            help='save the data to this file (pickle)')
    
    args = parser.parse_args()
    argsdict = vars(args)

    if len(args.files) == 1:
        argsdict['files'] = sorted(glob.glob(args.files[0]))
 
    if args.outfname is not None:
        for key in ['savefig','figname']:
            try:
                del argsdict[key]
            except KeyError:
                pass
        save_amoc_time_series(**argsdict)
    else:
        del argsdict['outfname']
        plot_amoc_time_series(**argsdict)
 
