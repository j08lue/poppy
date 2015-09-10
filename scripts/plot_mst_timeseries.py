#!/usr/bin/env python

import matplotlib.pyplot as plt ; plt.close('all')
import cPickle as pickle
import glob

import poppy.metrics
poppy.metrics.use_pandas = False

def save_mst_time_series(pattern, outfname, lat0=55, component=0):
    """Save MST time series from CESM/POP data"""
    ncfiles = sorted(glob.glob(pattern))
    meannsalt, timeax = poppy.metrics.get_mst(ncfiles, lat0=55, component=0)
    dataout = dict(
            component=component,lat0=lat0,timeax=timeax,
            meannsalt=meannsalt,ncfiles=ncfiles)
    with open(outfname,'wb') as fout:
        pickle.dump(dataout,fout)


def plot_mst_time_series(pattern, lat0=55, component=0, savefig=False, figname=''):
    """Plot MST time series from CESM/POP data"""
    ncfiles = sorted(glob.glob(pattern))
    meannsalt, timeax = poppy.metrics.get_mst(ncfiles, lat0=55, component=0)
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(timeax,meannsalt,'r')
    ax.grid(True)
    ax.set_ylabel('Meridional salt transport (Sv PPT)')
    ax.set_xlabel('integration year')
    ax.set_title('{} meridional salt transport at {}N'.format(
        poppy.metrics.componentnames[component],lat0))
    ax.xaxis_date()
    #fig.autofmt_xdate()

    if savefig or figname:
        figname = figname or 'MST_time_series.png'
        fig.savefig(figname,dpi=300)
    else:
        plt.show()


if __name__ == "__main__":

    defaults = dict(
        pattern = '*.pop.h.????-??.nc',
        component = 0,
        lat0 = 55,
        )

    def parse_coords(s):
        result = map(float,(s.replace('m','-')).split(','))
        return result[0] if len(result) == 0 else result

    import argparse
    parser = argparse.ArgumentParser(description="Plot N_SALT time series from CESM/POP data")
    parser.add_argument('pattern', nargs='+',
            help='file search pattern. Provide inside "" from command line!')
    parser.add_argument('-c','--component',type=int,choices=[0,1,2,3,4],help='transport component')
    parser.add_argument('-l','--lat0',type=parse_coords,help='latitude where to extract salt transport')
    parser.add_argument('-s','--savefig',action='store_true',help='set this to save the figure with default filename')
    parser.add_argument('-f','--figname',type=str,help='save the figure to the given filename (implies -s)')
    parser.add_argument('--savetofile',dest='outfname',type=str,help='save the data to this file (pickle)')
    
    parser.set_defaults(**defaults)
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
        save_mst_time_series(**argsdict)
    else:
        del argsdict['outfname']
        plot_mst_time_series(**argsdict)
