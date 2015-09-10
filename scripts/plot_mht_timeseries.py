#!/usr/bin/env python

import matplotlib.pyplot as plt ; plt.close('all')
import cPickle as pickle
import glob

import poppy.metrics
poppy.metrics.use_pandas = False

def save_mht_time_series(pattern,outfname,latlim=(30,60),component=0):
    """Store MHT time series from CESM/POP data"""
    ncfiles = sorted(glob.glob(pattern))
    maxmeannheat, timeax = poppy.metrics.get_mht(ncfiles,latlim=latlim,component=component)
    dataout = dict(
            component=component,latlim=latlim,timeax=timeax,
            maxmeannheat=maxmeannheat,ncfiles=ncfiles)
    with open(outfname,'wb') as fout:
        pickle.dump(dataout,fout)


def plot_mht_time_series(pattern,latlim=(30,60),component=0,savefig=False,figname=''):
    """Plot MHT time series from CESM/POP data"""
    ncfiles = sorted(glob.glob(pattern))
    maxmeannheat, timeax = poppy.metrics.get_mht(ncfiles,latlim=latlim,component=component)
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(timeax,maxmeannheat,'g')
    ax.grid(True)
    ax.set_ylabel('Meridional heat transport (PW)')
    ax.set_xlabel('integration year')
    ax.set_title('{0} MHT maximum between {1[0]}N and {1[1]}N'.format(
        poppy.metrics.componentnames[component],latlim))
    ax.xaxis_date()
    #fig.autofmt_xdate()
    
    if savefig or figname:
        figname = figname or 'MHT_time_series.png'
        fig.savefig(figname,dpi=300)
    else:
        plt.show()


if __name__ == "__main__":

    defaults = dict(
        pattern = '*.pop.h.????-??.nc',
        component = 0,
        latlim=(30,60),
        )

    def parse_coords(s):
        result = map(float,(s.replace('m','-')).split(','))
        return result[0] if len(result) == 0 else result

    import argparse
    parser = argparse.ArgumentParser(description="Plot N_HEAT time series from CESM/POP data")
    parser.add_argument('pattern', nargs='+',
            help='file search pattern. Provide inside "" from command line!')
    parser.add_argument('-c','--component',type=int,choices=[0,1,2,3,4],help='transport component')
    parser.add_argument('-s','--savefig',action='store_true',help='set this to save the figure with default filename',default=False)
    parser.add_argument('-l','--latlim',type=parse_coords,help='latitude limits where to search for max mht')
    parser.add_argument('--figname',type=str,help='save the figure to the given filename (implies -s)')
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
        save_mht_time_series(**argsdict)
    else:
        del argsdict['outfname']
        plot_mht_time_series(**argsdict)
