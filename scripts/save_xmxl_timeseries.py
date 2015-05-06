#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse
import glob

import poppy.metrics

regionlims = {
        'Global' : dict(lonlim=(-90,90),latlim=(-180,180)),
        'Atlantic' : dict(lonlim=(-80,40),latlim=(-70,70)),
        'PolarNorthAtlantic': dict(lonlim=(-80,60),latlim=(60,90)),
        'LabradorSea': dict(latlim=(50,60), lonlim=(-50,-40)),
        'NorthAtlantic' : dict(lonlim=(-80,20),latlim=(0,65)),
        'SubpolarNorthAtlantic': dict(lonlim=(-60,0),latlim=(40,65)),
        'SubtropicalNorthAtlantic': dict(lonlim=(-100,0),latlim=(10,30)),
        'TreguierNorthAtlantic': dict(lonlim=(-100,20),latlim=(10,50)),
        'SubtropicalSouthAtlantic': dict(lonlim=(-50,20),latlim=(-40,-10)),
        'BrazilEastCoast20S40W': dict(lonlim=(-45,-20),latlim=(-30,-10)),
        'EquatorialAtlantic': dict(lonlim=(-55,15),latlim=(-10,10)),
        'SubtropicalSouthPacific': dict(lonlim=(150,280),latlim=(-40,-10)),
        }

def get_annual_max_xmxl(files, region):
    ts = poppy.metrics.get_timeseries(
            files,
            varn='XMXL',
            grid='T',
            reducefunc=np.max,
            **regionlims[region])
    ts /= 100.
    tssmooth = pd.rolling_max(ts, window=12, center=True)
    #poppy.metrics._pandas_copy_meta_data(ts, tssmooth, addmeta=dict(
    #    name = 'XMXL_'+region,
    #    units = 'm'
    #    ))
    return tssmooth


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description="Save annual maximum XMXL time series")
    parser.add_argument('files', type=str, nargs='+',
            help='Files to read and concatenate')
    parser.add_argument('-r', '--region', type=str,
            help='Region name', choices=regionlims.keys())
    parser.add_argument('-o', '--outfile', type=str,
            help='Output file',default='xmxl.h5')
    args = parser.parse_args()

    if isinstance(args.files, basestring):
        files = [args.files]
    
    files = sorted(args.files)

    if len(files) == 1:
        files = sorted(glob.glob(files[0]))

    ts = get_annual_max_xmxl(files, region=args.region)

    ts.to_hdf(args.outfile, key='XMXL_'+args.region, mode='w', format='table')

    print 'Data saved to {}'.format(args.outfile)
