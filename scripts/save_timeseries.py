#!/usr/bin/env python

from __future__ import print_function
import numpy as np
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


def get_timeseries(files, varn, grid, region, func, outfile):
    if len(files) == 1:
        files = glob.glob(files[0])
    ts = poppy.metrics.get_timeseries(
            sorted(files),
            varn=varn,
            grid=grid,
            reducefunc=getattr(np, func),
            **regionlims[region])
    ts.to_hdf(outfile, key=varn+'_'+region, mode='w', format='table')
    print('Data saved to {}'.format(outfile))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Save time series of a given POP variable")
    parser.add_argument('files', nargs='+', help='Files to read and concatenate')
    parser.add_argument('-v', '--varn', help='Variable name')
    parser.add_argument('-g', '--grid', help='Grid', choices=['T', 'U'])
    parser.add_argument('-r', '--region', help='Region name', choices=regionlims.keys())
    parser.add_argument('-f', '--func', help='Name of NumPy function to reduce data along time axis (e.g. min,max,mean)')
    parser.add_argument('-o', '--outfile', help='Output file', default='max.h5')
    args = parser.parse_args()
    get_timeseries(**vars(args))

