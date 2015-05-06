#!/usr/bin/env python

import numpy as np
import argparse
import glob
import netCDF4

import poppy.metrics
import poppy.grid

import meta


def _get_grid_cell_area(fname, grid, lonlim, latlim):
    with netCDF4.Dataset(fname) as ds:
        mask = poppy.grid.get_mask_lonlat(ds, lonlim=lonlim, latlim=latlim, grid=grid)
        mask &= ds.variables['KM'+grid][:]>0
        jj,ii = np.where(mask)
        return ds.variables[grid+'AREA'][jj,ii]


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description="Save time series of a variable, averaged over a region")
    parser.add_argument('files', type=str, nargs='+',
            help='Files to read and concatenate')
    parser.add_argument('--varn', type=str,
            help='Variable name')
    parser.add_argument('--grid', type=str,
            help='Grid', choices=['T','U'])
    parser.add_argument('--region', type=str,
            help='Region name', choices=meta.regionlims.keys())
    parser.add_argument('--metric', type=str,
            help='Metric', choices=['mean', 'integral'], default='mean')
    parser.add_argument('-o', '--outfile', type=str,
            help='Output file')
    args = parser.parse_args()

    if isinstance(args.files, basestring):
        files = [args.files]
    
    files = sorted(args.files)

    if len(files) == 1:
        files = sorted(glob.glob(files[0]))

    area = _get_grid_cell_area(files[0], grid=args.grid, **meta.regionlims[args.region])
    def _weightedmean(data):
        return np.sum(data*area, axis=-1) / np.sum(area)
    def _integral(data):
        return np.sum(data*area*1e-4, axis=-1) 

    if args.metric == 'mean':
        reducefunc = _weightedmean
    elif args.metric == 'integral':
        reducefunc = _integral
    
    ts = poppy.metrics.get_timeseries(
            files,
            varn=args.varn,
            grid=args.grid,
            reducefunc=reducefunc,
            **meta.regionlims[args.region])
    
    ts.to_hdf(args.outfile, key='{0.varn}_{0.region}'.format(args), mode='w', format='table')

    print 'Data saved to {}'.format(args.outfile)
