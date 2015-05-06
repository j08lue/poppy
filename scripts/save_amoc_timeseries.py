#!/usr/bin/env python

import argparse
import os.path
import glob
import cPickle as pickle

import poppy.metrics

if __name__ == "__main__":

    def parse_coords(s):
        result = map(float,(s.replace('m','-')).split(','))
        return result[0] if len(result) == 0 else result

    parser = argparse.ArgumentParser(
            description="Extract and save AMOC time series from CESM/POP data")
    parser.add_argument('files', type=str, nargs='+',
            help='Files to read and concatenate')
    parser.add_argument('-l', '--latlim', type=parse_coords,
            help=(
                'Latitude limits for AMOC region, e.g. 30,60.'
                ' Takes negative values as e.g. m30 for -30.'),
            default=(30,60))
    parser.add_argument('-z', '--zlim', type=lambda s: map(float, s.split(',')),
            help='Depth limits for AMOC region, e.g. 500,9999 for below 500 metres.',
            default=(500,9999))
    parser.add_argument('--nosort', action='store_true', 
            help='Disable alphabetic sorting')
    parser.add_argument('-o', '--outfile', type=str, 
            help='Output file')
    args = parser.parse_args()

    if not args.nosort:
        files = sorted(args.files)

    if len(args.files) == 1:
        args.files = sorted(glob.glob(args.files[0]))

    df = poppy.metrics.get_amoc(args.files, latlim=args.latlim, zlim=args.zlim)

    if os.path.splitext(args.outfile)[-1] == '.h5':
        if not poppy.metrics.use_pandas:
            raise NotImplementedError('Saving to HDF5 requires Pandas!')
        df.to_hdf(args.outfile,key='df',mode='w',format='table')
    else:
        with open(args.outfile,'wb') as fout:
            pickle.dump(df,fout)
