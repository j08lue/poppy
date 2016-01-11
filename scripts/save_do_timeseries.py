#!/usr/bin/env python

from __future__ import print_function
import argparse
import glob
try:
    import cPickle as pickle
except ImportError:
    import pickle

from poppy import do_reader

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Read and concatenate *.pop.do.* files and save data series to file')
    parser.add_argument('files', type=str, nargs='+', help='Files to read and concatenate')
    parser.add_argument('-o', '--outfile', type=str, help='Output file')
    args = parser.parse_args()

    if isinstance(args.files, basestring):
        files = sorted(glob.glob(args.files))
    
    files = sorted(args.files)
    print('Processing {} files ...'.format(len(files)))

    df = do_reader.read_do_multifile(files=files)

    if args.outfile.endswith('.h5'):
        df.to_hdf(args.outfile, key='df', mode='w', format='table')
    else:
        print('    Using pickle ...')
        with open(args.outfile,'wb') as fout:
            pickle.dump(df,fout)
