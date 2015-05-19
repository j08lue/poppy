import xray
import glob


def read_netcdfs(files, dim, transform_func=None, **kwarg):
    def process_one_path(path):
        # use a context manager, to ensure the file gets closed after use
        with xray.open_dataset(path, **kwarg) as ds:
            # transform_func should do some sort of selection or
            # aggregation
            if transform_func is not None:
                ds = transform_func(ds)
            # load all data from the transformed dataset, to ensure we can
            # use it after closing each original file
            ds.load_data()
            return ds
    
    if isinstance(files, list):
        paths = []
        for p in files:
            paths += sorted(glob.glob(p))
    else:
        paths = sorted(glob.glob(files))
    print paths
    datasets = [process_one_path(p) for p in paths]
    xray.concat(datasets, dim)

if __name__ == "__main__":

    import argparse
    import traceback

    parser = argparse.ArgumentParser(description="Open one or several netCDF files from command line")
    parser.add_argument('ncfiles', nargs='+', help='netCDF files')
    args = parser.parse_args()

    try:
        ds = read_netcdfs(args.ncfiles, dim='time', decode_times=False)
        print('opened {} netCDF file into namespace as \'ds\''.format(len(args.ncfiles)))
    except RuntimeError:
        traceback.print_exc()
        print('Warning: File(s) could not be opened: {}'.format(args.ncfiles))
