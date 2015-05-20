#!/usr/env python
"""
This script can be executed as is or via an alias, e.g.

    alias nclookipy='ipython -i ~/path/to/nclook.py'

to quickly open a NetCDF file in Python.
"""
import netCDF4
import argparse
import traceback
import os.path
import matplotlib.pyplot as plt
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description="Open one or several netCDF files from command line")
parser.add_argument('ncfiles', nargs='+', help='path to netCDF file')
args = parser.parse_args()

try:
    ds = netCDF4.MFDataset(args.ncfiles)
    dsvar = ds.variables
    print('opened {} netCDF file into namespace as \'ds\''.format(len(args.ncfiles)))
except RuntimeError:
    traceback.print_exc()
    print('Warning: File(s) could not be opened: {}'.format(args.ncfiles))

def quickplot(varn, k, ds=ds, t=0, cmap='cubehelix', **kwarg):
    dsvar = ds.variables
    fig = plt.figure()
    plt.title('{} at {:.0f} m'.format(varn, dsvar['z_t'][k]*1e-2))
    plt.pcolormesh(dsvar[varn][t,k], cmap=cmap, **kwarg)
    plt.colorbar(label='{} ({})'.format(dsvar[varn].long_name, dsvar[varn].units))
    plt.show()
    return fig

def save_figure(fname, fig=None):
    if fig is None:
        fig = plt.gca().figure
    fname = os.path.expanduser(fname)
    if fname.endswith('.pdf'):
        fig.savefig(fname, bbox_inches='tight')
    else:
        fig.savefig(fname, dpi=200)

def print_varmeta(varn):
    """Print variable meta information"""
    print(ds.variables[varn])

def list_variables(units=True,ndim=None,shapes=True):
    """List variables in netCDF *dataset* including units if *units* for variables in *ndim* dimensions only, if specified"""
    def _get_variable_list(ds):
        for varn in ds.variables.keys():
            if ndim is not None:
                if len(ds.variables[varn].shape) != ndim:
                    continue
            s = ''
            try:
                s += '{} :'.format(varn)
            except:
                pass
            try:
                s += ' {}'.format(ds.variables[varn].long_name)
            except:
                pass
            try:
                varunits = ds.variables[varn].units
                if units and varunits:
                    s += ' [{}]'.format(varunits)
            except:
                pass
            if shapes:
                try:
                    s += ' {}'.format(ds.variables[varn].shape)
                except:
                    pass
                # print output
            if s: 
                print(s)
    try:
        with netCDF4.Dataset(ds) as ds_new:
            _get_variable_list(ds_new)
    except:
        _get_variable_list(ds)

