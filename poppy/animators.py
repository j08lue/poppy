import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as mpatheffects
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
import time
import datetime

import pyutils.dates as pud
import pycpt.modify
import pycpt.load
pycpt.load.register_cptcity_cmaps('http://soliton.vm.bytemark.co.uk/pub/cpt-city/gmt/GMT_ocean.cpt')

def _update_time(ani,fname):
    with netCDF4.Dataset(fname) as ds:
        timevar = ds.variables['time']
        ani.date = pud.pseudo_to_proper_datetime(
                netCDF4.num2date(timevar[0],timevar.units,timevar.calendar))
        ani.date -= datetime.timedelta(days=1)
        newyear = pud.datetime_to_decimal_year(ani.date)
        try:
            ani.elapsed = newyear - ani.initialyear
        except AttributeError:
            ani.elapsed = 0.
            ani.initialyear = newyear
        ani.year = newyear

def _update_timestamp(ani, fname):
    _update_time(ani, fname)
    text = 'model year: {0.year:04d}-{0.month:02d}\nelapsed: {1:.2f} years'.format(ani.date,ani.elapsed)
    try:
        ani.timestamp.set_text(text)
    except AttributeError:
        ani.timestamp = ani.ax.text(0.95,0.05,text,color='k',
                ha='right',va='bottom',transform=ani.ax.transAxes,
                path_effects=([mpatheffects.withStroke(linewidth=3, foreground='white')]))

def _find_depth_level(fname, depth):
    with netCDF4.Dataset(fname) as ds:
        return np.argmin(np.abs(ds.variables['z_w'][:]*1e-2 - depth))

def _get_level_depth(fname, k):
    with netCDF4.Dataset(fname) as ds:
        return ds.variables['z_w'][k]*1e-2

class Layer:
    """Horizontal layer of a given 4D variable"""
    def __init__(self,
            ax,
            ncfiles,
            varname,
            scale=1.,
            levels=None,
            cmap = None,
            t=0,
            ii=None, jj=None,
            k=0,
            depth=None, 
            pause = 0,
            with_timestamp = True,
            ):
        self.ax = ax
        self.ncfiles = ncfiles
        self.varname = varname
        self.scale = scale
        self.t = t
        #depth
        self.depth = depth
        if depth is not None:
            self.k = _find_depth_level(ncfiles[0], depth)
        elif k is not None:
            self.k = k
        else:
            raise ValueError('Either k or depth must be specified.')
        self.depth_k = _get_level_depth(ncfiles[0], self.k)
        # data shape
        self.datashape = self._get_datashape()
        self.ndim = len(self.datashape)
        if self.ndim not in [3,4]:
            raise NotImplementedError('Only 4D and 3D variables supported!')
        self.ny,self.nx = ny,nx = self.datashape[-2:]
        self.ii = np.arange(nx) if ii is None else np.mod(ii,nx)
        self.ii_original = ii
        self.jj = np.arange(ny) if jj is None else jj
        self.fig = self.ax.get_figure()
        self.pause = pause
        self.with_timestamp = with_timestamp
        self._make_axes()
        self._update_long_name(ncfiles[0])
        self.ax.autoscale(axis='both',tight=True)
        # levels
        if levels is None:
            sample_data = self._get_data(self.ncfiles[-1])
            ticker = mticker.MaxNLocator(nbins=21, symmetric=True)
            levels = ticker.tick_values(sample_data.min(), sample_data.max()) 
        self.cmap, self.norm = pycpt.modify.generate_cmap_norm(levels=levels, cm=(cmap or 'RdYlBu_r'))
    
    def _update_long_name(self,fname):
        with netCDF4.Dataset(fname) as ds:
            self.long_name = '{0.long_name} ({0.units})'.format(ds.variables[self.varname])
            if len(self.datashape) == 4:
                self.long_name += ' at {:.0f} m'.format(self.depth_k)

    def init(self,cbarpos='right'):
        fname = self.ncfiles[0]
        self.data = self._get_data(fname)
        self.img = self.ax.pcolormesh(
                self.xax,self.yax,
                self.data,
                cmap=self.cmap,norm=self.norm)
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes(cbarpos, size="5%", pad=0.05)
        self.cb = self.fig.colorbar(self.img, cax=cax, orientation='vertical',
                #format='%.1e',
                #label = self.long_name,
                )
        self.cb.ax.yaxis.set_ticks_position(cbarpos)
        #self.cb.ax.yaxis.set_label_position(cbarpos)
        if self.with_timestamp: _update_timestamp(self,fname)
        self.ax.set_xticklabels(['{:.0f}'.format(np.mod(i,self.nx)) for i in self.ax.get_xticks()])
        self.ax.set_title(self.long_name)
        return self.img

    def _get_datashape(self, fname=None):
        fname = fname or self.ncfiles[0]
        with netCDF4.Dataset(fname) as ds:
            self.datashape = ds.variables[self.varname].shape
            return self.datashape

    def _get_data(self, fname):
        with netCDF4.Dataset(fname) as ds:
            if self.ndim == 4:
                return ds.variables[self.varname][self.t,self.k,self.jj,self.ii] * self.scale
            elif self.ndim == 3:
                return ds.variables[self.varname][self.t,self.jj,self.ii] * self.scale

    def _make_axes(self):
        try:
            self.xax = np.concatenate((self.ii_original,[self.ii_original[-1]+1]))-0.5
            self.yax = self.jj
        except TypeError:
            self.xax = np.arange(len(self.ii)+1,dtype=float)-0.5
            self.yax = np.arange(len(self.jj)+1,dtype=float)-0.5

    def __call__(self,i):
        fname = self.ncfiles[i]
        self.data = self._get_data(fname)
        time.sleep(self.pause)
        self.img.set_array(self.data.ravel())
        if self.with_timestamp: _update_timestamp(self,fname)
        return self.img


class VerticalProfile:
    """Vertical profile of a given variable"""
    def __init__(self,
            ax,
            ncfiles,
            ii,jj,
            t = 0,
            varname = 'IAGE',
            levels = np.linspace(0,300,21),
            pause = 0,
            with_timestamp = True,
            ):
        self.ax = ax
        self.t = t
        self.ncfiles = ncfiles
        self.varname = varname
        self._update_long_name(ncfiles[0])
        self.jj = jj
        self.datashape = self._get_datashape(ncfiles[0])
        self.ii = np.mod(ii,self.datashape[3])
        self.fig = self.ax.get_figure()
        self.pause = pause
        self.with_timestamp = with_timestamp
        self._make_axes(ncfiles[0])
        self.ax.autoscale(axis='x',tight=True)
        self.ax.invert_yaxis()
        self.cmap,self.norm = pycpt.modify.generate_cmap_norm(levels,'RdYlBu_r')
    
    def _update_long_name(self,fname):
        with netCDF4.Dataset(fname) as ds:
            self.long_name = '{0.long_name} ({0.units})'.format(ds.variables[self.varname])

    def init(self,cbarpos='right'):
        fname = self.ncfiles[0]
        self.data = self._get_data(fname)
        self.img = self.ax.pcolormesh(
                self.xax,self.zax,
                self.data,
                cmap=self.cmap,norm=self.norm)
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes(cbarpos, size="5%", pad=0.05)
        self.cb = self.fig.colorbar(self.img,cax=cax,orientation='vertical')
        self.cb.set_label(self.long_name)
        if self.with_timestamp: _update_timestamp(self,fname)
        return self.img

    def _get_datashape(self,fname=None):
        fname = fname or self.ncfiles[0]
        with netCDF4.Dataset(fname) as ds:
            return ds.variables[self.varname].shape

    def _get_data(self,fname):
        with netCDF4.Dataset(fname) as ds:
            return ds.variables[self.varname][self.t,:,self.jj,self.ii]

    def _make_axes(self,fname):
        with netCDF4.Dataset(fname) as ds:
            dsvar = ds.variables
            self.zax = dsvar['z_w_bot'][:] * 1e-2
            if isinstance(self.ii,int):
                self.xax = np.arange(len(self.jj))
            elif isinstance(self.jj,int):
                self.xax = np.arange(len(self.ii))
            else:
                self.xax = np.arange(len(self.ii)+1,dtype=float)-0.5

    def __call__(self,i):
        fname = self.ncfiles[i]
        self.data = self._get_data(fname)
        time.sleep(self.pause)
        self.img.set_array(self.data.ravel())
        if self.with_timestamp: _update_timestamp(self,fname)
        return self.img

    def plot_map(self):
        self.mapfig = plt.figure()
        self.mapax = self.mapfig.add_subplot(111)
        with netCDF4.Dataset(self.ncfiles[0])as ds:
            depth = ds.variables['HU'][:]*1e-2
            depth = np.ma.masked_where(depth<=0,depth)
        self.mapimg = self.mapax.pcolormesh(depth,cmap='GMT_ocean_r')
        self.mapax.plot(self.ii,self.jj,'m',lw=2)
        return self.mapfig


