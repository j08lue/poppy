#!/usr/bin/env python
"""
Plot AMOC time series from h5 files

"""
from __future__ import print_function
import matplotlib.pyplot as plt
import pandas as pd
import os.path

# some style settings
plt.style.use('ggplot') ; plt.rc('lines', linewidth=2)
try:
    import palettable as pal
    colors = pal.tableau.Tableau_10.mpl_colors
    plt.rc('axes', color_cycle=colors)
except ImportError:
    pass


def main(files, maxyear=None, savefig=None):
    """Plot AMOC time series from h5 files
    """

    fig, ax = plt.subplots(figsize=(10,6))
    ax.set_ylabel('AMOC maximum (Sv)')

    kwargs = dict(key='df')
    if maxyear is not None:
        kwargs['where'] = 'index<={}'.format(maxyear)

    labels = []
    for fname in sorted(files):
        case = os.path.basename(fname).split('.AMOC.h5')[0]
        labels.append(case)
        print(case)
        df = pd.read_hdf(fname, **kwargs) 
        df.plot(label=case)#, ls=('--' if '_br' in fname else '-'))
    ax.legend(labels=labels,
              #loc='center left', bbox_to_anchor=(1, 0.5),
              loc=0,
             )
    ax.axhline(17, ls='--', color='grey')

    plt.show()

    if savefig is not None:
        if savefig.lower().endswith('.pdf'):
            fig.savefig(savefig, bbox_inches='tight')
        else:
            fig.savefig(savefig, dpi=200, bbox_inches='tight')

    return fig


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
            description="Plot AMOC time series from post-processed CESM/POP data")
    parser.add_argument('files', nargs='+', help='post-processed HDF5 files .h5')
    parser.add_argument('--maxyear', type=float, help='read only data until this model year')
    parser.add_argument('--savefig', type=str, help='save figure to this file')
    args = parser.parse_args()
 
    main(**vars(args))
