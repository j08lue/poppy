import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os.path
import pandas as pd
from pandas.stats.moments import rolling_mean
import gsw

try:
    import seaborn as sns
except ImportError:
    plt.style.use('ggplot')


def plot_TS(df):
    fig = plt.figure()
    ax = fig.gca()
    df.plot(x='Ti', y='Si', style='o', ax=ax)
    df.plot(x='Te', y='Se', style='o', ax=ax)
    df.plot(x='Tp', y='Sp', style='o', ax=ax)
    plt.show()



fieldsets = [
        ['Ts','Te','Tp'],
        ['Ss','Se','Sp'],
        ['rho_s','rho_e','rho_p'],
        ['Ms','Me','Mp'],
        ]

# observations as compiled in Danabasoglu et al., 2010 (DOI:10.1029/2010JC006243)
observations = {
        'DS' : {
            'Ts': 0.25, 'Tp': 2.1,
            'Ss': 34.81, 'Sp': 34.84,
            'Ms': (2.7, 3.8), 'Me': 2.3, 'Mp': 5.2,
            'rho_s': 27.94, 'rho_p': 27.85,
            'zt': (1600, 9999) },
        'FBC': {
            'Ts': 0, 'Tp': 3.3,
            'Ss' : 34.92, 'Sp' : (34.9, 35.15),
            'Ms': (1.5, 3.5), 'Me': 1.5, 'Mp': (2.7, 3.6),
            'rho_s': 28.07, 'rho_p': 27.9,
            'zt': (3000, 9999) }}

obs_props = dict(facecolor='0.6666', alpha=0.5)


def plot_do(files, strait='DS', maxyear=None, savefig=False, figname=None):
    files = sorted(files)

    # make sure all files exist
    for fname in files:
        if not os.path.isfile(fname):
            raise IOError('File not found: ' + fname)
    
    # keyword arguments for reading the files
    kwargs = dict(key='df')
    if maxyear:
        kwargs['where'] = 'ModelYear<={}'.format(maxyear)

    fig,axx = plt.subplots(
            nrows=len(fieldsets), ncols=len(fieldsets[0]), 
            sharex='all', sharey='row', figsize=(16,9))
    spt = fig.suptitle('Strait: {}'.format(strait))

    cases = [os.path.basename(fname).split('.do.h5')[0] for fname in files]

    # loop through files
    for nf, fname in enumerate(files):
        label = cases[nf]

        df = pd.read_hdf(fname, **kwargs)
        df = df.loc[strait]

        for i, fields in enumerate(fieldsets):
            for j, varn in enumerate(fields):
                ax = axx[i,j]
                ax.set_title(varn)
                if varn.startswith('rho_'):
                    # compute density from T,S
                    salt = rolling_mean(df[varn.replace('rho_','S')],365).values
                    temp = rolling_mean(df[varn.replace('rho_','T')],365).values
                    series = gsw.rho(salt, temp, 0)-1e3
                else:
                    # get series directly from file
                    series = rolling_mean(df[varn],365).values
                x = df.index.get_level_values('ModelYear')
                ax.plot(x, series, label=label)
            
                # only once
                if nf == 0:
                    # plot observations
                    try:
                        values = observations[strait][varn]
                        try:
                            obs_handle = ax.axhspan(values[0], values[1], **obs_props)
                        except TypeError:
                            ax.axhline(values, color='0.6666', linewidth=2)
                    except KeyError:
                        pass

    # legend with patch for observations
    handles, labels = axx.flat[0].get_legend_handles_labels()
    obs_handle = mpatches.Patch(**obs_props)
    handles += [obs_handle]
    labels += ['observations']
    fig.subplots_adjust(right=0.8)
    lgd = fig.legend(handles, labels, 
            bbox_to_anchor=(0.82,0.5), loc='center left',
            bbox_transform=fig.transFigure)

    # save figure to file
    if savefig or figname:
        figname = figname or 'ovf_props_{}_{}.pdf'.format(strait,'_'.join(cases))
        fig.savefig(figname, bbox_extra_artists=(lgd,spt,), bbox_inches='tight')
    else:
        plt.show()



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(
            description="Plot overflow diagnostics from post-processed CESM/POP data")
    parser.add_argument('files', nargs='+', help='post-processed HDF5 files .h5')
    parser.add_argument('--strait', type=str, default='DS', help='strait abbreviation (DS, FBC, ...)')
    parser.add_argument('--maxyear', type=float, help='read only data until this model year')
    parser.add_argument('--figname', type=str, help='figure name to use instead of default')
    parser.add_argument('-s', '--savefig', action='store_true', help='set this to save the figure')
    args = parser.parse_args()
 
    plot_do(**vars(args))
