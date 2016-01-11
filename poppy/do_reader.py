from __future__ import print_function
import pandas as pd
from io import StringIO
import itertools
import numpy as np
import datetime
import glob
from collections import defaultdict


def _add_index(df,year=1):
    df['decyear'] = year + np.arange(len(df),dtype='f8')/len(df)
    df.set_index('decyear',inplace=True)
    df.index.name = 'ModelYear'
    return df


def _yr_from_fname(fname):
    # e.g. SOMECASENAME.pop.do.0529-01-01-00000
    try:
        return int(fname[-16:-12])
    except:
        print('Warning: Unable to get year from file name: '+fname)
        return 0


def _date_from_fname(fname):
    # e.g. SOMECASENAME.pop.do.0529-01-01-00000
    try:
        return datetime.datetime.strftime(fname[-16:],'%Y-%m-%d-%H%M%S')
    except:
        print('Warning: Unable to get date from file name: '+fname)
        return datetime.datetime(year=0,month=1,day=1)


def read_do_file(fname,straits=['DS','FBC','RossSea','WeddellSea'],year=None):
    """Read a POP diagnostic overflow output (do) file
    
    These files are usually located in the `hist` directory when
    archived and in the `run` directory otherwise and match the 
    pattern `*pop.do.*`.

    Parameters
    ----------
    fname : str
        file name
    straits : list of str
        straits in the file (can be obtained from the POP overflow
        namelist for the respective grid)
        The length of this list must match the number of respective
        columns in the do file
    """
    nstraits = len(straits)
    contents = ['TS','tr']
    
    ios = {}
    lines = {}
    for key in contents:
        ios[key] = [StringIO() for s in straits]
        lines[key] = []

    # split lines into TS and tr
    with open(fname,'r') as fin:
        for line in fin:
            if line.startswith(' ovf_TS'):
                lines['TS'].append(line[18:])
            elif line.startswith(' ovf_tr'):
                lines['tr'].append(line[18:])

    # separate out straits
    for key in contents:
        for i in range(nstraits):
            # take every (nstraits)th line starting from 0,1,2,3...
            ios[key][i].writelines(itertools.islice(lines[key],i,None,nstraits))
            # rewind StringIOs
            ios[key][i].seek(0)

    # get start day
    year = year or _yr_from_fname(fname)

    cols = dict(
            TS = ['n','Ti','Si','Ts','Ss','Te','Se','Tp','Sp'],
            tr = ['n','phi','Ms','Me','Mp','m','zt'])
    usecols = dict(
            TS = np.arange(1,len(cols['TS'])),
            tr = np.arange(1,len(cols['tr'])))

    dfs = defaultdict(list)
    for key in contents:
        # select column names from usecols
        names = np.array(cols[key])[usecols[key]]
        # read separate straits into panda dfs
        for i,strait in enumerate(straits):
            df = pd.read_csv(ios[key][i],
                    names=names,
                    usecols=usecols[key],
                    delim_whitespace=True)
            _add_index(df,year)
            dfs[key].append(df)

    return pd.concat(
            [pd.concat(dfs[key],keys=straits,names=('Strait','ModelYear')) 
                for key in contents],axis=1)


def read_do_multifile(files,**kwargs):
    """Read multiple POP diagnostic overflow output (do) files 
    and concatenate their data.
    
    Parameters
    ----------
    files : list of str
        files to read and concatenate
    """
    if isinstance(files, str):
        files = [files]

    if len(files) == 1:
        files = sorted(glob.glob(files[0]))

    dfs = []
    for fname in files:
        df = read_do_file(fname,**kwargs)
        dfs.append(df)
    return pd.concat(dfs)
