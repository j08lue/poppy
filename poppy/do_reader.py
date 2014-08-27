import pandas as pd
from cStringIO import StringIO
import itertools
import numpy as np
import datetime
import glob



def _add_index(df,year=1):
    df['decyear'] = year + np.arange(len(df),dtype='f8')/len(df)
    df.set_index('decyear',inplace=True)
    df.index.name = 'Model year'
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

    ovf_TS = [StringIO() for s in straits]
    ovf_tr = [StringIO() for s in straits]
    lines_TS = []
    lines_tr = []

    # split lines into TS and tr
    with open(fname,'r') as fin:
        for line in fin:
            if line.startswith(' ovf_TS'):
                lines_TS.append(line[18:])
            elif line.startswith(' ovf_tr'):
                lines_tr.append(line[18:])

    # separate out straits
    for i in xrange(nstraits):
        # take every (nstraits)th line starting from 0,1,2,3...
        ovf_TS[i].writelines(itertools.islice(lines_TS,i,None,nstraits))
        ovf_tr[i].writelines(itertools.islice(lines_TS,i,None,nstraits))
        # rewind StringIOs
        ovf_TS[i].seek(0)
        ovf_tr[i].seek(0)

    # get start day
    year = year or _yr_from_fname(fname)

    # Read separate straits into panda dfs
    dfs_TS = []
    dfs_tr = []
    for i,strait in enumerate(straits):
        # TS
        df = pd.read_csv(ovf_TS[i],
                names=['n','Ti','Si','Ts','Ss','Te','Se','Tp','Sp'],
                delim_whitespace=True)
        _add_index(df,year)
        dfs_TS.append(df)
        # tr
        df = pd.read_csv(ovf_tr[i],
                names=['n','phi','Ms','Me','Mp','m','zt'],
                delim_whitespace=True)
        _add_index(df,year)
        dfs_tr.append(df)

    return dict(
            TS = pd.concat(dfs_TS,keys=straits,names=('Strait',dfs_TS[0].index.names[0])),
            tr = pd.concat(dfs_tr,keys=straits,names=('Strait',dfs_tr[0].index.names[0])),
        )


def read_do_multifile(files,datatype='TS',**kwargs):
    """Read multiple POP diagnostic overflow output (do) files 
    and concatenate their data.
    
    Parameters
    ----------
    files : list of str
        files to read and concatenate
    datatype : str in ['TS','tr']
        data type to return
    """
    if isinstance(files, basestring):
        files = [files]

    if len(files) == 1:
        files = sorted(glob.glob(files[0]))

    dfs = []
    for fname in files:
        data = read_do_file(fname,**kwargs)[datatype]
        dfs.append(data)
    df = pd.concat(dfs)
    df.files = files
    return df
