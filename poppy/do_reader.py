import pandas
from cStringIO import StringIO
import itertools
import numpy as np


def _add_tday_index(df,start=1):
    df['tday'] = np.arange(start,len(df)+start)
    df.set_index('tday')
    df.index.name = 'Day of year'
    return df

def _yr_from_fname(fname):
    # e.g. SOMECASENAME.pop.do.0529-01-01-00000
    try:
        return int(fname[-16:-12])
    except:
        print('Warning: Unable to get year from file name: '+fname)
        return 0


def read_do_file(fname,straits=['DS','FBC','RossSea','WeddellSea']):
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

    start = (_yr_from_fname(fname)-1)*365

    # Read separate straits into panda dfs
    data_TS = {}
    data_tr = {}
    for i,strait in enumerate(straits):
        # TS
        df = pandas.read_csv(ovf_TS[i],
                names=['n','Ti','Si','Ts','Ss','Te','Se','Tp','Sp'],
                delim_whitespace=True)
        df = _add_tday_index(df,start)
        data_TS[strait] = df
        # tr
        df = pandas.read_csv(ovf_tr[i],
                names=['n','phi','Ms','Me','Mp','m','zt'],
                delim_whitespace=True)
        df = _add_tday_index(df,start)
        data_tr[strait] = df

    return dict(TS=data_TS,tr=data_tr)

