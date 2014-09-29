"""
This module contains a routine for parsing POP overflow (ovf)
parameterization namelist files.

These files exist for each ocean grid and contain information
about the location of grid points used by the ovf parameterization.
"""
from collections import defaultdict
import numpy as np
import itertools

verbose = False

class Struct:
    def __init__(self, entries={}, **kwargs):
        self.__dict__.update(entries)
        self.__dict__.update(kwargs)
    def __setitem__(self, k, v):
        self.__dict__[k] = v
    def __getitem__(self, k):
        return self.__dict__[k]
    def __repr__(self):
        return '\n'.join(
            ['{} : {}'.format(
                k,repr(v)) for k,v in sorted(self.__dict__.items())])
    def keys(self):
        return list(self.__dict__.keys())
    def iteritems(self):
        return sorted(self.__dict__.items())
    def items(self):
        return sorted(self.__dict__.items())


def parse_ovf_file(fname, outof_into_offset=False, output_zerobased=True):
    """Parse a POP overflow parameterization namelist file
    and return the grid indices of the cells involved
    
    Parameters
    ----------
    fname : str
        path to input file
    outof_into_offset : bool
        whether to add the offset to the indices that is stated in the
        direction bit:
                                             2
                                          ________
                      y ^                |        |
                        |               3|   ijk  |1
                        +--->            |        |
                            x            |________|
                                             4
    output_zerobased : bool
        whether to return the indices zero-based (Python)
        if False, they are returned 1-based (original, Fortran)
    """
    overflows = {}
    regionsdict = defaultdict(Struct)
    with open(fname,'r') as f:
        at_beginning = True
        regions_def = False
        nsets = 0
        for line in f:
            # beginning
            if at_beginning:
                try:
                    novf = int(line[1])
                    if verbose: print 'Number of overflows: {}'.format(novf)
                    at_beginning = False
                    continue
                except ValueError:
                    pass
            elif not line[1].isspace():
                # meta info start
                iovf = int(line[1])
                ovfname = ''.join(line[2:32].split()).replace('\'','')
                if verbose: print 'Parsing overflow {}'.format(ovfname)
                overflows[ovfname] = {}
                continue
            elif line[:5].isspace() and '! regional' in line:
                regions_def = True
                regkeys = itertools.cycle(['inflow', 'src', 'ent'])
                continue
            elif regions_def:
                key = next(regkeys)
                idx = map(int, line[5:30].split())
                for i, ind in zip(idx, ['imin', 'imax', 'jmin', 'jmax', 'kmin', 'kmax']):
                    if output_zerobased:
                        i -= 1
                    regionsdict[key][ind] = i
                if key == 'ent':
                    # end of regions definitions
                    regions_def = False
                    overflows[ovfname]['regions'] = regionsdict
                    regionsdict = defaultdict(Struct)
                    continue
            elif not line[3].isspace():
                if 'number of kmt changes' in line:
                    continue
                # start of set
                if '# prd sets' in line:
                    settype = 'prd'
                    nsets = int(line[3])
                    setdict = defaultdict(list)
                    continue
                elif 'prd' in line:
                    npts = int(line[3])
                    continue
                else:
                    nsets = 1
                    npts = int(line[3])
                    settype = line[34:37]
                    setdict = defaultdict(list)
                    continue
            elif nsets > 0:
                # in set
                i,j,k,d = map(int, line[4:17].split())
                if outof_into_offset:
                    if d == 1:
                        i+=1
                    elif d == 3:
                        i-=1
                    elif d == 2:
                        j+=1
                    elif d == 4:
                        j-=1
                setdict['ii'].append(i)
                setdict['jj'].append(j)
                setdict['kk'].append(k)
                npts -= 1
                if npts > 0:
                    continue
                else:
                    nsets -= 1
                    if nsets == 0:
                        for key in setdict.keys():
                            ind = np.array(setdict[key])
                            if output_zerobased: 
                                ind -= 1
                            setdict[key] = ind
                        overflows[ovfname][settype] = setdict
                        setdict = defaultdict(list)
                    continue
            elif nsets == 0:
                continue                
            else:
                continue

        if iovf != novf:
            print ('Warning: Something might have gone wrong. Initial '
                   'number of overflows not equal to the number of parsed overflows.')

    return overflows


