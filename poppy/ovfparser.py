"""
This module contains a routine for parsing POP overflow (ovf)
parameterization namelist files.

These files exist for each ocean grid and contain information
about the location of grid points used by the ovf parameterization.
"""
from collections import defaultdict
import numpy as np

verbose = False

def parse_ovf_file(fname, zerobased=True):
    """Parse a POP overflow parameterization namelist file
    and return the grid indices of the cells involved"""
    overflows = {}
    with open(fname,'r') as f:
        at_beginning = True
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
            elif line[1] != ' ':
                # meta info start
                iovf = int(line[1])
                ovfname = ''.join(line[2:32].split()).replace('\'','')
                if verbose: print 'Parsing overflow {}'.format(ovfname)
                overflows[ovfname] = {}
                continue
            elif line[3] != ' ':
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
                indx = line[4:17].split()
                i,j,k,d = [int(v) for v in indx]
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
                            if zerobased: 
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


