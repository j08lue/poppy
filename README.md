poppy
=====

Python tools for processing output data from the Parallel Ocean Program (POP)

Strong as [Popeye the Sailor Man](http://en.wikipedia.org/wiki/Popeye)

Installation
------------

    git clone https://github.com/j08lue/poppy.git
    cd poppy
    python setup.py install

Modules
-------
The `poppy` package contains modules for various data postprocessing tasks. Of particular interest might be the `metrics` package, that contains routines for extracting large-scale ocean metrics such as the maximum AMOC strength or meridional heat transport from an arbitrary number of files. The metrics are stored as Pandas dataframes in HDF5 (if Pandas is available) or Pickled NumPy arrays and can be easily plotted together.

Scripts
-------
The `scripts` directory contains mainly command-line interfaces for the different `metrics` functions, e.g. to plot the AMOC strength evolution directly from the model output files:

    cd <CASENAME>/ocn/hist
    plot_amoc_timeseries.py *.pop.h.????-??.nc

nclookipy
---------
The `nclook.py` script offers a quicklook into netCDF files. If you create an alias like

    alias nclookipy='ipython -i /path/to/nclook.py'
    
you can call e.g.

    nclookipy *.pop.h.0100-??.nc

and have the files opened and ready to explore as a multi-file `xray` or `netCDF4` Dataset.


Requirements
------------
Please see the `setup.py` for dependencies.
Written in python2 but works well with `2to3`. One method to make this work is to rename the `poppy` folder to e.g. `poppy_py2`, and then calling

    2to3 poppy_py2 -o poppy -W -n
    
It should then be possible to install the package with `python setup.py install`.
