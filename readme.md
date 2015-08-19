# Crane

Pre-/post processing tools for ocean modeling products.

## Roadmap

- fix all imports
    - should use absolute: import crane.data.timeArray
    - or explicit relatime: import .timeArray
    - avoid from x import * to keep namespaces intact
- clean unused functions
- pep8 all source
- data module is unfortunately named... change
- fix setup.py so that scripts get copied to correct bin dir
- revisit docstrings (use numpy convetion?)
- treat tappy better
    - fork tappy, morph to a generic installable package
    - or find alternative
- revisit plotting routines:
    - inherit from plot base classes
    - remove obsolete stackPlot routines
- re-implement extraction routines
    - process all variables simultaneously
    - save each day to disk immediately;
        i.e. do not store whole time period in memory
- remove all old station files, use only csv implementation
- use a database in stationCollection (in-memory): PyDbLite?
- move definitions to separate file
    - plot variable names/units
    - tracer names/file name conventions
- once relatively stable, progressively change naming conventions
- revisit dataContainer
    - transect data shape should be intuitive [nx, nz, nComp, nTime]
    - requires generization of data array?
    - use existing data model implementation? hdf5 object?
    - could meshContainer and dataContainer be merged?
        - dataContainer is a meshContainer with degenerate mesh
            (e.g. a point)

## Wishlist

- make quick plot routine for any data: plotDataContainer
- use matplotlib styles
- support interactive plots in ipython (via bokeh?)
    - time series and maps etc
- add support for different time zones
    - epoch time stamps are OK
    - all datetime objs are currently assumed to be in PST(!)
    - same for simulation time/matlab datenums/plot time
