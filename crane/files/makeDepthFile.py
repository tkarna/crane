#!/usr/bin/env python
"""Create a depth intregtated NetCDF *.61.nc files from *.63.nc files.

jesse.e.lopez
"""
#------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------
import os
import shutil
import argparse

import numpy as np
import netCDF4 as nc

import data.ncExtract as nce
import data.dirTreeManager as dtm
import data.meshContainer as mc

#------------------------------------------------------------------------------
# Constants
#------------------------------------------------------------------------------
DEPTH_INT_DESC = 'DI'
DEPTH_AVG_DESC = 'DAVG'
FILL_VALUE = -99.0


#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
def createDerivedNCFile(old, old_var, new_path, new_var, data, verbose=False):
    """Create a netcdf file using another as a template but with new data.

    Changes dimensions if going from (ntime, nvrt, nodes) to (ntime, nodes).

    Derived from: http://bit.ly/1vqox0k

    Params:
    -------
    old - ncfile
      NetCDF file to copy metadata from for new file
    old_var - string
      Variable name
    new_path - string
      Path to new NetCDF file to create
    new_var - string
      New variable name
    data - ndarray
      Array of data to save to new file in place of original
    verbose - bool (default: false)
      Print details about file creation
    """
    new = nc.Dataset(new_path, 'w', format='NETCDF4')
    # Copy dimensions
    for dim_name, dim in old.dimensions.iteritems():
        if verbose:
            print dim_name, len(dim)
        new.createDimension(dim_name, len(dim))

    # Copy variables except old_var
    for name, var in old.variables.iteritems():
        if name != old_var:
            out_var = new.createVariable(name, var.dtype, var.dimensions)
            out_var[:] = var[:]
            if verbose:
                print name, var.dimensions
        else:
            # 3D
            if len(var.dimensions) == len(data.shape):
                out_var = new.createVariable(new_var, var.dtype,
                                             var.dimensions)
            # 2D
            else:
                new_dims = (u'time', u'node')
                out_var = new.createVariable(new_var, var.dtype, new_dims)
            out_var[:] = data[:]
            if verbose:
                print new_var, new_dims
        # Copy variable attributes
        for name in vars(var):
            # Skip private
            if not name.startswith('_'):
                out_var.setncattr(name, var.getncattr(name))

    # Copy other attributes
    if hasattr(old, 'description'):
        new.description = old.description
    if hasattr(old, 'dataType'):
        new.dataType = old.dataType
    if hasattr(old, 'cruise'):
        new.cruise = old.cruise
    if hasattr(old, 'institute'):
        new.institue = old.institute
    if hasattr(old, 'location'):
        new.location = old.location
    if hasattr(old, 'ship'):
        new.ship = old.ship
    if hasattr(old, 'variable'):
        new.variable = old.variable
    if hasattr(old, 'tag'):
        new.tag = old.tag

    new.sync()
    new.close()


def makeDepthFile(data_file_path, method='avg'):
    """Create depth averaged or integrated file from 3D input fie.

    Params:
    -------
    data_file_path - string
      Path to 3D file to integrated
    Method - string
      Method to use over the vertical.  Either 'avg' or 'int'.
    """
    try:
        f = nce.selfeNCFile(data_file_path)
    except:
        print 'Problem opening file {0}'.format(data_file_path)
        raise
    f.readHeader()
    data_dir = os.path.dirname(data_file_path)
    filename = os.path.basename(data_file_path)
    split = filename.split('.')[0].split('_')
    stack = split[0]
    if len(split[1:]) > 1:
        var = '_'.join(split[1:])
    else:
        var = split[1]

    data = f.ncfile.variables[var][:]
    elev = f.ncfile.variables['elev'][:]
    bath = f.ncfile.variables['depth'][:]
    di_data = np.zeros_like(elev)
    ntimes = f.ncfile.variables['time'][:].shape[0]
    nnodes = f.ncfile.variables['elev'][:].shape[1]

    # Only integrate non-dry values (not the -99.0 fill value)
    good_3d_data = np.ma.masked_equal(data, FILL_VALUE)

    # Operate over values for every time step
    for t in np.arange(ntimes):
        # Get height of prisms
        Z, kbp2, iwet = f.vCoords.computeVerticalCoordinates(elev[t, :], bath)
        z_top = Z[1:, :]
        z_bot = Z[:-1, :]
        elem_height = z_top - z_bot
        total_height = np.ma.empty((nnodes))
        # Fancy indexing is too slow, loop is faster
        for n in np.arange(nnodes):
            total_height[n] = Z[-1, n] - Z[kbp2[n], n]

        # Average top and bottom nodes for element value
        top = good_3d_data[t, 1:, :]
        bot = good_3d_data[t, :-1, :]
        if method == 'int':
            di_data[t, :] = np.sum((top + bot) / 2.0 * elem_height, axis=0)
        else:
            height_weights = elem_height / total_height
            di_data[t, :] = np.sum((top + bot) / 2.0 * height_weights, axis=0)

    # Create NetCDF file
    if method == 'int':
        di_var = var + '_' + DEPTH_INT_DESC
    else:
        di_var = var + '_' + DEPTH_AVG_DESC
    di_filename = '{0}_{1}.61.nc'.format(stack, di_var)
    print ' - Creating {0}'.format(di_filename)
    di_file_path = os.path.join(os.path.dirname(data_file_path), di_filename)

    createDerivedNCFile(f.ncfile, var, di_file_path, di_var, di_data)
    f.ncfile.close()


def process3DFiles(data_dir, var, first_stack, last_stack, method='avg'):
    """Process 3D *.63.nc files to create depth integrated or averaged  files.

    Params:
    -------
    data_dir : string
      Path to combined model files
    var : string
      Name of variable (e.g. salt, hvel, trcr_1, etc.)
    first_stack : int
      First stack to combine (e.g. 1_*, 7_*, 100_*, etc.)
    last_stack : int
      Last stack to combine (e.g. 101_*, 4_*, 3_*, etc.)
    method : string
      Method to operate over vertical.  Either 'avg' or 'int'.
    """
    print ('Creating {0}_{1}.61.nc depth integrated files for stacks {2} to '
           '{3} in {4}'.format(var, method.upper(), first_stack, last_stack,
                               os.path.abspath(data_dir)))
    for s in np.arange(first_stack, last_stack + 1):
        data_file = '%s_%s.63.nc' % (s, var)
        data_file_path = os.path.join(data_dir, data_file)
        makeDepthFile(data_file_path, method)


#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
def parseCommandLine():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='path to combined data directory')
    parser.add_argument(
        'var', help='variable to operate on (salt, temp, etc.)')
    parser.add_argument('first_stack', type=int, help='first stack to process')
    parser.add_argument('last_stack', type=int, help='last stack to process')
    parser.add_argument(
        '-m',
        '--method',
        type=str,
        choices=[
            'avg',
            'int'],
        default='avg',
        help=("operation over vertical (default 'avg')"))

    args = parser.parse_args()

    process3DFiles(args.data_dir, args.var, args.first_stack, args.last_stack,
                   args.method)

if __name__ == '__main__':
    parseCommandLine()
