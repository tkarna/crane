#!/usr/bin/env python
"""Create a depth intregtated NetCDF *.61.nc files from *.63.nc files.

Jesse Lopez
"""
# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
import os
import shutil
import argparse
import glob

import numpy as np
import netCDF4 as nc

import crane.data.ncExtract as nce
import crane.data.dirTreeManager as dtm
import crane.data.meshContainer as mc
from crane.physicalVariableDefs import NODE_FILES, ELEM_FILES, WHOLE_LEVELS, SCALAR_FILES, VECTOR_FILES, FILES_2D


# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
DEPTH_INT_DESC = 'DI'
DEPTH_AVG_DESC = 'DAVG'
FILL_VALUE = -99.0


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------
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
                if len(old.dimensions['node']) == data.shape[-1]:
                    new_dims = (u'time', u'node')
                elif len(old.dimensions['face']) == data.shape[-1]:
                    new_dims = (u'time', u'face')
                else:
                    new_dims = (u'time', u'edge')
                out_var = new.createVariable(new_var, var.dtype, new_dims)
            # Push the data to the file
            out_var[:] = data[:]
            new.sync()
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

    data_dir = os.path.dirname(data_file_path)
    filename = os.path.basename(data_file_path)
    split = filename.split('.')[0].split('_')
    stack = split[0]
    if len(split[1:]) > 1:
        var = '_'.join(split[1:])
    else:
        var = split[1]
    ftype = filename.split('.')[1]

    f.readHeader()

    data = f.ncfile.variables[var][:]
    elev = f.ncfile.variables['elev'][:]
    bath = f.ncfile.variables['depth'][:]
    di_data = np.zeros_like(elev)
    ntimes = f.ncfile.variables['time'][:].shape[0]
    npoints = f.ncfile.variables['elev'][:].shape[1]

    # Only integrate non-dry values (not the -99.0 fill value)
    good_3d_data = np.ma.masked_equal(data, FILL_VALUE)

    # Operate over values for every time step
    for t in np.arange(ntimes):
        # Get height of prisms
        Z, kbp2, iwet = f.vCoords.computeVerticalCoordinates(elev[t, :], bath)
        z_top = Z[1:, :]
        z_bot = Z[:-1, :]
        elem_height = z_top - z_bot
        total_height = np.ma.empty((npoints))
        # Fancy indexing is too slow, loop is faster
        for n in np.arange(npoints):
            total_height[n] = Z[-1, n] - Z[kbp2[n], n]

        if ftype in WHOLE_LEVELS:
            # Average top and bottom nodes for element value
            top = good_3d_data[t, 1:, :]
            bot = good_3d_data[t, :-1, :]
            if method == 'int':
                di_data[t, :] = np.sum((top + bot) / 2.0 * elem_height, axis=0)
            else:
                height_weights = elem_height / total_height
                di_data[t, :] = np.sum((top + bot) / 2.0 * height_weights, axis=0)
        # Bottom vertical levels are duplicated in combine and skipped here
        else:
            if method == 'int':
                di_data[t, :] = np.sum(good_3d_data[t, 1::, :] * elem_height, axis=0)
            else:
                height_weights = elem_height / total_height
                di_data[t, :] = np.sum(good_3d_data[t, 1::, :] * height_weights, axis=0)

    # Create NetCDF file
    if method == 'int':
        di_var = var + '_' + DEPTH_INT_DESC
    else:
        di_var = var + '_' + DEPTH_AVG_DESC
    if ftype in NODE_FILES:
        if ftype in SCALAR_FILES:
            new_ftype = '61'
        else:
            new_ftype = '62'
    elif ftype in ELEM_FILES:
        new_ftype = '66'
    else:
        if ftype in SCALAR_FILES:
            new_ftype = '65'
        else:
            new_ftype = '71'
    di_filename = '{0}_{1}.{2}.nc'.format(stack, di_var, new_ftype)
    print ' - Creating {0} from {1}'.format(di_filename, filename)
    di_file_path = os.path.join(os.path.dirname(data_file_path), di_filename)

    createDerivedNCFile(f.ncfile, var, di_file_path, di_var, di_data)
    f.ncfile.close()


def process3DFiles(data_dir, var, first_stack, last_stack, method='avg'):
    """Process 3D *.nc files to create depth integrated or averaged files.

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
    print('Creating depth integrated files for stacks {0} to '
          '{1} in {2}'.format(first_stack, last_stack,
                              os.path.abspath(data_dir)))

    for s in np.arange(first_stack, last_stack + 1):
        g = os.path.join(data_dir, '%s_%s.*' % (s, var))
        for f in glob.glob(g):
            ftype = f.split('.')[1]
            if ftype in FILES_2D or ftype in VECTOR_FILES:
                print('Cannot operate on 2D or vector files: {0}'.format(f))
                continue
            makeDepthFile(f, method)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def parseCommandLine():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='path to combined data directory')
    parser.add_argument('var', help='variable to operate on (salt, temp, etc.)')
    parser.add_argument('first_stack', type=int, help='first stack to process')
    parser.add_argument('last_stack', type=int, help='last stack to process')
    parser.add_argument('-m', '--method', type=str, choices=['avg', 'int'],
                        default='avg',
                        help=("operation over vertical (default 'avg')"))
    args = parser.parse_args()

    process3DFiles(args.data_dir, args.var, args.first_stack, args.last_stack,
                   args.method)

if __name__ == '__main__':
    parseCommandLine()
