#!/usr/bin/env python
"""Create bed_flux.61.nc NetCDF files from bed_depth.61.nc files.

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
MISSING_VALUE = -9999.0
# Need to have the model output this
POROSITY = 0.5


#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
def elementToNodalValues(elem_data, connectivity):
    """Calculates and resturns nodal value from surrounding elements.

    Nodal values are average of surrounding elements.

    Params:
    -------
    elem_data - ndarray (ntimes, nElems)
      Element centered data
    connectiviy - ndarray (nElems, 3)
      Connectivity table (0-based)

    Returns:
    --------
    nodal_data - ndarray (nNodes)
      Nodal centered data
    """
    ntimes = elem_data.shape[0]
    nelems = elem_data.shape[1]
    nnodes = connectivity.max()+1

    n_data = np.zeros((ntimes, nnodes))
    n_count = np.zeros((nnodes,))

    for t in np.arange(ntimes):
        n_count = n_count*0
        for e in np.arange(nelems):
            # 3 nodes for each element
            for n in np.arange(3):
                node = connectivity[e, n]
                n_data[t, node] += elem_data[t, e]
                n_count[node] += 1
        n_data[t, :] = n_data[t, :] / n_count

    return n_data


def calculateBedFlux(bed_depth, start_depth=None):
    """Calculate and return the bed flux from changes in bed depth.

    Erosion and deposition calculated from change in bed over 1 time step.

    Params:
    -------
    bed_depth : netCDF4 file
      bed_depth file
    start_depth : ndarray(1, nNodes)
      bed_depth for time immediately before that included in bed_depth

    Returns:
    --------
    bed_flux : ndarray(nTimes-1, nNodes)
      Bed flux array
    """
    # Calculate erosion/deposition at elements
    bd = bed_depth.variables['bed_depth'][:]
    bed_flux_nodes = bd[1:, :] - bd[:-1, :]
    # Change from 1-base to 0-base
    face_nodes = bed_depth.variables['face_nodes'][:]-1
    bed_flux_elems = bed_flux_nodes[:, face_nodes].mean(axis=2)
    elem_areas = mc.computeAreas(face_nodes,
                                 bed_depth.variables['node_x'],
                                 bed_depth.variables['node_y'])
    bed_flux_elems = bed_flux_elems * elem_areas * POROSITY

    # Convert element values back to nodes for plotting and so I don't have to populate
    # a new NetCDF file.
    # TODO: Make a NetCDF file for this.
    bed_flux_nodes = elementToNodalValues(bed_flux_elems, face_nodes)

    # Try to calculate dDepth/dt else fill with nan.  Maintains array shape.
    if start_depth is not None:
        first_bed_flux = bed_depth.variables['bed_depth'][0, :] - start_depth
        # Match shape (ntimes, nnodes)
        first_bed_flux = np.expand_dims(first_bed_flux, 0)
        bed_flux = np.concatenate([first_bed_flux, bed_flux_nodes])
    else:
        fill = np.ones((1, bed_flux_nodes[0, :].shape[0]))*MISSING_VALUE
        bed_flux = np.concatenate([fill, bed_flux_nodes])

    return bed_flux


def makeBedFlux(bed_depth_path, start_depth=None):
    """Create bed_flux.61.nc NetCDF file.

    Params:
    -------
    bed_depth : string
      Path to bed_depth file
    start_depth : ndarray(1, nNodes)
      bed_depth for time immediately before that included in bed_depth

    Returns:
    --------
    last_depth : ndarray(1, nNodes)
      bed_depth for last time step
    """
    bed_depth = nc.Dataset(bed_depth_path)
    bed_flux = calculateBedFlux(bed_depth, start_depth)

    # Copy NetCDF file
    stack = bed_depth_path.split('/')[-1].split('_')[0]
    bed_flux_file = '%s_bed_flux.61.nc' % stack
    print ' - Creating %s' % bed_flux_file
    bed_flux_path = os.path.join(os.path.dirname(bed_depth_path), bed_flux_file)
    shutil.copyfile(bed_depth_path, bed_flux_path)

    # Alter copied file for deposition and erosion
    ef_file = nc.Dataset(bed_flux_path, 'a')
    ef_file.renameVariable('bed_depth', 'bed_flux')
    # (ntimes, nvrt, nnodes)
    ef = np.expand_dims(bed_flux, 1)
    ef_file.variables['bed_flux'][:] = ef
    ef_file.sync()
    ef_file.close()

    # Only return last time step
    return bed_flux[-1, :]


def processBedDepth(data_dir, first_stack, last_stack):
    """Process bed_depth files to create bed_flux files.

    Params:
    -------
    data_dir : string
      Path to combined model files
    first_stack : int
      First stack to combine (1_*, 7_*, 100_*, etc.)
    last_stack : int
      Last stack to combine (101_*, 4_*, 3_*, etc.)
    """
    print ' Creating bed_flux.61.nc files for stacks %s to %s in %s' % (
        first_stack, last_stack, os.path.abspath(data_dir))
    last_depth = None
    for s in np.arange(first_stack, last_stack+1):
        bed_depth_file = '%s_bed_depth.61.nc' % s
        bed_depth_path = os.path.join(data_dir, bed_depth_file)
        last_depth = makeBedFlux(bed_depth_path, last_depth)


#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
def parseCommandLine():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='Path to combined data directory')
    parser.add_argument('first_stack', type=int, help='First stack to processes')
    parser.add_argument('last_stack', type=int, help='Last stack to processes')

    args = parser.parse_args()

    processBedDepth(args.data_dir, args.first_stack, args.last_stack)


if __name__ == '__main__':
    parseCommandLine()
