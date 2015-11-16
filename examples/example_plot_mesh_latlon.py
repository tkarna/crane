#!/usr/bin/python
"""

Plots grid for publications

Tuomas Karna 2015-09-18
"""
import os
from files.gmshInterface import gmshMesh
import files.gr3Interface as gr3Interface

from plotting.plotBase import *
from plotting.slabPlot import slabSnapshotDC


def readAnyMeshFile(inFile, dataFile=None):
    """Reads mesh data in a meshContainer. Supported formats are GR3 (SELFE), MSH (GMSH) and meshContainer netCDF."""
    fname, ext = os.path.splitext(inFile)
    if ext in ['.gr3', '.ll']:
        mc = gr3Interface.readGR3FileToMC(inFile)
    elif ext == '.msh':
        gmsh = gmshMesh.fromFile(inFile)
        if dataFile:
            mc = gmsh.getMeshContainer(fromFile=dataFile)
        else:
            mc = gmsh.getMeshContainer(fromFile=inFile)
    elif ext == '.nc':
        mc = meshContainer.loadFromNetCDF(inFile)
    else:
        raise Exception('Unknown file extension: ' + ext + '/n' + inFile)
    return mc


def addChar(ax, char=None):
    if char:
        ax.text(-0.04, 1.01, char + ')', fontsize=fontsize + 0,
                verticalalignment='bottom', horizontalalignment='left',
                transform=ax.transAxes)

# --------------------------------------------------------------

# tune font size
fontsize = 14
matplotlib.rcParams['font.size'] = fontsize

# read mesh
mc = readAnyMeshFile('/home/workspace/ccalmr/hindcasts/db33/hgrid.gr3')

imgDir = 'tmp'
fPrefix = 'mesh'
filetype = 'png'
imgName = 'combo'
createDirectory(imgDir)

# bounding box for image
# this is in the same coordinates as the mesh
#bBox = [-50000, 490000, -500000, 740000]  # whole domain
bBox = [330000, 403000, 275000, 312000]  # estuary
mc = mc.cropGrid(bBox)

# figure size is 15 by 15 inch
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111)

coordSys = 'latlon'  # will plot lat/lon coords
#coordSys = 'spcs'    # will plot distances in km

dia = slabSnapshotDC(clabel='', unit='-', coordSys=coordSys)
dia.setAxes(ax)
dia.addSample(mc, 0, bbox=bBox, plotType='mesh', linewidth=0.15, color='k')
# add title (if needed)
dia.addTitle('Hello I\'m Mojy, I\'m from Iran')

# fine tune map bounds
# with latlon bounds may be too large
# None will plot all in the bounding box
ylim = None  # [46.08, 46.36]
xlim = None

# tune number of ticks of x/y axis (if neeeded)
dia.updateXAxis(maxticks=12, xlim=xlim)
dia.updateYAxis(maxticks=3, ylim=ylim)

# save image
# bbox_tight will crop the figure to content
file = 'my_example_grid_figure'
saveFigure(imgDir, file, filetype, verbose=True, dpi=200, bbox_tight=True)
plt.close()
