"""
Frequently used functions.

Tuomas Karna 2015-08-21
"""
import datetime
import os
import string
from crane import plt, matplotlib
from crane.plotting.plotDataContainer import plotDataContainer
from crane.data import meshContainer
from crane.files import gmshInterface
from crane.files import gr3Interface


def readAnyMeshFile(inFile, dataFile=None):
    """
    Reads mesh data in a meshContainer.

    Supported formats are GR3 (SELFE), MSH (GMSH) and meshContainer netCDF.

    Examples:

    mc = readAnyMeshFile('mymesh.gr3')
    mc = readAnyMeshFile('myrun/data/slab/slab_salt_s0_2012-04-08_2012-04-18.nc')
    mc = readAnyMeshFile('mymesh.msh', dataFile='001_salt.msh')
    """
    fname, ext = os.path.splitext(inFile)
    if ext == '.gr3':
        mc = gr3Interface.readGR3FileToMC(inFile)
    elif ext == '.msh':
        gmsh = gmshInterface.gmshMesh.fromFile(inFile)
        if dataFile:
            mc = gmsh.getMeshContainer(fromFile=dataFile)
        else:
            mc = gmsh.getMeshContainer(fromFile=inFile)
    elif ext == '.nc':
        mc = meshContainer.meshContainer.loadFromNetCDF(inFile)
    else:
        raise Exception('Unknown file extension: ' + ext + '\n' + inFile)
    return mc


def parseTimeStr(timeStr):
    """
    Parses time string to datetime object.

    Parameters
    ----------
    timeStr : string
        string to parse. Following formats are supported:
        %Y-%m-%d-%H-%M-%S'
        %Y-%m-%d-%H-%M
        %Y-%m-%d-%H
        %Y-%m-%d

    Returns
    -------
    dt : datetime
        datetime object
    """
    if len(timeStr.split('-')) == 3:
        dt = datetime.datetime.strptime(timeStr, '%Y-%m-%d')
    elif len(timeStr.split('-')) == 4:
        dt = datetime.datetime.strptime(timeStr, '%Y-%m-%d-%H')
    elif len(timeStr.split('-')) == 5:
        dt = datetime.datetime.strptime(timeStr, '%Y-%m-%d-%H-%M')
    elif len(timeStr.split('-')) == 6:
        dt = datetime.datetime.strptime(timeStr, '%Y-%m-%d-%H-%M-%S')
    else:
        raise Exception('Could not parse date string:', timeStr)
    return dt


def createDirectory(path):
    """
    Creates given directory if it does not exist already.

    Raises an exception if a file with the same name exists.
    """
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise Exception('file with same name exists', path)
    else:
        os.makedirs(path)
    return path


def saveFigure(
        path,
        filename,
        extensions,
        verbose=False,
        dpi=200,
        bbox_tight=False):
    """
    Saves currently open figure in the given format.
    If extensions is a list, figure is saved in multiple formats

    Parameters
    ----------
    path : string
        Directory where image is stored
    filename : string
        Output file name without an extensions
    extensions : string or a list of strings
        File extension(s) that define the image file format. Must be supported by matplotlib.
    verbose : bool
        If True prints file name on stdout
    dpi : int
        Output file dot density (default 200)
    bbox_tight : bool
        If set will autocrop the figure to fit content (even if larger than original figure).
    """
    kw = {}
    if bbox_tight:
        kw['bbox_inches'] = 'tight'
    if not isinstance(extensions, list):
        extensions = [extensions]
    for ext in extensions:
        f = os.path.join(path, filename + '.' + ext)
        if verbose:
            print 'saving to', f
        plt.savefig(f, dpi=dpi, **kw)


def add_colorbar(parent_ax, p, pos='left', label=None, pad=None, aspect=None, **kwargs):
    """Adds a colorbar to parent axis without reshaping the parent axes"""
    parent_pos = parent_ax.get_position()
    fig = parent_ax.get_figure()
    parent_dx = parent_pos.x1 - parent_pos.x0
    parent_dy = parent_pos.y1 - parent_pos.y0
    vertical_pos = ['left', 'right']
    if pos in vertical_pos:
        orientation = 'vertical'
        y0 = parent_pos.y0
        dy = parent_dy
        if aspect is None:
            aspect = 1.0/6.0
        dx = aspect*dy*parent_dx
        if pos == 'right':
            if pad is None:
                pad = 0.02
            offset = pad*parent_dx
            x0 = parent_pos.x1 + offset
        else:
            if pad is None:
                pad = -0.30
            offset = pad*parent_dx
            x0 = parent_pos.x0 + offset
    else:
        orientation = 'horizontal'
        x0 = parent_pos.x0
        dx = parent_dx
        if aspect is None:
            aspect = 1.0/7.0
        dy = aspect*dx*parent_dy
        if pos == 'bottom':
            if pad is None:
                pad = -0.25
            offset = pad*parent_dy
            y0 = parent_pos.y0 + offset
        else:
            if pad is None:
                pad = 0.02
            offset = pad*parent_dy
            y0 = parent_pos.y1 + offset
    cb_ax = fig.add_axes((x0, y0, dx, dy))
    kwargs.setdefault('orientation', orientation)
    cb = plt.colorbar(p, cax=cb_ax, **kwargs)
    if pos == 'left':
        cb_ax.yaxis.tick_left()
        cb_ax.yaxis.set_label_position('right')
    if pos == 'top':
        cb_ax.xaxis.tick_top()
        cb_ax.xaxis.set_label_position('top')
    if label:
        cb.set_label(label)
    return cb, cb_ax


class AxisLabeler(object):
    """
    Adds identification characters such as a), b) etc to axes
    """
    def __init__(self, xloc=0.0, yloc=1.02, **kwargs):
        self.kwargs = kwargs
        self.kwargs.setdefault('horizontalalignment', 'left')
        self.kwargs.setdefault('verticalalignment', 'bottom')
        self.kwargs.setdefault('fontsizeincrement', 2)
        self.fontsizeincrement = self.kwargs.pop('fontsizeincrement')
        fontsize = matplotlib.rcParams['font.size']
        self.kwargs.setdefault('fontsize', fontsize + self.fontsizeincrement)
        self.xloc = xloc
        self.yloc = yloc
        self.char_iter = iter(string.ascii_lowercase)  # iterator for alphabets

    def addlabel(self, ax, char=None):
        if char is None:
            char = self.char_iter.next()
        ax.text(self.xloc, self.yloc, char + ')', transform=ax.transAxes,
                **self.kwargs)
