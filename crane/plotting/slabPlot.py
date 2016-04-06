#!/usr/bin/python

"""
A class for plotting horizontal slices of 3D data (slabs).

Tuomas Karna 2012-11-27
"""

import numpy as np
import datetime
from crane.data import timeArray
from crane.plotting import plotBase
from matplotlib.tri import Triangulation
from crane.data import coordSys as coordSysModule
from scipy.spatial import cKDTree


def generateTriangulation(x, y, connectivity):
    """Return a Triangulation from a connectivity array"""
    return Triangulation(x.copy(), y.copy(), connectivity.copy())


def transformBBox(bbox, trans):
    """Apply coordinate transformation on bounding box"""
    new_bbox = [0, 0, 0, 0]
    new_bbox[0], new_bbox[2] = trans.transform(bbox[0], bbox[2])
    new_bbox[1], new_bbox[3] = trans.transform(bbox[1], bbox[3])
    return new_bbox


def convertTriCoords(tri, bbox=None, scalar=1, coordSys='spcs',
                     flipXCoord=False):
    if not bbox:
        bbox = [tri.x.min(), tri.x.max(), tri.y.min(), tri.y.max()]

    if coordSys.lower() == 'none':
        return tri, unitTransform, bbox
    elif coordSys == 'latlon':
        trans = latlonTransformer(bbox[0], bbox[2], scalar, flipXCoord)
    elif coordSys == 'utm':
        trans = utmTransformer(bbox[0], bbox[2], scalar)
    else:
        trans = coordTransformer(bbox[0], bbox[2], scalar)

    new_x, new_y = trans.transform(tri.x, tri.y)
    new_tri = generateTriangulation(new_x, new_y, tri.triangles)
    new_tri.set_mask(tri.mask)
    new_bbox = transformBBox(bbox, trans)
    return new_tri, trans, new_bbox


class coordTransformer(object):
    """Simple linear coordinate transformer class for slabs."""

    def __init__(self, biasX, biasY, scale):
        self.biasX = biasX
        self.biasY = biasY
        self.scale = scale

    def transform(self, x, y):
        newX = (x - self.biasX) * self.scale
        newY = (y - self.biasY) * self.scale
        return newX, newY


class latlonTransformer(object):
    """Transforms state plane coordinates to lat lon"""

    def __init__(self, biasX, biasY, scale, flipXCoord):
        self.biasX = biasX
        self.biasY = biasY
        self.scale = scale
        self.flipXCoord = flipXCoord

    def transform(self, x, y):
        newX, newY = coordSysModule.spcs2lonlat(x, y)
        if self.flipXCoord:
            newX *= -1.0
        return newX * self.scale, newY * self.scale


class utmTransformer(object):
    """Transforms state plane coordinates to UTM coordinates"""

    def __init__(self, biasX, biasY, scale):
        self.biasX = biasX
        self.biasY = biasY
        self.scale = scale

    def transform(self, x, y):
        newX, newY = coordSysModule.spcs2utm(x, y)
        return newX * self.scale, newY * self.scale

unitTransform = coordTransformer(0.0, 0.0, 1.0)


class slabSnapshot(plotBase.colorPlotBase):
    """Plots a horizontal map of a triangular mesh"""

    def __init__(self, **kwargs):
        defaultArgs = {}
        defaultArgs.update(kwargs)
        # default plot options for all diagrams
        self.unit = defaultArgs.pop('unit', '-')
        self.coordSys = defaultArgs.pop('coordSys', 'spcs')
        self.clabel = defaultArgs.pop('clabel', '')
        self.flipXCoord = defaultArgs.pop('flipXCoord', False)
        if self.coordSys == 'latlon':
            xunit = u'\u00B0E'
            if self.flipXCoord:
                xunit = u'\u00B0W'
            yunit = u'\u00B0N'
            xlabel = 'Longitude'
            ylabel = 'Latitude'
            coordScalar = 1.0
            aspect = 1.3
        elif self.coordSys == 'utm':
            xunit = u'kmW'
            yunit = u'kmN'
            xlabel = 'UTM 10 N'
            ylabel = 'UTM 10 N'
            coordScalar = 1e-3
            aspect = 1.0
        elif self.coordSys == 'spcs':
            xunit = u'km'
            yunit = u'km'
            xlabel = 'x'
            ylabel = 'y'
            coordScalar = 1e-3
            aspect = 1.0
        self.coordScalar = defaultArgs.pop('coordscalar', coordScalar)
        self.aspect = defaultArgs.pop('aspect', aspect)
        self.xlabel = defaultArgs.pop('xlabel', xlabel)
        self.ylabel = defaultArgs.pop('ylabel', ylabel)
        self.xunit = defaultArgs.pop('xunit', xunit)
        self.yunit = defaultArgs.pop('yunit', yunit)
        defaultArgs.setdefault('clim', [])
        defaultArgs.setdefault('climIsLog', False)
        self.plotType = defaultArgs.pop('plotType', 'color')
        self.logScale = defaultArgs.pop('logScale', False)
        if self.logScale:
            self.unit = 'log10(' + self.unit + ')'

        super(slabSnapshot, self).__init__(**defaultArgs)
        self.ax = None
        self.cax = None
        self.coordTransformer = None

    def setAxes(self, ax):
        """Set axes for the diagram. All data will be plotted in these axes.
        """
        self.ax = ax
        # ax.grid()
        yStr = self.ylabel
        if self.yunit:
            yStr += ' [' + self.yunit + ']'
        ax.set_ylabel(yStr)
        xStr = self.xlabel
        if self.xunit:
            xStr += ' [' + self.xunit + ']'
        ax.set_xlabel(xStr)
        if self.coordSys == 'latlon':
            # prevent constant notation in ticklabels
            ax.ticklabel_format(axis='both', useOffset=False)

    def prepDataForPlotting(self, tri, variable=None, **kwargs):
        """Generates data array that to be sent to tripcolor/tricontour function"""
        kw = {}
        kw.update(self.defaultArgs)
        kw.update(kwargs)
        clim = kw.pop('clim', [])
        climIsLog = kw.pop('climIsLog', False)

        kw.pop('clabel', None)
        kw.pop('unit', None)
        kw.pop('logScale', None)

        plotType = kw.get('plotType', 'color')

        if variable is not None:
            data = variable.copy() if not self.logScale else np.log10(variable)
            # data for each element instead of each node
            dataByElems = data is not None and len(
                data) == tri.triangles.shape[0]
            if dataByElems:
                plotType = 'color'  # only one that supports element-wise data
                kw.pop('N', 20)
                tri.set_mask(np.isnan(data))
            else:
                tri.set_mask(np.max(np.isnan(data[tri.triangles]), axis=1))
            if clim:
                cmin, cmax = clim
                if self.logScale and not climIsLog:
                    cmin, cmax = np.log10(cmin), np.log10(cmax)

                # saturate limits
                threshold = 1e-9
                goodIx = np.logical_not(np.isnan(data))
                data[
                    np.logical_and(
                        ~np.isnan(data),
                        data < cmin)] = cmin + threshold
                data[
                    np.logical_and(
                        ~np.isnan(data),
                        data > cmax)] = cmax - threshold
            else:
                cmin, cmax = np.nanmin(data), np.nanmax(data) + 1e-9
            if plotType in ['contour', 'contourf']:
                N = kw.pop('N', 20)
                kw.setdefault('levels', np.linspace(cmin, cmax, N))
            if plotType == 'color':
                kw.setdefault('vmin', cmin)
                kw.setdefault('vmax', cmax)
            if dataByElems:
                kw['facecolors'] = data
                data = None
            kw['dataByElems'] = dataByElems
        else:
            data = None
            kw['dataByElems'] = False
        return data, kw

    def addSample(self, tri, variable=None, **kwargs):
        """Add scalar field to the diagram"""
        if kwargs.get('plotType') == 'quiver':
            self.plot_quiver(tri, variable, **kwargs)
            return

        # NOTE only plots the first component
        data, kw = self.prepDataForPlotting(tri, variable[:, 0], **kwargs)

        plotType = kw.pop('plotType', 'color')
        dataByElems = kw.pop('dataByElems', False)
        bbox = kw.pop('bbox', [])
        if self.coordSys == 'none':
            tri_trans, trans, convertedBBox = tri, unitTransform, bbox
        else:
            tri_trans, trans, convertedBBox = convertTriCoords(
                tri, bbox, scalar=self.coordScalar, coordSys=self.coordSys,
                flipXCoord=self.flipXCoord)

        self.coordTransformer = trans

        if plotType == 'mesh':
            kw.setdefault('linewidth', 0.3)
            self.ax.triplot(tri_trans, **kw)
        elif plotType == 'contour':
            self.cax = self.ax.tricontour(tri_trans, data, **kw)
        elif plotType == 'contourf':
            self.cax = self.ax.tricontourf(tri_trans, data, **kw)
        else:  # color
            if not dataByElems:
                kw.setdefault('shading', 'gouraud')
            kw.pop('coordsys', None)
            self.cax = self.ax.tripcolor(tri_trans, data, **kw)
        self.ax.set_aspect(self.aspect)
        if bbox:
            self.ax.set_xlim(convertedBBox[:2])
            self.ax.set_ylim(convertedBBox[2:])

    def plot_quiver(self, tri, var, **kwargs):
        """Does a 2D quiver vector plot.
        Must first plot a scalar field on the diagram.

        Quiver specific kwargs:
        nquiverpoints : int
            number of points to draw along the longest edge of the plot
        maxmagnitude : float
            draw only vectors whose lenght is less than this value
        scale : float
            define arrow scaling factor, smaller value yields longer arrows
        """
        assert var.shape[1] == 2, 'quiver plot requires 2d vector data'

        kw = dict(self.defaultArgs)
        kw.update(kwargs)

        bbox = kw.pop('bbox', None)
        if self.coordTransformer is None:
            tri_trans, trans, convertedBBox = convertTriCoords(
                tri, bbox, scalar=self.coordScalar, coordSys=self.coordSys,
                flipXCoord=self.flipXCoord)
            self.coordTransformer = trans
        else:
            if bbox is not None:
                convertedBBox = transformBBox(bbox, self.coordTransformer)

        kwargs.pop('plotType')
        npoints = int(kwargs.pop('nquiverpoints', 120))
        maxmag = float(kwargs.pop('maxmagnitude', 100.0))
        kwargs.setdefault('units', 'dots')
        kwargs.setdefault('width', 1.7)
        kwargs.setdefault('headlength', 6)
        # get x,y,u,v arrays
        uv_x = tri.x.ravel()
        uv_y = tri.y.ravel()
        uv_x, uv_y = self.coordTransformer.transform(uv_x, uv_y)
        u = var[:, 0]
        v = var[:, 1]
        # filter too long vectors
        mag = np.hypot(u, v)
        ix = mag < maxmag
        u = u[ix]
        v = v[ix]
        uv_x = uv_x[ix]
        uv_y = uv_y[ix]
        # generate regular grid that spans the plot
        xlim = uv_x.min(), uv_x.max()
        ylim = uv_y.min(), uv_y.max()
        if bbox is not None:
            xlim = max(convertedBBox[0], xlim[0]), min(convertedBBox[1], xlim[1])
            ylim = max(convertedBBox[2], ylim[0]), min(convertedBBox[3], ylim[1])
        lx = xlim[1] - xlim[0]
        ly = ylim[1] - ylim[0]
        aspect = ly/lx
        if lx > ly:
            nx = npoints
            ny = aspect*npoints
        else:
            nx = npoints/aspect
            ny = npoints
        xgrid = np.linspace(xlim[0], xlim[1], nx)
        ygrid = np.linspace(ylim[0], ylim[1], ny)
        # filter by taking nearest nodes to the regular grid
        tree = cKDTree(np.concatenate((uv_x[:, None], uv_y[:, None]), axis=1))
        xx, yy = np.meshgrid(xgrid, ygrid)
        p = np.vstack((xx.ravel(), yy.ravel())).T
        nearestnodes = tree.query(p)[1]
        ix = np.unique(nearestnodes)
        q = self.ax.quiver(uv_x[ix], uv_y[ix], u[ix], v[ix], **kwargs)
        self.quiverobj = q

        if bbox is not None:
            self.ax.set_xlim(convertedBBox[:2])
            self.ax.set_ylim(convertedBBox[2:])
        return q

    def addQuiverLegend(self, x, y, u, label, **kwargs):
        """Adds a quiver legend to plot at location (x,y) and lenght u"""
        assert hasattr(self, 'quiverobj'), 'quiver object has not been created'
        coords = kwargs.pop('coordinates', 'axes')
        if coords == 'data':
            xx, yy = self.coordTransformer.transform(x, y)
        else:
            xx, yy = x, y
        self.ax.quiverkey(self.quiverobj, xx, yy, u, label, **kwargs)

    def updatePlot(self, tri, variable=None, **kwargs):
        """Replaces the plot with new data without re-drawing the whole axes."""
        # As dry mask changes in time, essentially requires a re-draw of the
        # mesh
        if self.ax is not None:
            self.ax.collections.remove(self.cax)
            self.addSample(tri, variable, **kwargs)

    def updateColorBar(self):
        """Updates the color bar. This is meant to be called when the associated
           plot has been updated trough updatePlot routine."""
        raise NotImplementedError(
            'updateColorBar not implemented for slab plots')

    def addTransectMarker(self, x, y, label=None, textkw={}, **kwargs):
        """Adds a line to the plot, defined by the coordinates.
        If label is given, adds text to the end point. kwargs are passed to ax.plot command."""
        # convert coordinates
        x_tr, y_tr = self.coordTransformer.transform(x, y)
        # defargs
        kwargs.setdefault('color', 'k')
        kwargs.setdefault('linewidth', 1.3)
        self.ax.plot(x_tr, y_tr, **kwargs)
        if label is not None:
            tkw = {}
            tkw['color'] = kwargs['color']
            tkw['verticalalignment'] = 'bottom'
            tkw['horizontalalignment'] = 'left'
            tkw.update(textkw)
            off = min(0.01 * (self.ax.dataLim.xmax - self.ax.dataLim.xmin),
                      0.01 * (self.ax.get_xlim()[1] - self.ax.get_xlim()[0]))
            self.ax.text(x_tr[-1] + off, y_tr[-1] + off, label, **tkw)

    def addStationMarker(self, x, y, printLabel=False, textkw={}, **kwargs):
        """Adds a marker in the given position with a label text.
        kwargs are passed to ax.plot routine."""
        # convert coordinates
        x, y = self.coordTransformer.transform(x, y)
        # defargs
        kwargs.setdefault('linestyle', 'none')
        kwargs.setdefault('markersize', 8)
        kwargs.setdefault('color', 'k')
        kwargs.setdefault('marker', 'o')
        xoffset = kwargs.pop('xoffset', 0.01)
        yoffset = kwargs.pop('yoffset', 0.01)
        txt = None
        if printLabel:
            label = kwargs.pop('label')
            xoff = min(
                xoffset * (self.ax.dataLim.xmax - self.ax.dataLim.xmin),
                xoffset * (self.ax.get_xlim()[1] - self.ax.get_xlim()[0]))
            yoff = min(
                yoffset * (self.ax.dataLim.ymax - self.ax.dataLim.ymin),
                yoffset * (self.ax.get_ylim()[1] - self.ax.get_ylim()[0]))
            tkw = {}
            tkw['color'] = kwargs['color']
            tkw['verticalalignment'] = 'bottom'
            tkw['horizontalalignment'] = 'left'
            tkw.update(textkw)
            if tkw['horizontalalignment'] == 'right':
                xoff *= -1
            if tkw['verticalalignment'] == 'top':
                yoff *= -1
            txt = self.ax.text(x + xoff, y + yoff, label, **tkw)
        marker = self.ax.plot(x, y, **kwargs)
        return marker, txt

    def addStationLegend(self, **kwargs):
        kw = {}
        kw['prop'] = {'size': 12}
        kw['numpoints'] = 1
        if 'loc' not in kwargs and 'bbox_to_anchor' not in kwargs:
            kw['loc'] = 'upper left'
        #kw['bbox_to_anchor'] = (-0.2, 1.0)
        kw['ncol'] = 1
        kw.update(kwargs)
        self.ax.legend(**kw)

    def calcColorTicks(self, tri, **kwargs):
        """Calculates a reasonable amount of ticks for colorbar based on slab w/h ratio"""
        height = tri.y.max - tri.y.min
        width = tri.x.max - tri.x.min
        whRatio = float(height) / width


def generateSlabFromMeshContainer(mc, timeStamp):
    """Creates triangulation of the grid and retrieves data at correct time.
    """
    if isinstance(timeStamp, datetime.datetime):
        # interpolate to correct time
        t_epoch = timeArray.datetimeToEpochTime(timeStamp)
        newTime = timeArray.timeArray(np.array([t_epoch]), 'epoch')
        mc = mc.interpolateInTime(newTime, acceptNaNs=True)
        x = mc.x.flatten()
        y = mc.y.flatten()
        z = mc.z.flatten()
        data = mc.data[:, :, 0]  # shape (npts,nComp)
        time = mc.time.getDatetime(0)
    else:
        x = mc.x[:, timeStamp] if mc.xDependsOnTime else mc.x
        y = mc.y[:, timeStamp] if mc.yDependsOnTime else mc.y
        z = mc.z[:, timeStamp] if mc.zDependsOnTime else mc.z
        data = mc.data[:, :, timeStamp]  # shape (npts,nComp)
        time = mc.time.getDatetime(timeStamp)

    triang = generateTriangulation(x, y, mc.connectivity)
    return triang, data, time


class slabSnapshotDC(slabSnapshot):
    """Slab plot class that uses dataContainer as input."""

    def __init__(self, **defaultArgs):
        slabSnapshot.__init__(self, **defaultArgs)

    def addSample(self, mc, timeStamp, **kwargs):
        """Add slab from dataContainer wih appropriate data restructuring.
           Args:
           mc        -- (meshContainer) data
           timeStamp -- (int or datetime) time index/stamp to plot.
                        If datetime object, temporal interpolation is performed.
           kwargs    -- (**) passed to slabSnapshot.addSample
        """
        tri, var, time = generateSlabFromMeshContainer(mc, timeStamp)
        slabSnapshot.addSample(self, tri, var, **kwargs)
        titleStr = mc.description.split(
            '.')[0] + ' ' + time.strftime('%Y-%m-%d %H:%M:%S') + ' (PST)'
        titleStr = kwargs.pop('title', titleStr)
        self.addTitle(titleStr)

    def updatePlot(self, mc, timeStamp, **kwargs):
        """Updates the most recent plotted data array with new values."""
        triang, var, time = generateSlabFromMeshContainer(mc, timeStamp)
        slabSnapshot.updatePlot(self, triang, var, **kwargs)
        titleStr = mc.description.split(
            '.')[0] + ' ' + time.strftime('%Y-%m-%d %H:%M:%S') + ' (PST)'
        titleStr = kwargs.pop('title', titleStr)
        self.updateTitle(titleStr)


class stackSlabPlot(plotBase.stackPlotBase):
    """A class for stacking multiple slabs in the same plot."""

    def addPlot(self, tag, **kwargs):
        """Adds a new subplot to the diagram"""
        kw = {}
        kw.update(self.defArgs)
        kw.update(kwargs)
        plot = slabSnapshot(**kw)
        super(stackSlabPlot, self).addPlot(plot, tag)

    def addSample(self, tag, tri, variable, **kwargs):
        self.plots[tag].addSample(tri, variable, **kwargs)
        if kwargs.get('draw_cbar', True):
            self.plots[tag].showColorBar()

    def updatePlot(self, tag, tri, variable, **kwargs):
        self.plots[tag].updatePlot(tri, variable, **kwargs)

    def addTransectMarker(self, tag, x, y, label=None, **kwargs):
        """Adds a marker to the plot, defined by the coordinates.
        If label is given, adds text to the end point. kwargs are passed to ax.plot command."""
        if tag != 'all':
            self.plots[tag].addStationMarker(x, y, label, **kwargs)
        else:
            for t in self.plots:  # recursion
                self.addStationMarker(t, x, y, label, **kwargs)

    def addStationMarker(self, tag, x, y, label, textkw={}, **kwargs):
        """Adds a filled circle in the given position with a label text.
        If tag is 'all', marker is added to all subplots.
        kwargs are passed to ax.plot routine."""
        if tag != 'all':
            self.plots[tag].addStationMarker(x, y, label, **kwargs)
        else:
            for t in self.plots:  # recursion
                self.addStationMarker(t, x, y, label, **kwargs)


class stackSlabPlotDC(stackSlabPlot):
    """stackSlabPlot class that uses dataContainer as an input"""

    def __init__(self, **defArgs):
        stackSlabPlot.__init__(self, **defArgs)

    def addSample(self, tag, mc, timeStamp, **kwargs):
        """Add transect from meshContainer wih appropriate data restructuring.
           Args:
           tag       -- (string) tag for identifiying the subplot
           mc        -- (meshContainer) data
           timeStamp -- (int or datetime) time index/stamp to plot.
                        If datetime object, temporal interpolation is performed.
           kwargs    -- (**) passed to transectSnapshot.addSample
        """
        triang, var, time = generateSlabFromMeshContainer(mc, timeStamp)
        stackSlabPlot.addSample(self, tag, triang, var, **kwargs)

    def updatePlot(self, tag, mc, timeStamp, **kwargs):
        triang, var, time = generateSlabFromMeshContainer(mc, timeStamp)
        stackSlabPlot.updatePlot(self, tag, triang, var, **kwargs)

    def addTransectMarker(self, tag, x, y, label=None, **kwargs):
        """Adds a marker to the plot, defined by the coordinates.
        If label is given, adds text to the end point. kwargs are passed to ax.plot command."""
        if tag != 'all':
            self.plots[tag].addTransectMarker(x, y, label, **kwargs)
        else:
            for t in self.plots:  # recursion
                self.addTransectMarker(t, x, y, label, **kwargs)

    def addStationMarker(self, tag, x, y, label, textkw={}, **kwargs):
        """Adds a filled circle in the given position with a label text.
        If tag is 'all', marker is added to all subplots.
        kwargs are passed to ax.plot routine."""
        kwargs.setdefault('label', label)
        kwargs.setdefault('printLabel', True)
        if tag != 'all':
            self.plots[tag].addStationMarker(x, y, **kwargs)
        else:
            for t in self.plots:  # recursion
                self.addStationMarker(t, x, y, **kwargs)


if __name__ == '__main__':
    from crane import plt
    # generate dummy data
    # First create the x and y coordinates of the points.
    n_angles = 36
    n_radii = 8
    min_radius = 0.25
    radii = np.linspace(min_radius, 0.95, n_radii)

    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    angles = np.repeat(angles[..., np.newaxis], n_radii, axis=1)
    angles[:, 1::2] += np.pi / n_angles

    x = (radii * np.cos(angles)).flatten()
    y = (radii * np.sin(angles)).flatten()
    z = (np.cos(radii) * np.cos(angles * 3.0)).flatten()

    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triang = tri.Triangulation(x, y)

    # Mask off unwanted triangles.
    xmid = x[triang.triangles].mean(axis=1)
    ymid = y[triang.triangles].mean(axis=1)
    mask = np.where(xmid * xmid + ymid * ymid < min_radius * min_radius, 1, 0)
    triang.set_mask(mask)

    fig = plt.figure()
    dia = slabSnapshot(clabel='Salinity', unit='psu')
    dia.setupAxes(fig)
    dia.addSample(triang, z, linewidth=0.3, edgecolors=[0.5, 0.5, 0.5])
    dia.addSample(triang, plotType='mesh')
    dia.showColorBar()
    plt.show()

    # example with meshContainer
    from crane.data import meshContainer
    d0 = meshContainer.meshContainer.loadFromNetCDF(
        '/home/tuomas/workspace/cmop/utils/processing/data/tmp/test_hvel_400_2010-06-14_2010-06-15.nc')
    d0 = d0.cropGrid([4000, 1e5, 200, 1e5])
    print d0
    print d0.x.shape, d0.y.shape, d0.connectivity.shape, d0.data[:, 0, -1].shape
    print d0.connectivity.min(), d0.connectivity.max()

    fig = plt.figure()
    dia = slabSnapshotDC(clabel='Salinity', unit='psu')
    dia.setupAxes(fig)
    dia.addSample(d0, 20, linewidth=0.3, edgecolors='b')
    dia.setupAxes(fig)
    dia.addSample(d0, 20, linewidth=0.3, edgecolors='b')
    dia.showColorBar()
    plt.show()
