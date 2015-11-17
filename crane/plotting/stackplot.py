"""
NOTE this module is now OBSOLETE and should be removed ASAP
NOTE time series plots should be handled by timeSeriesPlotDC2 derived methods
"""
import matplotlib

import numpy as np
import numpy.ma as ma

from crane import plt
from matplotlib import cm
from matplotlib import dates

from pylab import find

import pylabutil


def getvalidplotargs():
    return [
        'alpha',
        'animated',
        'antialiased',
        'axes',
        'clip_box',
        'clip_on',
        'clip_path',
        'color',
        'contains',
        'dash_capstyle',
        'dash_joinstyle',
        'dashes',
        'data',
        'drawstyle',
        'figure',
        'fillstyle',
        'gid',
        'label',
        'linestyle',
        'linewidth',
        'lod',
        'marker',
        'markeredgecolor',
        'markeredgewidth',
        'markerfacecolor',
        'markersize',
        'markevery',
        'picker',
        'pickradius',
        'rasterized',
        'snap',
        'solid_capstyle',
        'solid_joinstyle',
        'transform',
        'url',
        'visible',
        'xdata',
        'ydata',
        'zorder']


def addGaps(t, arr, **kwargs):
    '''Take a timeseries of data with presumed even timestep.
    Find gaps in the data that are longer than
     the timestep * gapfactor, where timestep can be given as an argument
    or calculated from the data as median difference of t.
    A NaN value is inserted between the to sides of the gap to cause a
    gap in a timeseries plot. The data can contain NaN values initially.
     NaN values are stripped out before the calculation of the timestep.
     Optional arguments: timestep , gapfactor
    '''
    timestep = kwargs.get('timestep', None)
    factor = kwargs.get('gapfactor', 5)
    t = ma.array(t, mask=np.equal(t, None))
    arr = ma.array(arr, mask=np.equal(arr, None) | np.equal(arr, np.nan))
    # ensure that the arrays are masked arrays, with None values masked
    # find the overlap of the masks of both the time and the data
    t.mask = arr.mask | t.mask
    arr.mask = t.mask  # set both arrays to use the overlap of the masks
    t = t.compressed()  # compress the time to remove the masked values
    a = ma.array(arr.compressed(), mask=False)  # compress the data
    dt = np.diff(t)  # find the actual timesteps between data values
    if timestep is None:  # if no user supplied expected timestep was supplied
        timestep = np.median(dt)  # calculate the median timestep
    # find gaps larger than some factor times the expected timestep
    gaps = find(np.abs(dt) > timestep * factor)
    a = np.insert(a, gaps + 1, np.nan)   # insert np.nan values in gaps
    t = np.insert(t, gaps + 1, np.nan)
    return (t, a)


def ismasked(x):
    '''Return True if the argument is a MaskedArray (via duck typing), False otherwise'''
    return hasattr(x, 'mask')


def maskNulls(arr, maskval=None):
    '''Return a masked array masked on values = mask, extending the mask if necessary.'''

    if ismasked(arr):  # retain existing mask
        arr.mask = np.equal(np.array(arr, dtype='object'), maskval) | arr.mask
        return arr
    else:
        return ma.masked_array(
            arr,
            mask=np.equal(
                np.array(
                    arr,
                    dtype='object'),
                maskval))


def GetTZCorrection():
    """ Returns -8*60*60 (8 hours, 60 minutes, 60 seconds) """
    return -8.0 * 3600.0


class dataset:

    def __init__(self, lbl, X=None, Y=None, Z=None, **kwargs):
        '''
        additional arguments include:
           xlim         : (low x, high x)
           ylim         : (low y, high y)
           zlim         : (low z, high z)
           zlabel: label for a colorbar for the 3rd dimension
           color: color for plotting dataset
        '''
        self.label = lbl
        self.kwargs = kwargs
        try:
            X = np.array([todate(t) for t in X])
            self.x_is_time = True
        except:
            self.x_is_time = False
            pass
        self.X = X
        self.Y = Y
        self.Z = Z
        if self.X is not None:
            self.X = maskNulls(self.X)
        if self.Y is not None:
            self.Y = maskNulls(self.Y)
        if self.Z is not None:
            self.Z = maskNulls(self.Z)


def stackplot(data, **kwargs):
    '''
  Creates a matplotlib figure (via pylab) displaying a vertical stack of
  scatter plots with a shared x-axis.

  Arguments:
    data     : a list of (yaxis, overlays) pairs,
               where overlays is a list of (label, X, Y) pairs, and
               label and yaxis are strings and array is a numpy numeric array
               Examples:
                  [('salinity', [dataset('my data', aX, asY,**kwargs), dataset('your data', bX, bsY, **kwargs)])]
                  [
                   ('salinity', [
                                 dataset('my data', aX, atY),
                                 dataset('your data', bX, btY)
                                ]
                   ),
                   ('temperature', [
                                 dataset('my data', aX, atY),
                                 dataset('your data', bX, btY)
                                ]
                   ),
                  ]

  Optional keyword arguments are:
    xlabel       : label to place on the x-axis.
    x_is_time    : Interpret x-axis as time strings
    dateformat   : xtick format for time axis,
                   defaults to '%m-%d %H:%M'
    sourceformat : source format for time string data,
                   defaults to '%Y-%m-%d %H:%M:%S'
    color_order  : order to rotate through colors, defaults to 'rgbcmykw'
    width        : width of the figure, defaults to 6in
    title        : title of the figure
    nolegend     : suppress the legend
    xlim         : (low x, high x)
    ylim         : (low y, high y)
    style        : style of plot marker {. - -- : , o }  plus others
    plottype     : 'scatter' or '', defaults to ''
    handlegaps   : useful for timeseries non-scatter plots with a line style,
                   defaults to False
    timestep     : used with handlegaps, expected time step,
                   defaults to the median timestep of the data
    gapfactor    : used with handlegaps, this times the timestep determines
                   the minimum size of gap to handle,
                   defaults to 5
  (see http://matplotlib.sourceforge.net/matplotlib.pyplot.html#-plot)
  '''
    #cmop.debug('starting stacked plot')
    # prepare the figure layout values
    n = len(data)
    w_scale = 1.0
    topmargin = 0.15
    plotheight = 1.5
    max_legend_width = 0.2  # maximum legend width needs to be carried across all subplots
    titlestr = kwargs.get('title', '')
    plottype = kwargs.get('plottype', '')

    if titlestr:
        # increase size of top marign to handle multiline titles
        topmargin = topmargin + 0.15 * len(titlestr.split('\n'))
    thefig = plt.gcf()
    figure_width = float(kwargs.get('width', thefig.get_figwidth())) * w_scale
    thefig.set_figwidth(figure_width)
    if n * plotheight + topmargin > thefig.get_figheight():
        thefig.set_figheight(n * plotheight + topmargin)

    if not data:
        return

    sps = []

    # get the xlabel and xdata
    xlbl = kwargs.get('xlabel', '')
    x_is_time = kwargs.get('x_is_time', xlbl == 'time')
    y_is_time = kwargs.get('y_is_time', False)
    # print 'xlabel = %s' % xlbl
    colors = kwargs.get('color_order', 'rgbcmyk')

    # figure out the right locator and formatter

    def analyzeX(data):
        ''' Return the min and max X values from the set of datasets.'''
        minx = None
        maxx = None
        for var, overlays in data:
            for D in overlays:
                X = D.X
                if X is not None:
                    X = np.array([X])
                    lox, hix = np.nanmin(X), np.nanmax(X)
                    if minx is None:
                        minx, maxx = lox, hix
                    if lox < minx:
                        minx = lox
                    if hix > maxx:
                        maxx = hix

        return minx, maxx

    # get the x range
    if 'xlim' in kwargs:
        (minx, maxx) = kwargs['xlim']
    else:
        minx, maxx = analyzeX(data)
    # get the global y range, if there is one
    miny = None
    maxy = None
    if 'ylim' in kwargs:
        (miny, maxy) = kwargs['ylim']
    if x_is_time:
        # just min and max x if x_is_time
        if minx is not None and minx == maxx:
            maxx = maxx + 1 / 24.

    lastplot = None
    #foo = data[0]
    #data[0] = data[1]
    #data[1] = foo
    dxlblflag = False
    lmaxx = None
    for i, (var, overlays) in enumerate(data):
        empty = True
        if lastplot:
            p = plt.subplot(n, 1, i + 1, sharex=lastplot)
        else:
            p = plt.subplot(n, 1, i + 1)
        lastplot = p
        for j, D in enumerate(overlays):

            # get axis values out of dataset
            X = D.X
            Y = D.Y
            Z = D.Z
            Dkwargs = kwargs.copy()
            # print 'default:' , Dkwargs
            Dkwargs.update(D.kwargs)
            # print 'modified:' , Dkwargs
            Dkwargs['label'] = D.label
            y_is_time = Dkwargs.get('y_is_time', False)
            clabel = Dkwargs.get('zlabel', '')
            if 'ylim' in Dkwargs:
                (lminy, lmaxy) = Dkwargs['ylim']
            else:
                (lminy, lmaxy) = (miny, maxy)
            if 'xlim' in Dkwargs:
                (lminx, lmaxx) = Dkwargs['xlim']
            else:
                (lminx, lmaxx) = (minx, maxx)
            if 'xlabel' in D.kwargs:
                dxlbl = D.kwargs['xlabel']
                dxlblflag = True
            else:
                dxlbl = xlbl
            # strip out Ys that do not correspond to a valid X
            if X is not None:
                ids = np.where(X.mask == False)
                X = X[ids]
                if Y is not None:
                    Y = Y[ids]
                if Z is not None:
                    Z = Z[ids]

            #  check to see if we have at least one tuple
            if Y is not None:
                tuples = Y.count()
            else:
                tuples = 0
            if tuples != 0:
                empty = False

            plottype = Dkwargs.get('plottype', '')
            cticks = Dkwargs.get('cticks', None)
            xformatstring = Dkwargs.get('xformatstring', None)
            yformatstring = Dkwargs.get('yformatstring', None)
            xformatstringm = Dkwargs.get('xformatstringm', None)
            yformatstringm = Dkwargs.get('yformatstringm', None)
            xyearticks = Dkwargs.get('xyearticks', None)
            yyearticks = Dkwargs.get('yyearticks', None)
            yearstr = Dkwargs.get('yearstr', '%Y')
            if plottype == 'arbitrary':
                plotfunc = Dkwargs.get('plotfunc', None)
                plotvars = []
                if X is not None:
                    plotvars.append(X)
                if Y is not None:
                    plotvars.append(Y)
                if Z is not None:
                    plotvars.append(Z)
                plotvars.extend(Dkwargs.get('plotvars', []))
                plotargs = Dkwargs.get('plotargs', {})
                if plotfunc is not None:
                    im = plotfunc(*plotvars, **plotargs)
            elif plottype == 'scattermap':
                scatargs = dict()
                for d in Dkwargs:
                    if d in getvalidscatterargs():
                        scatargs[d] = Dkwargs[d]
                # if data has Z values, then create a scatter plot, with Z as color
                #cmop.debug('starting scatter plot function.')
                im = scattermap(X, Y, color=Z, **scatargs)
                #cmop.debug('finished scatter plot function.')
                #norm = matplotlib.colors.Normalize(vmin = min(Z), vmax = max(Z), clip = False)
                plt.axes(p)
            elif plottype == 'scatter':
                scatargs = dict()
                for d in Dkwargs:
                    if d in getvalidscatterargs():
                        scatargs[d] = Dkwargs[d]
                # if data has Z values, then create a scatter plot, with Z as color
                #cmop.debug('starting scatter plot function.')

                im = scatter3(X, Y, color=Z, **scatargs)

                #cmop.debug('finished scatter plot function.')
                #norm = matplotlib.colors.Normalize(vmin = min(Z), vmax = max(Z), clip = False)
                plt.axes(p)
            elif plottype == 'vlines':
                vlinesargs = dict()
                for d in Dkwargs:
                    if d in getvalidvlinesargs():
                        vlinesargs[d] = Dkwargs[d]
                im = matplotlib.pyplot.vlines(X, Y, Z, **vlinesargs)
            elif plottype == 'fill':
                fillargs = dict()
                for d in Dkwargs:
                    if d in getvalidfillargs():
                        fillargs[d] = Dkwargs[d]
                im = matplotlib.pyplot.fill(X, Y, **fillargs)
            elif plottype == 'pcolor':
                # binsize should come from DB
                binsize = Dkwargs.get('binsize', 45)
                vmin = Dkwargs.get('vmin', -1)
                vmax = Dkwargs.get('vmax', 1)
                k = len(X) / binsize
                x1 = X.reshape((k, binsize))[:, 0]
                y1 = Y.reshape((k, binsize))[0, :]
                (x2, y2) = np.meshgrid(x1, y1)
                z1 = np.reshape(Z, (k, binsize)).transpose()
                cdict = {'red': ((0.0, 0, 1),
                                 (0.01, 0, 0),
                                 (0.35, 0, 0),
                                 (0.66, 1, 1),
                                 (0.89, 1, 1),
                                 (1, 0.5, 0.5)),
                         'green': ((0.0, 0, 1),
                                   (0.01, 0, 0),
                                   (0.125, 0, 0),
                                   (0.375, 1, 1),
                                   (0.640, 1, 1),
                                   (0.910, 0, 0),
                                   (1, 0, 0)),
                         'blue': ((0.0, 0, 1),
                                  (0.01, 0.5, 0.5),
                                  (0.11, 1, 1),
                                  (0.34, 1, 1),
                                  (0.65, 0, 0),
                                  (0.75, 0, 0),
                                  (1, 0, 0))}
                my_cmap = matplotlib.colors.LinearSegmentedColormap(
                    'my_colormap',
                    cdict, 256)
                cdict = cm.jet._segmentdata

                im = plt.pcolor(x2, y2, z1, cmap=my_cmap, vmin=vmin, vmax=vmax)
            else:
                # construct the marker string, using specified color or default
                # color
                if 'color' in Dkwargs:
                    dcolor = Dkwargs['color']
                else:
                    dcolor = colors[(j + 1) % len(colors)]
                marker = Dkwargs.get('style', '.')
                marker = Dkwargs.get('marker', marker)
                if marker == '-':
                    Dkwargs['marker'] = 'None'
                    Dkwargs['linestyle'] = '-'
                else:
                    Dkwargs['marker'] = marker
                    if 'linestyle' not in Dkwargs:
                        Dkwargs['linestyle'] = 'None'
                plotargs = dict()
                for d in Dkwargs:
                    if d in getvalidplotargs():
                        plotargs[d] = Dkwargs[d]
            # print plotargs
                # if data has no ZZ values, use plot
                # if x_is_time & Dkwargs.get('handlegaps', False):
                if Dkwargs.get('handlegaps', False):
                    X, Y = addGaps(X, Y, **Dkwargs)
                #im = plt.gca().plot(X,Y,marker=marker,linestyle=linestyle,color=dcolor,label=lbl)
                im = plt.gca().plot(X, Y, **plotargs)

            # set the x-axis limits
            if lminx is not None:
                plt.gca().set_xlim(lminx, lmaxx)
            else:
                # do nothing
                pass

            # set the y-axis limits
            if lminy is not None:
                plt.gca().set_ylim(lminy, lmaxy)
            else:
                # do nothing
                pass

        # make a FontProperties object
        props = matplotlib.font_manager.FontProperties(size='x-small')

        plottype = kwargs.get('plottype', '')
        col = None
        # make the legend or colorbar
        if plottype == 'scatter' or plottype == 'scattermap':
            colorb_width = 0
            if not empty and not kwargs.get('nocolorbar', False):
                # estimate colorbar height.
                col = colorbarhelper(plt.colorbar(), clabel, cticks=cticks)
                # workaround for matplotlib bug with colorb outside of axes.
                # '.' markers are clipped by the axes bounding box
                # patch submitted to matplotlib-devel for a thorough fix
                if col is not None:
                    # grrr....there seems to be no way to get the final legend width
                    # even forcing a redraw of the whole figure doesn't work
                    # draw()
                    # legend_width = leg.legendPatch.get_width()
                    # instead, estimate width by length of text
                    #tw = max([len(t.get_text()) for t in leg.texts])
                    legend_width = 1.3
                    if legend_width > max_legend_width:
                        max_legend_width = legend_width
                    #thefig.set_figwidth(thefig.get_figwidth() + legend_width)
        elif plottype == 'pcolor':
            colorb_width = 0
            if not empty and not kwargs.get('nocolorbar', False):
                # estimate colorbar height.
                col = colorbarhelper(plt.colorbar(), clabel, cticks=cticks)
                # workaround for matplotlib bug with colorb outside of axes.
                # '.' markers are clipped by the axes bounding box
                # patch submitted to matplotlib-devel for a thorough fix
                if col is not None:
                    # grrr....there seems to be no way to get the final legend width
                    # even forcing a redraw of the whole figure doesn't work
                    # draw()
                    # legend_width = leg.legendPatch.get_width()
                    # instead, estimate width by length of text
                    #tw = max([len(t.get_text()) for t in leg.texts])
                    legend_width = 1.3
                    if legend_width > max_legend_width:
                        max_legend_width = legend_width
                    #thefig.set_figwidth(thefig.get_figwidth() + legend_width)
        else:
            leg = None
            legend_width = 0
            if not empty and not kwargs.get('nolegend', False):
                # estimate legend height.
                legend_height = 0.2 + len(
                    overlays) * 0.35 / thefig.get_figheight()
                if legend_height > plotheight:
                    thefig.set_figheight(topmargin + legend_height * n * 1.5)
                    plotheight = legend_height
                bottom = 0.5 - (legend_height / plotheight) / 2
                leg = plt.legend(
                    prop=props,
                    labelspacing=0.02,
                    borderpad=0.5,
                    loc=(
                        1.01,
                        bottom),
                    numpoints=1)
                plt.axes(p)
                # workaround for matplotlib bug with legend outside of axes.
                # '.' markers are clipped by the axes bounding box
                # patch submitted to matplotlib-devel for a thorough fix
                if leg is not None:
                    for x in leg.legendHandles:
                        try:
                            x._legmarker.set_clip_on(False)
                        except:
                            pass
                    # grrr....there seems to be no way to get the final legend width
                    # even forcing a redraw of the whole figure doesn't work
                    # draw()
                    # legend_width = leg.legendPatch.get_width()
                    # instead, estimate width by length of text
                        tw = 0
                        for t in leg.texts:
                            txt = t.get_text()
                            ttt = txt.split('\n')
                            ttw = max([len(tt) for tt in ttt])
                            if ttw > tw:
                                tw = ttw
                    #tw = max([len(t.get_text()) for t in [t.split('\n') for t in leg.texts]])
                        legend_width = (
                            0.45 * (figure_width / 5.5) + 0.07 * tw)
                    # print "legend width hack: %s" % legend_width
                    # draw_if_interactive()
                    # axes(p)
                    # get_frame().get_window_extent().width does vary with the legend width, but in an extremely opaque way
                    #legend_width = leg.get_frame().get_window_extent().width/thefig.dpi
                    # print "function value: %s " % legend_width
                        if legend_width > max_legend_width:
                            max_legend_width = legend_width
                    #thefig.set_figwidth(thefig.get_figwidth() + legend_width)

        # add this subplot to a list for adjustment later
        sps.append((p, plt.gci(), col))

        if x_is_time:
            # set the formatter and tick locator if the x axis is time data
            if kwargs.get('timetick', '') == '':
                if lminx == lmaxx:
                    lmaxx = lmaxx + 1 / 24.
                loc, fmt = dates.date_ticker_factory(lmaxx - lminx, numticks=6)
                plt.gca().xaxis_date()
                plt.gca().xaxis.set_major_formatter(fmt)
                plt.gca().xaxis.set_major_locator(loc)
            else:
                # determine xaxis range, used to determine tick labels if
                # x-axis is time
                if lmaxx is not None:
                    dtrange = lmaxx - lminx
                    ltime = (lmaxx + lminx) / 2
                    if dtrange < 1827 and xyearticks:
                        dtrange = 1827  # force y-axis tick labels to be at least years
                    pylabutil.adaptive_date_ticks(
                        plt.gca().xaxis, dtrange=dtrange, ltime=ltime,
                        debug=False, formatm=xformatstringm,
                        format=xformatstring, lformat=dxlbl, yearstr=yearstr)

                else:
                    dtrange = None
                    ltime = None

        if y_is_time:
            # set the formatter and tick locator if the x axis is time data
            if kwargs.get('timetick', '') == '':
                if lminy == lmaxy:
                    lmaxy = lmaxy + 1 / 24.
                loc, fmt = dates.date_ticker_factory(lmaxy - lminy, numticks=6)
                plt.gca().yaxis_date()
                plt.gca().yaxis.set_major_formatter(fmt)
                plt.gca().yaxis.set_major_locator(loc)
            else:
                # determine yaxis range, used to determine tick labels if
                # y-axis is time
                if lmaxy is not None:
                    dtrange = lmaxy - lminy
                    if dtrange < 1827 and yyearticks:
                        dtrange = 1827  # force y-axis tick labels to be at least years
                    ltime = (lmaxy + lminy) / 2
                    pylabutil.adaptive_date_ticks(
                        plt.gca().yaxis,
                        dtrange=dtrange,
                        ltime=ltime,
                        debug=False,
                        formatm=yformatstringm,
                        format=yformatstring,
                        lformat=r'%s' %
                        var)
                else:
                    dtrange = None
                    ltime = None
                plt.setp(plt.gca().get_yticklabels(), fontsize=8)
                plt.setp(plt.gca().get_yticklabels(minor=True), fontsize=8)

        # set the ylabel
        plt.ylabel(r'%s' % (var,), size='smaller')
        # set the subplot xlabel if it is specified

        # handle time tick labels using a helper function in pylabutil
        if i + 1 == n or dxlblflag:
            if x_is_time:
                if kwargs.get('timetick', '') == '':
                    plt.setp(plt.gca().get_xticklabels(), 'rotation', 20,
                             'horizontalalignment', 'right', fontsize=8)
                else:
                    # set the formatter and tick locator
                    # adjust the font size
                    plt.setp(plt.gca().get_xticklabels(), fontsize=8)
                    plt.setp(plt.gca().get_xticklabels(minor=True), fontsize=8)
                    # get_xlabel returns a string, rather than a handle object
                plt.xlabel(r'%s' % (plt.gca().get_xlabel()), size='smaller')
            elif dxlblflag:
                plt.xlabel(r'%s' % (dxlbl,), size='smaller')
                plt.setp(plt.gca().get_xticklabels(), fontsize=8)
            else:
                plt.xlabel(r'%s' % (xlbl,), size='smaller')
                plt.setp(plt.gca().get_xticklabels(), fontsize=8)

        else:
            if x_is_time:
                # set the formatter and tick locator
                plt.setp(plt.gca().get_xticklabels(minor=True), visible=False)
                plt.xlabel('')
            plt.setp(plt.gca().get_xticklabels(), visible=False)

        # turn on the grid
        plt.grid(True)

        # shrink the size of the tick labels
        labels = plt.gca().get_xticklabels() + plt.gca().get_yticklabels()
        for label in labels:
            label.set_size(8)

    # set the title
    if titlestr:
        sps[0][0].set_title(titlestr, size='x-small')

    # adjust the spacing and margins for the plots
    figure_width = thefig.get_figwidth()
    if dxlblflag:
        hspace = 0.3
    else:
        hspace = 0.1
    ###############################################
    if 'legendwidth' in kwargs:
        max_legend_width = kwargs['legendwidth']
        while 5.5 * 0.15 / figure_width > 1 - max_legend_width / figure_width:
            max_legend_width -= 1
    plt.subplots_adjust(
        left=5.5 *
        0.15 /
        figure_width,
        right=1 -
        max_legend_width /
        figure_width,
        top=1 -
        topmargin /
        thefig.get_figheight(),
        bottom=0.5 /
        thefig.get_figheight(),
        hspace=hspace)
    #plt.subplots_adjust(left=5.5*0.15/figure_width, right=1-max_legend_width/figure_width, top=1-topmargin/thefig.get_figheight(), bottom=0.5/thefig.get_figheight(), hspace=0.1)

    # Refresh the axes.  Needs to be done after everything else is ready
    for p, imgi, c in sps:
        plt.axes(p)
        if c is not None:
            lims = p.get_position().get_points()

            clims = [
                lims[1][0] + 0.05, lims[0][1],
                1 - (lims[1][0] + 0.05),
                lims[1][1] - lims[0][1]]
            c.ax.set_position(clims)

    # return the figure
    return thefig
