"""
Implements quick plot function that plots any dataContainer using the available metadata.

Tuomas Karna 2016-02-09
"""
from crane import plt, matplotlib
from crane.plotting.timeSeriesPlot import timeSeriesPlotDC2
from crane.physicalVariableDefs import UNITS, VARS
import crane.data.dataContainer as dataContainer

# TODO add support for all plot types
# TODO must be able to plot multiple dcs at once

TIMESERIES_TYPES = ['timeseries', 'sil', 'volintegral']

def plotDataContainer(dclist0, ax=None, show=False, dia=None, **kwargs):
    """
    Makes a default plot of the given dataContainer.
    """
    if isinstance(dclist0, dataContainer.dataContainer):
        dclist = [dclist0]
    else:
        dclist = dclist0

    # create axes if necessary
    create_new_axes = ax is None
    if create_new_axes:
        # create new axis
        fig = plt.figure(figsize=(12, 6))

    # get all datatypes
    datatypelist = []
    for dc in dclist:
        datatype = dc.getMetaData('dataType')
        if datatype not in datatypelist:
            datatypelist.append(datatype)

    # create axes
    if create_new_axes:
        # default: create new axes for each dataContainer
        nb_axes = len(dclist)
        axes_dict = {}
        for i, dc in enumerate(dclist):
            axes_dict[dc] = fig.add_subplot(nb_axes, 1, i+1)
        
        def get_axes(dc):
            return axes_dict[dc]

    for datatype in datatypelist:

        dc = dclist[0]
        # Figure out correct metadata
        meta = dc.getMetaData()
        datatype = meta.get('dataType')
        if datatype is None:
            raise Exception('dataContainer does not have dataType attribute. I do not know how to continue')
        var = meta.get('variable', 'variable')
        varstr = VARS.get(var, var)
        unit = UNITS.get(var, '-')
        tag = meta.get('tag', 'tag')
        ylabel = kwargs.pop('ylabel', varstr)
        unit = kwargs.pop('unit', unit)
        label = kwargs.pop('label', None)
        showLegend = kwargs.pop('showLegend', True)

        # create a diagram
        if datatype in TIMESERIES_TYPES:
            if create_new_axes:
                # default: each dc in its own subplot
                for dc in dclist:
                    ax = get_axes(dc)
                    tag = dc.getMetaData('tag')
                    var = dc.getMetaData('variable')
                    varstr = VARS.get(var, var)
                    unit = UNITS.get(var, '-')
                    if label is None:
                        label = tag
                    ylabel = varstr
                    if dia is None:
                        dia = timeSeriesPlotDC2(ylabel=ylabel, unit=unit)
                        dia.setAxes(ax)
                    dia.addSample(dc, label=label, **kwargs)
                    dia.addTitle(tag)
            else:
                # default: plot all to same axes and diagram object
                if dia is None:
                    dia = timeSeriesPlotDC2(ylabel=ylabel, unit=unit)
                    dia.setAxes(ax)
                for dc in dclist:
                    tag = dc.getMetaData('tag')
                    if label is None:
                        label = tag
                    dia.addSample(dc, label=label, **kwargs)
                if showLegend:
                    dia.showLegend()
        else:
            raise NotImplementedError('plotting dataType "{:}" is not currenlyt supported'.format(datatype))

    if show:
        plt.show()

    return dia
