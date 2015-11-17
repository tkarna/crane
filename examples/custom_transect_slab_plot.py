"""
Example script for doing custom plots that combine transects slabs.

Tuomas Karna 2014-09-03
"""
import data.dirTreeManager as dtm
from data.dataContainer import *
from plotting.plotBase import *
from plotting.slabPlot import *
from plotting.transectPlot import *
from plotting.timeSeriesPlot import *
from files.csvStationFile import *
#-------------------------------------------------------------------------
# set metadata and load data sets
#-------------------------------------------------------------------------

# set global font size
fontsize = 16
matplotlib.rcParams['font.size'] = fontsize

# output images
imgDir = createDirectory('example_images')
imgPrefix = 'comboplot_example'
imgFiletype = 'png'

# root for the model data tree
basepath = '/home/workspace/ccalmr53/karnat/runs/newBBL/'
# this rule is used to read files stored in monthly files e.g.
# run007/data/stations/tpoin/tpoin.0.A.elev/201205.nc
rule_monthly = dtm.defaultTreeRule(path=basepath)
# this rule is used to read files stored in single files e.g.
# run007/data/transect/nchannel_salt_0_2012-05-01_2012-05-14.nc
rule = dtm.oldTreeRule(path=basepath)
# metadata for run
runTag = 'run007'

# load transects
# reads file run007/data/transect/nchannel_salt_0_2012-05-01_2012-05-14.nc
transectDC_salt = dtm.getDataContainer(
    rule=rule,
    tag=runTag,
    dataType='transect',
    location='nchannel',
    variable='salt')
transectDC_kine = dtm.getDataContainer(
    rule=rule,
    tag=runTag,
    dataType='transect',
    location='nchannel',
    variable='kine')
# equivalent to (in case you don't want to use dirTreeManager)
#transectDC = dataContainer.loadFromNetCDF(basepath+'/run007/data/transect/nchannel_salt_0_2012-05-01_2012-05-14.nc')

# load slab
# read file run007/data/slab/slab_salt_s1_2012-05-01_2012-05-18.nc
slabDC = dtm.getDataContainer(rule=rule, tag=runTag, dataType='slab',
                              location='slab', variable='salt', slevel=1)

# load time series
elevDC = dtm.getDataContainer(
    rule=rule_monthly,
    tag=runTag,
    dataType='timeseries',
    location='tpoin',
    variable='elev',
    msldepth='0')

# load station file for plotting station locations
stationFile = '/home/users/karnat/share/real_stations.csv'
staFileObj = csvStationFile().readFromFile(stationFile)
stationsToPlot = ['ncbn1', 'saturn01']


def makePlot(timeStamp):
    """
    Function that will make a plot for one time step.

    NOTE timeStamp doesn't need coinside with model exports: the data will be
    interpolated in time if it doesn't coinside with model export.
    """

    #-------------------------------------------------------------------------
    # Create figure and axes
    #-------------------------------------------------------------------------

    width = 12  # inches
    height = 8  # inches
    fig = plt.figure(figsize=(width, height))

    # We'll use gridspec to create axes in rectangular 6-by-5 lattice
    import matplotlib.gridspec as gridspec
    nrows = 6
    ncols = 5
    Grid = gridspec.GridSpec(nrows, ncols)

    # axis for elevation time series
    axElev = fig.add_subplot(Grid[:2, :2])  # first 2 rows, first 2 columns
    # axis for slab
    axSlab = fig.add_subplot(Grid[:2, 2:])  # first 2 rows, columns > 2
    # and the transects
    axTran1 = fig.add_subplot(Grid[2:4, :])  # rows 2,3,4, all columns
    # rows 5,6,7, all columns, share x/y axis with previous (sets same ticks
    # etc)
    axTran2 = fig.add_subplot(Grid[4:6, :], sharex=axTran1, sharey=axTran1)

    # gridspec allows to tune the spacing between plots (unit is fraction of
    # font size)
    boundary_pad = 3.5
    horizontal_pad = 0.2
    vertical_pad = 1.0
    # figure area left,bottom,right,top in normalized coordinates [0,1]
    bounds = [0, 0, 1, 1]
    Grid.tight_layout(
        fig,
        pad=boundary_pad,
        w_pad=horizontal_pad,
        h_pad=vertical_pad,
        rect=bounds)

    #-------------------------------------------------------------------------
    # Create plots
    #-------------------------------------------------------------------------

    # for all avaiable colormaps see ( '_r' reverses the colormap )
    # http://matplotlib.org/examples/color/colormaps_reference.html
    colormap = plt.get_cmap('Spectral_r')
    colormap_kine = plt.get_cmap('gist_heat')

    # slab
    salt_clim = [0, 32]
    ncontours = 16
    # bouding box for slab [xmin,xmax,ymin,ymax] in model x,y coordinates
    estuarybbox = [330000, 360000, 284500, 297500]
    dia = slabSnapshotDC(
        clabel='Salinity',
        unit='psu',
        clim=salt_clim,
        cmap=colormap)
    dia.setAxes(axSlab)
    dia.addSample(slabDC, timeStamp=timeStamp, plotType='contourf',
                  bbox=estuarybbox, N=ncontours)
    # overrides default format for colorbar floats
    dia.showColorBar(format='%.2g')
    #dia.addTitle('in case you want a custom title')
    # get transect (x,y) coordinates from the transectDC
    transectXYCoords = generateTransectFromDataContainer(transectDC_salt, 0)[4]
    # plot transect on the map (thin black on thick white)
    dia.addTransectMarker(transectXYCoords[:, 0], transectXYCoords[:, 1],
                          color='w', linewidth=2.0)
    dia.addTransectMarker(transectXYCoords[:, 0], transectXYCoords[:, 1],
                          color='k', linewidth=1.0)
    # plot station markers
    for station in stationsToPlot:
        staX = staFileObj.getX(station)
        staY = staFileObj.getY(station)
        dia.addStationMarker(
            staX,
            staY,
            label=station,
            printLabel=True,
            marker='*')
    # add text to plot. x,y are in normalized axis coordinates [0,1]
    dia.ax.text(0.05, 0.98, 'custom text', fontsize=fontsize,
                verticalalignment='top', horizontalalignment='left',
                transform=dia.ax.transAxes)

    # elevation time series
    # define the time range to plot
    elevStartTime = datetime.datetime(2012, 5, 4, 0, 0)
    elevEndTime = datetime.datetime(2012, 5, 5, 0, 15)
    elevMeanTime = elevStartTime + (elevEndTime - elevStartTime) / 2
    elevLim = [-1.5, 2.5]
    dia = timeSeriesPlotDC2(
        xlabel=elevMeanTime.strftime('%Y %b %d'),
        ylim=elevLim)
    dia.setAxes(axElev)
    #dia.addShadedRange( timeStamp, timeStamp+datetime.timedelta(seconds=30), facecolor='IndianRed')
    dia.addShadedRange(
        timeStamp,
        timeStamp,
        edgecolor='IndianRed',
        facecolor='none',
        linewidth=2)
    tag = elevDC.getMetaData('tag')
    dia.addSample(
        elevDC.timeWindow(
            elevStartTime,
            elevEndTime),
        label=tag,
        color='k')
    dia.addTitle('Elevation ({0:s}) [m]'.format(
        elevDC.getMetaData('location').upper()))
    # adjust the number of ticks in x/y axis
    dia.updateXAxis(maxticks=5)
    dia.updateYAxis(maxticks=3, prune='lower')

    # transects
    dia = transectSnapshotDC(
        clabel='Salinity',
        unit='psu',
        cmap=colormap,
        clim=salt_clim)
    dia.setAxes(axTran1)
    #transectDC_salt.data *= 1e-3
    dia.addSample(transectDC_salt, timeStamp, N=ncontours)
    dia.addTitle('')
    dia.showColorBar()
    # plot station markers
    for station in stationsToPlot:
        staX = staFileObj.getX(station)
        staY = staFileObj.getY(station)
        dia.addStationMarker(staX, staY, label=station, color='k',
                             linewidth=1.5, linestyle='dashed')
    # do not show x axis ticks and label for this plot
    dia.hideXTicks()

    dia = transectSnapshotDC(clabel='TKE', unit='m2s-1', logScale=True,
                             clim=[-7, -2], climIsLog=True, cmap=colormap_kine)
    dia.setAxes(axTran2)
    dia.addSample(transectDC_kine, timeStamp, N=ncontours)
    # plot station markers
    for station in stationsToPlot:
        staX = staFileObj.getX(station)
        staY = staFileObj.getY(station)
        dia.addStationMarker(staX, staY, label=station, color='k',
                             linewidth=1.5, linestyle='dashed')
    dia.addTitle('')
    dia.showColorBar()
    dia.updateXAxis(maxticks=15)
    dia.updateYAxis(maxticks=6)

    #-------------------------------------------------------------------------
    # Save to disk
    #-------------------------------------------------------------------------
    dateStr = timeStamp.strftime('%Y-%m-%d_%H-%M')
    filename = '_'.join([imgPrefix, dateStr])
    saveFigure(
        imgDir,
        filename,
        imgFiletype,
        verbose=True,
        dpi=200,
        bbox_tight=True)
    plt.close()

#-------------------------------------------------------------------------
# Loop over a range of time instances
#-------------------------------------------------------------------------

# loop over a time range with 15 min inrements
loopStart = datetime.datetime(2012, 5, 4)
loopEnd = datetime.datetime(2012, 5, 5)
from dateutil import relativedelta
loopStep = relativedelta.relativedelta(minutes=15)

# plot a specific time stamp
timeStamp = loopStart
while timeStamp <= loopEnd:
    makePlot(timeStamp)
    timeStamp += loopStep
