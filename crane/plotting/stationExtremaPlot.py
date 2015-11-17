#!/usr/bin/python

"""
A class for plotting min/max values along the estuary.

Tuomas Karna 2012-09-25
"""
import numpy as np
# TODO import only modules
from crane.data.dataContainer import dataContainer
from crane.plotting import plotBase
from crane import plt

# class that takes a list of station names and distance coordinate, and time series
# TODO inherit from plotBase?


class stationExtremaPlot(object):

    def __init__(self, fieldName, stationCoords, unit='', **defArgs):
        """Creates a new object for given stations and (along river) coordinates

        Args:
        stationNames  --  (list of string) station names
        stationCoords --  (list of float) station along river coordinates
        """
        self.fieldName = fieldName
        self.stationNames = stationCoords.keys()
        self.stationCoords = stationCoords
        self.unit = unit
        self.xRelOffset = 0.05
        # default color sequence for samples (typ 'r' for observation)
        self.defColors = ['r', 'b', 'g', 'k', 'm', 'c']
        self.ax = None
        self.figsize = (10, 5)
        self.data = []
        self.labelStrings = []
        self.labelSymbols = []
        self.defaultArgs = {}
        self.defaultArgs['linestyle'] = 'none'
        self.defaultArgs.update(defArgs)
        self.defaultMeanArgs = dict(self.defaultArgs)
        self.defaultMeanArgs['markersize'] = 6
        #self.defaultMeanArgs['markersize'] = 8
        self.defaultMeanArgs['markeredgewidth'] = 2.5
        self.defaultMeanArgs['markerfacecolor'] = 'none'
        self.defaultMeanArgs['marker'] = 'o'
        self.defaultMinMaxArgs = dict(self.defaultArgs)
        self.defaultMinMaxArgs['markersize'] = 12
        self.defaultMinMaxArgs['markeredgewidth'] = 2.5
        self.defaultMinMaxArgs['markerfacecolor'] = 'none'
        self.defaultMinMaxArgs['marker'] = '_'
        self.minSeparation = 0.04  # percentage of original x-axis length
        self.fixCoordinateSpacing()

    def createDefaultFigure(self):
        """Creates default figure, create default axes and assing them.
        Will be called automatically if user defined axes are not present."""
        fig = plt.figure(figsize=self.figsize)
        self.setupAxes(fig)

    def fixCoordinateSpacing(self):
        tickLabels = np.array(self.stationNames)
        x = np.array(self.stationCoords.values(), dtype=float)
        sortedIx = x.argsort()
        x = x[sortedIx]
        tickLabels = tickLabels[sortedIx]
        dx = np.diff(x)
        dxThreshold = self.minSeparation * (x[-1] - x[0])
        dx[dx < dxThreshold] = dxThreshold
        x = np.hstack(([0], np.cumsum(dx))) + x[0]
        self.stationCoords = dict(zip(tickLabels, x))

    def createAxes(self, fig, rect=111):
        """Create rectangular axes for the diagram. Returns the axes.
        """
        ax = fig.add_subplot(rect)
        fig.subplots_adjust(top=0.93, bottom=0.14, left=0.08, right=0.98)
        ax.grid()
        varStr = self.fieldName
        if self.unit:
            varStr += ' [' + self.unit + ']'
        # ax.set_xlabel(errStr)
        ax.set_ylabel(varStr)
        tickLabels = np.array([sta.replace('saturn', 'sat')
                               for sta in self.stationCoords.keys()])
        x = self.stationCoords.values()
        ax.set_xticks(x)
        ax.set_xticklabels(tickLabels, rotation=90)
        xmin = min(self.stationCoords.values())
        xmax = max(self.stationCoords.values())
        xRange = xmax - xmin
        offset = 0.05
        xlim = [xmin - offset * xRange, xmax + offset * xRange]
        ax.set_xlim(xlim)

        return ax

    def setAxes(self, ax, parent_ax=None):
        """Set axes for the diagram. All data will be plotted in these axes.
        """
        self.ax = ax

    def setupAxes(self, fig, rect=111):
        """Creates new axes and set them.
        A shorthand for calling createAxes and setAxes.
        Returns the axes."""
        ax = self.createAxes(fig, rect)
        self.setAxes(ax)
        return ax

    def plotSampleSet(self, label, stationDict, **kwargs):
        """Plots a set of station in one command

        Args:
          label -- (string) label to identify the dataset (e.g. 'obs', 'db22')
          stationDict -- (dict) {'stationName':dataArray}, where dataArray is singleton numpy array
          kwargs -- key word argument that will be passed to ax.plot command
        """
        for sta in stationDict:
            if not 'color' in kwargs:
                kwargs['color'] = self.defColors[len(self.labelStrings)]
            self.plotSample(label, sta, stationDict[sta], **kwargs)

    def plotSample(self, label, station, sample, **kwargs):
        """Plots one sample in the diagram.

        Args:
          label -- (string) label to identify the dataset (e.g. 'obs', 'db22')
          stationName -- (sting) must be one of the names in self.stationNames, otherwise ignored
          sample -- (np.array) time series data to plot
          kwargs -- key word argument that will be passed to ax.plot command
        """
        if not self.ax:
            self.createDefaultFigure()
        sMin = sample.min()
        sMax = sample.max()
        sMean = sample.mean()
        if station not in self.stationNames:
            print 'station ', station, 'not defined'
            return
        minMaxArgs = dict(self.defaultMinMaxArgs)
        minMaxArgs.update(kwargs)
        meanArgs = dict(self.defaultMeanArgs)
        meanArgs.update(kwargs)
        meanArgs['markeredgecolor'] = meanArgs['color']
        self.ax.plot(self.stationCoords[station], sMin, **minMaxArgs)
        p = self.ax.plot(self.stationCoords[station], sMean, **meanArgs)
        self.ax.plot(self.stationCoords[station], sMax, **minMaxArgs)
        if label not in self.labelStrings:
            self.labelStrings.append(label)
            self.labelSymbols.append(p[0])

    def isEmpty(self):
        return len(self.labelStrings) == 0

    def addTitle(self, titleStr, **kwargs):
        """Adds a title in the axes. Optional arguments are passed to set_title routine"""
        if self.ax:
            self.ax.set_title(titleStr, **kwargs)

    def showLegend(self, **kwargs):
        """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
        if self.ax:
            self.ax.legend(
                self.labelSymbols,
                self.labelStrings,
                numpoints=1,
                **kwargs)


class stationExtremaPlotDC(stationExtremaPlot):
    """stationExtremaPlot class that uses dataContainers"""

    def plotSample(self, label, station, dataCont, **kwargs):
        """Plots one sample in the digram.

        Args:
          label -- (string) label to identify the dataset (e.g. 'obs', 'db22')
          stationName -- (sting) must be one of the names in self.stationNames, otherwise ignored
          sample -- (dataContainer) time series data to plot
          kwargs -- key word argument that will be passed to ax.plot command
        """
        sample = np.squeeze(dataCont.data)
        stationExtremaPlot.plotSample(self, label, station, sample, **kwargs)


if __name__ == '__main__':

    from crane.data import timeArray
    from datetime import datetime

    # examples with numpy array inputs
    # generate data
    startTime = datetime(2010, 1, 12, 0, 0, 0)
    endTime = datetime(2010, 2, 13, 3, 30, 0)
    dt = 900.0
    ta = timeArray.generateSimulationTimeArray(
        startTime, endTime, dt).asEpoch()
    t = ta.array

    T = 44714
    ref = np.sin(t / T) + 0.95 * np.sin(0.95 * t / T)
    m1 = 0.90 * np.sin(t / T) + 0.7 * np.sin(0.85 * t / T + 0.3) - 0.08
    m2 = 0.80 * np.sin(t / T + 0.8) + 0.9 * np.sin(0.95 * t / T + 0.5) - 0.12
    m3 = 0.78 * np.sin(t / T) + 0.82 * np.sin(0.90 * t / T) + 0.03

    stationX = [10, 20, 21, 80, 90]
    stationNames = ['mouth', 'sta01', 'sta02', 'sta03', 'river']
    stationCoords = dict(zip(stationNames, stationX))

    dia = stationExtremaPlot('Elevation', stationCoords, unit='m')
    dia.plotSample('obs', 'mouth', 1.0 * ref, color='r')
    dia.plotSample('obs', 'sta01', 0.9 * ref, color='r')
    dia.plotSample('obs', 'sta02', 0.8 * ref, color='r')
    dia.plotSample('obs', 'sta03', 0.5 * ref, color='r')
    dia.plotSample('obs', 'river', 0.4 * ref, color='r')
    dia.plotSample('model 1', 'mouth', 1.0 * m1, color='b')
    dia.plotSample('model 1', 'river', 0.35 * m1 + 0.04, color='b')
    dia.addTitle('stationExtremaPlot example')
    dia.showLegend()
    plt.show()

    dia = stationExtremaPlot(
        'Elevation',
        stationCoords,
        unit='m',
        markersize=20)
    dataset = dict()
    dataset['mouth'] = 1.0 * ref
    dataset['sta01'] = 0.9 * ref
    dataset['sta02'] = 0.8 * ref
    dataset['sta03'] = 0.5 * ref
    dataset['river'] = 0.4 * ref
    dia.plotSampleSet('obs', dataset)
    dataset = dict()
    dataset['mouth'] = 1.0 * m3
    dataset['sta01'] = 0.9 * m3
    dataset['sta02'] = 0.8 * m3
    dataset['sta03'] = 0.5 * m3
    dataset['river'] = 0.4 * m3
    dia.plotSampleSet('model 3', dataset)
    dia.addTitle('stationExtremaPlot example')
    dia.showLegend()
    plt.show()

    # example with dataContainer objects
    # generate data
    d0 = dataContainer.fromTimeSeries('Observation', ta, ref, ['elev'])
    d1 = dataContainer.fromTimeSeries('model Eins', ta, m1, ['elev'])
    d2 = dataContainer.fromTimeSeries('model Zwei', ta, m2, ['elev'])
    d3 = dataContainer.fromTimeSeries('model Drei', ta, m3, ['elev'])

    dia = stationExtremaPlotDC('Elevation', stationCoords, unit='m')
    dataset = dict()
    for i, station in enumerate(dia.stationNames):
        dd = d0.copy()
        dd.data *= 1 - 0.1 * i
        dataset[station] = dd
    dia.plotSampleSet('obs', dataset)
    for i, station in enumerate(dia.stationNames):
        dd = d1.copy()
        dd.data *= 1 - 0.1 * i
        dataset[station] = dd
    dia.plotSampleSet('model 1', dataset)
    dia.addTitle('stationExtremaPlotDC example')
    dia.showLegend()
    plt.show()
