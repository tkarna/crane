#!/usr/bin/python

"""
A class for plotting tidal constituents.
Uses matplotlib bar plotting routine.


Tuomas Karna 2012-09-17
"""
import numpy as np
from plotting.plotBase import *

# class that takes the data (error time series) and plots the histogram
class barPlotBase(plotBase) :
  """A base class for doing bar plots"""
  def __init__(self, **defArgs) :
    # default plot options for all diagrams
    self.labels = []
    self.defaultArgs = {}
    self.defaultArgs['unit'] = ''
    self.defaultArgs['bins'] = 30
    self.defaultArgs['normed'] = False
    self.defaultArgs['barGap'] = 0.2
    self.defaultArgs['xlabel'] = ''
    self.defaultArgs['ylabel'] = ''
    # override with user input
    self.defaultArgs.update( defArgs )
    # remove unit from the arg list (doesn't work with plotting commands)
    self.xlabel = self.defaultArgs.pop('xlabel')
    self.ylabel = self.defaultArgs.pop('ylabel')
    self.unit = self.defaultArgs.pop('unit')
    self.barGap = self.defaultArgs.pop('barGap')
    self.ax = None
    # set of all constituents
    self.constituents = set()
    self.data = []
    self.nDataSets = 0
    self.rects = []
    self.sampleArgs = []
    # default color sequence for samples (typ 'r' for observation)
    self.colors = ['r', 'b','g','k','m','c']

  def setAxes(self, ax) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    self.ax = ax
    ax.grid()
    yStr = self.ylabel
    if self.unit :
      yStr += ' [' + self.unit + ']'
    ax.set_ylabel(yStr, multialignment='center')

  def addSample(self, constDict, label,**kwargs) :
    """Adds constituents to the diagram, identified with the given label string."""
    constituents = constDict.keys()
    self.constituents.update( constituents )
    self.data.append( constDict )
    self.labels.append( label )
    self.sampleArgs.append( kwargs )
    
  def makePlot(self, include=[], exclude=[] ) :
    """Plots all the constituents in the data sets.
    Use include (list of const. names) to include only certain constituents (if present).
    Use exclude (list of const. names) to exclude certain constituents. 
    """
    N = len(self.data)
    width = 1.0-self.barGap
    constituents = list(self.constituents)
    constituents.sort()
    if include :
      constituents = include
    if exclude :
      constituents = constituents.difference(set(exclude))
    constituents = list(constituents)
    # x coord of left corner of each box
    xLocBase = np.arange(len(constituents))
    xLocDict = dict(zip(constituents, xLocBase))
    for i,d in enumerate(self.data) :
      restrictedConsts = set(constituents).intersection(set(d.keys()))
      xLoc = np.array( [ xLocDict[const] for const in restrictedConsts ], dtype=float )
      vals = np.array( [ d[const] for const in restrictedConsts ], dtype=float )
      # x coord of left corner for i-th data set
      xLoc += self.barGap/2 + i*width/N
      sa = dict(self.sampleArgs[i])
      if not 'color' in sa :
        # use standard sequence
        sa['color'] = self.colors[i]
      rects = self.ax.bar(xLoc, vals, width/N, **sa)
    self.ax.set_xticks(xLocBase+0.5)
    self.ax.set_xticklabels( constituents, rotation=90 )
    # add black zero line
    self.ax.plot( [self.ax.dataLim.xmin, self.ax.dataLim.xmax], [0,0], 'k-', label='_nolegend_' )

  def addTitle(self, titleStr, **kwargs) :
    """Adds a title in the axes. Optional arguments are passed to set_title routine"""
    self.ax.set_title(titleStr, **kwargs)

  def showLegend(self, *args, **kwargs) :
    """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
    self.ax.legend(self.labels, numpoints=1, *args, **kwargs)

class amplitudePlot(barPlotBase) :
  """Amplitude bar plot"""
  def __init__(self, **defArgs) :
    if not 'ylabel' in defArgs.keys() :
      defArgs['ylabel'] = 'Amplitude'
    barPlotBase.__init__(self,**defArgs)
  
class phasePlot(barPlotBase) :
  """Phase bar plot"""
  def __init__(self, **defArgs) :
    if not 'ylabel' in defArgs.keys() :
      defArgs['ylabel'] = 'Phase'
    barPlotBase.__init__(self,**defArgs)
  def makePlot(self, include=[], exclude=[]) :
    barPlotBase.makePlot(self,include,exclude)
    self.ax.set_yticks( np.arange(0,380,60) )
  
class amplitudePhasePlot(object) :
  """Combined phase and amplitude plot """
  def __init__(self, **defArgs) :
    self.fig = plt.figure(figsize=(9,5))
    self.amp = amplitudePlot(**defArgs)
    defArgs['unit'] = 'deg'
    self.pha = phasePlot(**defArgs)
    
    self.amp.setupAxes(self.fig, 211)
    self.pha.setupAxes(self.fig, 212)
  
  def addSample(self, tidalConsts, label, **kwargs) :
    """Adds tidal constituents to the plot.
    tidalConsts is a tidalConstituents object.
    """
    self.amp.addSample( tidalConsts.amplitude, label, **kwargs )
    self.pha.addSample( tidalConsts.phase, label, **kwargs )
    
  def makePlot(self, include=[], exclude=[], ampThreshold=0.0, nConsts=0) :
    """Plots all the constituents in the data sets.
    Use include (list of const. names) to include only certain constituents (if present).
    Use exclude (list of const. names) to exclude certain constituents.
    Use ampThreshold to include constituents whose amplitude is above the given threshold.
    Use nConsts to include only nConst larges constituents (in amplitude).
    """
    constituents = filterConstituents( self.amp, include, exclude, ampThreshold, nConsts )
    self.amp.makePlot(include=constituents)
    self.pha.makePlot(include=constituents)
    # no xlabel for amp plot
    self.amp.ax.set_xticklabels([])
    self.amp.ax.set_xlabel('')
    
  def showLegend(self, *args, **kwargs) :
    """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
    if not 'loc' in kwargs.keys() :
      kwargs['loc'] = [0.93,0.3]
    self.amp.showLegend(*args, **kwargs)
    
  def addTitle(self, titleStr, **kwargs) :
    self.fig.suptitle(titleStr, **kwargs)

def filterConstituents( ampPlot, include=[], exclude=[], ampThreshold=0.0, nConsts=0 ) :
    constituents = ampPlot.constituents
    if include :
      constituents = constituents.intersection(set(include))
    if exclude :
      constituents = constituents.difference(set(exclude))
    constituents = list(constituents)
    if ampThreshold or nConsts :
      # sort magnitudes
      maxAmp = np.zeros(len(constituents))
      for i,d in enumerate(ampPlot.data) :
        for j,const in enumerate(constituents) :
          maxAmp[j] = max( maxAmp[j], d[const] )
      # find good ones
      goodIx = []
      if ampThreshold :
        goodIx = np.nonzero(maxAmp >= ampThreshold)[0]
      if nConsts and nConsts < len(goodIx) :
        orderedIx = np.argsort(maxAmp)
        goodIx = orderedIx[-nConsts:]
      # discard bad ones
      constituents = [ constituents[i] for i in goodIx ]
    if not constituents :
      print 'warning : no constituents match the criteria, plotting all'
    return constituents
  
    
class stackAmplitudePhasePlot(object) :
  def __init__(self, **defArgs) :
    self.defArgs = defArgs
    figsize = defArgs.pop( 'figsize' , (10,10) )
    self.fig = plt.figure( figsize=figsize )
    self.amplim = defArgs.pop( 'amplim', [] )
    self.phalim = defArgs.pop( 'phalim', [0,360] )
    self.plots = dict() # stores all diagrams
    self.tags = [] # to keep ordering
    self.legends = dict() # keep track of all legends
    self.axGrid = None

  def hasPlot( self, tag ) :
    return tag in self.tags

  def addPlot(self, tag, **kwargs) :
    """Adds a new subplot to the diagram"""
    self.tags.append( tag )
    kw = dict(self.defArgs)
    kw.update(kwargs)
    self.unit = kw.pop('unit')
    self.plots[tag] = {}
    self.plots[tag]['amp'] = amplitudePlot(**kw)
    kw['ylabel'] = ''
    self.plots[tag]['pha'] = phasePlot(**kw)

  def addSample( self, tag, sample, label, **kwargs ) :
    if not tag in self.tags :
      self.addPlot( tag )
    self.plots[tag]['amp'].addSample( sample.amplitude, label, **kwargs )
    self.plots[tag]['pha'].addSample( sample.phase, label, **kwargs )
    
  def makePlot(self, **kwargs) :
    constituents = set()
    include = kwargs.pop('include',[])
    exclude = kwargs.pop('exclude',[])
    ampThreshold = kwargs.pop('ampThreshold',0.0)
    nConsts = kwargs.pop('nConsts',0)
    nPlots = len(self.tags)
    if nPlots == 0 :
      return # no data
    self.axGrid = []
    for i,tag in enumerate(self.tags) :
      if i == 0 :
        axes = [ self.fig.add_subplot(nPlots,2,2*i+1),
                self.fig.add_subplot(nPlots,2,2*i+2) ]
      else :
        amp0 = self.axGrid[0][0]
        pha0 = self.axGrid[0][1]
        axes = [ self.fig.add_subplot(nPlots,2,2*i+1,sharey=amp0),
                 self.fig.add_subplot(nPlots,2,2*i+2,sharey=pha0) ]
      self.axGrid.append( axes )

      self.plots[tag]['amp'].setupAxes(self.fig, self.axGrid[i][0])
      self.plots[tag]['pha'].setupAxes(self.fig, self.axGrid[i][1])
      
      const = filterConstituents( self.plots[tag]['amp'], include, exclude, ampThreshold, nConsts )
      constituents = constituents.union( set(const) )
    constituents = list( constituents )
    for i,tag in enumerate(self.tags) :
      self.plots[tag]['amp'].makePlot(include=constituents)
      self.plots[tag]['pha'].makePlot(include=constituents)
    for i in range(len(self.axGrid)-1) :
      self.axGrid[i][0].set_xticklabels('')
      self.axGrid[i][1].set_xticklabels('')
    for i in range(len(self.axGrid)) :
      if self.phalim :
        self.axGrid[i][1].set_ylim(self.phalim)
        self.axGrid[i][1].set_yticks( np.arange(self.phalim[0],self.phalim[1]+10,60) )
      if self.amplim :
        self.axGrid[i][0].set_ylim(self.amplim)
    # set titles
    titleStr = 'Amplitude'
    if self.unit :
      titleStr += ' ['+self.unit+']'
    self.axGrid[0][0].set_title( titleStr )
    titleStr = 'Phase'
    unit = 'deg'
    titleStr += ' ['+unit+']'
    self.axGrid[0][1].set_title( titleStr )

  def isEmpty( self ) :
    for tag in self.plots :
      for d in self.plots[tag]['amp'].data :
        return False
    return True

  def addTitle(self, titleStr, **kwargs) :
    """Adds a title in the axes. Optional arguments are passed to set_title routine"""
    self.fig.suptitle(titleStr, **kwargs)

  def showLegend(self, **kwargs) :
    """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
    if not 'loc' in kwargs :
      kwargs['loc'] = [0.7,1.0]
    #for tag in self.tags :
    for tag in [self.tags[0]] :
      self.plots[tag]['pha'].showLegend(prop=dict(size='small'),**kwargs)

def test() :
  labels = ['O1', 'K1', 'K2', 'L2', 'N2', 'M2', 'M4', 'M6']
  amp    = [ 1.2, 0.23,  0.1, 0.13, .341, .512, .433,  .31]
  pha    = [ 212,   45,  124,   54,  192,   76,   85,   39]
  
  ampDict = dict(zip(labels, amp))
  phaDict = dict(zip(labels, pha))

  from data.harmonicAnalysis import tidalConstituents
  tc = tidalConstituents('description', ampDict, phaDict, 'elev', startTime=None, endTime=None )
  
  fig = plt.figure()
  dia = amplitudePlot(unit='m')
  ax = dia.setupAxes(fig)
  dia.addSample(ampDict,'one')
  dia.addSample(ampDict,'two')
  dia.addSample(ampDict,'Drei')
  dia.addSample(ampDict,'fyra')
  dia.makePlot(exclude=['O1'])
  dia.showLegend()
  dia.addTitle('amplitudePlot example')
  #plt.savefig('foo1.png')
  plt.show()

  fig = plt.figure()
  dia = phasePlot(unit='deg')
  ax = dia.setupAxes(fig)
  dia.addSample(phaDict,'one',color='b')
  dia.addSample(phaDict,'two',color='r')
  dia.addSample(phaDict,'Drei',color='g')
  dia.addSample(phaDict,'fyra',color='k')
  dia.makePlot(include=['L2', 'N2', 'M2', 'M4', 'M6'])
  dia.showLegend()
  dia.addTitle('phasePlot example')
  #plt.savefig('foo2.png')
  plt.show()
  
  dia = amplitudePhasePlot(unit='m')
  dia.addSample(tc,'one',color='b')
  dia.addSample(tc,'two',color='r')
  dia.addSample(tc,'Drei',color='g')
  dia.addSample(tc,'fyra',color='k')
  dia.makePlot(ampThreshold=0.11,nConsts=5)
  dia.showLegend()
  dia.addTitle('amplitudePhasePlot example')
  #plt.savefig('foo3.png')
  plt.show()

  dia = stackAmplitudePhasePlot(unit='m',ylim=[0,0.5])
  dia.addPlot('top',ylabel='ssss')
  dia.addSample('top',tc,'one',color='b')
  dia.addSample('mid',tc,'two',color='r')
  dia.addSample('bot',tc,'Drei',color='g')
  dia.addSample('top',tc,'fyra',color='k')
  dia.makePlot(ampThreshold=0.11)
  dia.showLegend()
  dia.addTitle('stackAmplitudePhasePlot example')
  #plt.savefig('foo3.png')
  plt.show()

  # examples with real data
  from data.dataContainer import dataContainer
  
  dataDir = '/home/karnat/temp/db29skill/data/obs/'
  d1 = dataContainer.loadFromNetCDF(dataDir+'tpoin_elev_0_2011-05-13_2011-07-19.nc')
  dataDir = '/home/karnat/temp/db29skill/data/run15_29/'
  d2 = dataContainer.loadFromNetCDF(dataDir+'tpoin_elev_0_2011-05-14_2011-07-19.nc')
  dia = amplitudePhasePlot(unit='m')
  tc1 = tidalConstituents.computeFromData( d1 )
  dia.addSample(tc1,'obs')
  tc2 = tidalConstituents.computeFromData( d2 )
  dia.addSample(tc2,'run15_29')
  dia.makePlot()
  dia.showLegend()
  dia.addTitle('elevation tpoin')
  plt.savefig('foo4.png')
  #plt.show()

if __name__=='__main__':
  test()

