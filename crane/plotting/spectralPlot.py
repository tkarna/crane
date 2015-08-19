#!/usr/bin/python

"""
A class for plotting signals in frequency domain.


Tuomas Karna 2012-10-04 17:17:23
"""
import numpy as np
# TODO import only modules
from crane.data.dataContainer import dataContainer
from crane.data.periodogram import *
from crane.plotting.plotBase import *

# TODO create a comprehensive class or use TAPPY
class tidalConstituents(object) :
  # angular speeds (deg/hour)
  T = 15.00 # the rotation of the Earth on its axis, with respect to the Sun
  h = .04106864 # the rotation of the Earth about the sun
  s = .54901653 #the rotation of the Moon about the Earth
  p = .00464183 #the precession of the Moon's perigee
  angularSpeeds = {
  'M2':     2*T-2*s+2*h,
  'N2':   2*T-3*s+2*h+p,
  'S2':           2*T,
  'K1':          T+h,
  'L2':    2*T-s+2*h-p,
  'O1':       T-2*s+h,
  'Sa':            h,
  'nu2':  2*T-3*s+4*h-p,
  'K2':        2*T+2*h,
  'Mm':          s-p,
  'P1':          T-h,
  'M4': 57.968208468,
  'M6': 86.952312720,
  'M8': 115.936416972,
  'S6': 90.0,
  'Mf': 1.0980331,
  'S4': 60.0,
  
  #'Mf':          2.0*(s - zeta),
                  }

  def getAllConstituents( self ) :
    return self.angularSpeeds.keys()
    
  def getPeriod( self, constituent ) :
    """Returns constituent's period in seconds"""
    if constituent in self.angularSpeeds :
      return 1./(self.angularSpeeds[constituent]/360.0/3600)

  def getAngularSpeed( self, constituent ) :
    """Returns constituent's angular speed in deg/hour"""
    if constituent in self.angularSpeeds :
      return self.angularSpeeds[constituent]
  
class spectralPlot(object) :
  """Plots signals in frequency domain"""
  def __init__(self, frequencies=None, TMin=3*3600, TMax=32*3600, fN=2000, **kwargs) :
    if frequencies == None :
      frequencies = np.linspace( 1./TMax,1./TMin, fN )
    self.defArgs = kwargs
    self.frequencies = frequencies
    self.xlabel = kwargs.pop( 'xlabel', 'Period' )
    self.xunit = kwargs.pop( 'xunit', 'hours' )
    #self.xlabel = 'Frequency'
    #self.xunit = 'Hz'
    self.ylabel = kwargs.pop( 'ylabel', '' )
    self.titleStr = kwargs.pop( 'title', '' )
    self.ax = None

  def createAxes(self, fig, rect=111) :
    """Create rectangular diagram axes for the diagram. Returns the axes.
    """
    ax = fig.add_subplot(rect)
    return ax

  def setAxes(self, ax) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    self.ax = ax
    ax.grid()
    xStr = self.xlabel
    if self.xunit :
      xStr += ' ['+self.xunit+']'
    ax.set_xlabel(xStr)
    ax.set_ylabel(self.ylabel)

  def setupAxes(self, fig, rect=111) :
    """Creates new axes and set them.
    A shorthand for calling createAxes and setAxes.
    Returns the axes."""
    ax = self.createAxes(fig, rect)
    self.setAxes(ax)
    return ax

  def addSample(self, t, y, **kwargs) :
    """Add a sample in the diagram. kwargs are passed to pyplot.plot routine."""
    if not self.ax :
      # create default figure
      fig = plt.figure(figsize=(9,4))
      self.setupAxes(fig)
    kw = dict(self.defArgs)
    kw.update(kwargs)
    if 'xlim' in kw :
      self.ax.set_xlim( kw.pop('xlim') )
    if 'ylim' in kw :
      self.ax.set_ylim( kw.pop('ylim') )
    # TODO detrend berofe computing spectrum
    spec = computeSpectrum(t, y, self.frequencies)
    xdata = self.convertToUnits( self.frequencies )
    self.ax.semilogx(xdata, spec, **kw)
    self.updateXAxis()
    
  def updateXAxis( self ) :
    xmax = self.ax.dataLim.xmax
    xmin = self.ax.dataLim.xmin
    if 'xlim' in self.defArgs :
      xmin, xmax = self.defArgs['xlim']
    #print xmin, xmax
    decimal = -int(floor(log10(xmax)))
    #print decimal
    tickMax = round(xmax, decimal )
    decimal = -int(floor(log10(xmin)))
    #print decimal
    tickMin = round(xmin, decimal )
    #print tickMin,tickMax
    xTicks = logspace(log10(tickMin),log10(tickMax),7)
    #print xTicks
    sigDigits = 2
    xTicks = [ round(x , -int(floor(log10(xmax)))+sigDigits-1 ) for x in xTicks ]
    self.ax.set_xticks( xTicks )
    self.ax.set_xticklabels( xTicks )
    self.ax.set_xlim( [xmin, xmax] )
    from matplotlib.ticker import MaxNLocator
    self.ax.yaxis.set_major_locator(MaxNLocator(4))

  def convertToUnits( self,freq ) :
    xdata = freq
    if self.xunit == 'hours' :
      xdata = 1.0/freq/3600.
    elif self.xunit == 's' :
      xdata = 1.0/freq
    return xdata
    
  def addConstituentIndicators( self, consts, **kwargs ) :
    if not self.ax :
      return
    tcs = tidalConstituents()
    defArgs = {}
    defArgs['linestyle'] = 'dashed'
    defArgs['color'] = 'k'
    defArgs.update(kwargs)
    # try to place marker text so that they don't overlap
    ymax0 = 1.3*self.ax.dataLim.ymax
    tmp = self.ax.transData.inverted().transform([[0,0],[0,12]])
    fontSizeInDataCoords = tmp[1,1]-tmp[0,1]
    yoffset = 1.2*fontSizeInDataCoords
    xmin, xmax = self.ax.dataLim.xmin, self.ax.dataLim.xmax
    for i,c in enumerate(consts) :
      period = tcs.getPeriod(c)
      if period :
        x = self.convertToUnits(1./period)
        if x < xmin or x > xmax :
          continue # skip if outside data range
        ymin,ymax = 0, max( ymax0 - (i+1)*yoffset, ymax0/3 )
        self.ax.semilogx( [x,x], [ymin,ymax], **defArgs)
        self.ax.text( x, ymax, c, verticalalignment='top' )
    self.updateXAxis()

  def addTitle(self, titleStr=None, **kwargs) :
    """Adds a title in the axes. Optional arguments are passed to set_title routine"""
    if not titleStr :
      titleStr = self.titleStr
    if self.ax :
      self.ax.set_title(titleStr, **kwargs)

  def showLegend(self, **kwargs) :
    """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
    if self.ax :
      if not 'loc' in kwargs :
        kwargs['loc'] = 'best'
      self.ax.legend(**kwargs)

class spectralPlotDC(spectralPlot) :
  
  def addSample( self, sample, **kwargs) :
    if not 'label' in kwargs :
      kwargs['label'] = samle.description
    spectralPlot.addSample( self, sample.time.array, np.squeeze(sample.data), **kwargs)

class stackSpectralPlot(object) :
  def __init__(self, **defArgs) :
    self.defArgs = defArgs
    figsize = defArgs.pop( 'figsize' , (10,10) )
    self.fig = plt.figure( figsize=figsize )
    self.plots = dict() # stores all diagrams
    self.tags = [] # to keep ordering
    self.data = dict()
    self.legends = dict() # keep track of all legends
    self.axGrid = None

  def addPlot(self, tag, **kwargs) :
    """Adds a new subplot to the diagram"""
    self.tags.append( tag )
    kw = dict(self.defArgs)
    kw.update(kwargs)
    self.plots[tag] = spectralPlot(**kw)
    self.data[tag] = []

  def addSample( self, tag, t, y, **kwargs ) :
    if not tag in self.tags :
      self.addPlot( tag )
    self.data[tag].append( (t, y, kwargs ) )

  def isEmpty( self ) :
    for tag in self.data :
      for d in self.data[tag] :
        return False
    return True
    
  def makePlot(self,**kwargs) :
    nPlots = len(self.tags)
    if nPlots == 0 :
      return # no data
    #from mpl_toolkits.axes_grid1 import AxesGrid
    #self.axGrid = AxesGrid(self.fig, 111, # similar to subplot(111)
                    #direction='column',
                    #add_all=True,
                    #nrows_ncols = (nPlots, 1), # creates 2x2 grid of axes
                    #axes_pad=0.4, # pad between axes in inch.
                    #share_all=False,
                    #aspect=False,
                    #label_mode = "L", # tick labels, "1" (only the lower left axes), "L" (left most and bottom most axes), or "all".
                    ##cbar_location = "right", # [right|top]
                    ##cbar_mode="none", #[None|single|each]
                    ##cbar_size="7%", # size of the colorbar
                    ##cbar_pad="3%", # pad between image axes and colorbar axes
                    #)
    self.axGrid = []
    for i,tag in enumerate(self.tags) :
      if i == 0 :
        self.axGrid.append( self.fig.add_subplot(nPlots,1,i+1) )
      else :
        ax0 = self.axGrid[0]
        self.axGrid.append( self.fig.add_subplot(nPlots,1,i+1,sharex=ax0) )
    # make plots
    for i,tag in enumerate(self.tags) :
      ax = self.axGrid[i]
      self.plots[tag].setAxes( ax )
      # add all data
      for (t,y,kw) in self.data[tag] :
        self.plots[tag].addSample( t,y, **kw )
      if self.plots[tag].titleStr :
        self.plots[tag].addTitle()
      #self.plots[tag].makePlot(**kwargs)
      if i != nPlots-1 :
        ax.set_xlabel('')
        plt.setp( ax.get_xticklabels(), visible=False)
    # NOTE gridspec can handle this
    # gs1 = gridspec.GridSpec(2, 2)
    # gs1.tight_layout(fig, rect=[0, 0.03, 1, 0.97])
    self.fig.tight_layout()
    self.fig.subplots_adjust(top=0.93)

  def addTitle(self, titleStr, **kwargs) :
    """Adds a title in the axes. Optional arguments are passed to set_title routine"""
    self.fig.suptitle(titleStr, **kwargs)

  def showLegend(self, **kwargs) :
    """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
    for tag in self.tags :
      self.plots[ tag ].showLegend(prop=dict(size='small'),**kwargs)

  def addConstituentIndicators(self, consts, **kwargs ) :
    for tag in self.tags :
      self.plots[ tag ].addConstituentIndicators(consts, **kwargs)

class stackSpectralPlotDC(stackSpectralPlot) :
  
  def addSample( self, tag, sample, **kwargs ) :
    if not tag in self.tags :
      self.addPlot( tag )
    if tag not in self.data :
      self.data[tag] = []
    self.data[tag].append( (sample.time.array, np.squeeze(sample.data), kwargs ) )


if __name__=='__main__':

  from crane.data import timeArray
  from datetime import datetime

  ### examples with numpy array inputs
  # generate data
  startTime = datetime(2010,1,12,0,0,0)
  endTime = datetime(2010,2,13,3,30,0)
  dt = 90.0
  ta = timeArray.generateSimulationTimeArray(startTime,endTime,dt).asEpoch()
  t = ta.array
  # Randomly select a fraction of an array with timesteps:
  frac_points = 0.6 # Fraction of points to discard
  t = t[ np.random.rand(len(t)) >= frac_points ]
  ta = timeArray.timeArray(t, 'epoch')

  tcs = tidalConstituents()
  TM2 = tcs.getPeriod('M2')
  TS2 = tcs.getPeriod('S2')
  ref = np.sin(2*np.pi*t/TM2) + 0.95*np.sin(2*np.pi*t/TS2)
  m1 = 0.90*np.sin(2*np.pi*t/TM2)*np.sin(2*np.pi*t/TS2+0.3) - 0.08
  m2 = 0.80*np.sin(2*np.pi*t/TM2+0.8) + 0.9*np.sin(2*np.pi*t/TS2+0.5) - 0.12
  m3 = 0.78*np.sin(2*np.pi*t/TM2) + 0.82*np.sin(2*np.pi*t/TS2) + 0.03

  dia = spectralPlot(xunit='hours',fN=1000)
  dia.addSample(t,m1,label='model 1')
  dia.addSample(t,m2,label='model 2')
  dia.addSample(t,m3,label='model 3')
  dia.addTitle('spectralPlot example')
  dia.addConstituentIndicators( ['M2','S2','N2','L2','K2','K1','O1','P1','M4','M6','M8','Mf'])
  #dia.addConstituentIndicators( tcs.getAllConstituents() )
  dia.showLegend()
  plt.show()

  dia = stackSpectralPlot(xunit='hours',fN=1000)
  dia.addPlot('top',title='station 1')
  dia.addSample('top',t,m1,label='model 1')
  dia.addPlot('mid',title='station 2')
  dia.addSample('mid',t,m2,label='model 2')
  dia.addPlot('bot',title='station 3')
  dia.addSample('bot',t,m3,label='model 3')
  dia.addTitle('stackSpectralPlot example')
  dia.makePlot()
  dia.addConstituentIndicators( ['M2','S2','N2','L2','K2','K1','O1','P1','M4','M6','M8','Mf'])
  dia.showLegend()
  plt.show()
  

  ### examples with dataContainers
  d0 = dataContainer.fromTimeSeries( 'model Eins', ta, ref, ['elev'] )
  d1 = dataContainer.fromTimeSeries( 'model Eins', ta, m1, ['elev'] )
  d2 = dataContainer.fromTimeSeries( 'model Eins', ta, m2, ['elev'] )
  d3 = dataContainer.fromTimeSeries( 'model Eins', ta, m3, ['elev'] )

  dia = spectralPlotDC(xunit='hours')
  dia.addSample(d1,label='model 1')
  dia.addSample(d2,label='model 2')
  dia.addSample(d3,label='model 3')
  dia.addTitle('spectralPlotDC example')
  dia.addConstituentIndicators( ['M2','S2','N2','L2','K2','K1','O1','P1','M4','M6','M8','Mf'])
  #dia.addConstituentIndicators( tcs.getAllConstituents() )
  dia.showLegend()
  plt.show()

  dia = stackSpectralPlotDC(xunit='hours',fN=1000)
  for i in range(8) :
    dia.addPlot(str(i),title='station '+str(i))
    dia.addSample(str(i),d2,label='model 2')
  dia.addTitle('stackSpectralPlotDC example')
  dia.makePlot()
  dia.addConstituentIndicators( ['M2','S2','K1','O1','M4','M6','M8','Mf'])
  dia.showLegend()
  plt.show()
