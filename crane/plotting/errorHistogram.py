#!/usr/bin/python

"""
A class for plotting error histograms.
Uses matplotlib hist plotting routine.


Tuomas Karna 2012-09-14
"""
import numpy as np
import traceback
import sys
import matplotlib
from matplotlib.ticker import FuncFormatter
# TODO import only modules
from crane.data.dataContainer import dataContainer
from crane.plotting.plotBase import *
import crane.data.statistics as statMod

def tickLabelToPercent(y, position):
  # Ignore the passed in position. This has the effect of scaling the default
  # tick locations.
  s = str(100 * y)
  # The percent symbol needs escaping in latex
  if matplotlib.rcParams['text.usetex'] == True:
    return s + r'$\%$'
  else:
    return s + '%'

# class that takes the data (error time series) and plots the histogram
class errorHistogramBase(object) :
  """Error histogram class"""
  #def __init__(self, error='', label='', **defArgs) :
  def __init__(self, **defArgs) :
    # default plot options for all diagrams
    self.dataRange = [0,0]
    self.data = []
    self.weights = []
    self.labels = []
    self.stats = []
    self.defaultArgs = {}
    self.defaultArgs['xlabel'] = 'Error'
    self.defaultArgs['ylabel'] = 'Frequency'
    self.defaultArgs['unit'] = ''
    self.defaultArgs['bins'] = 30
    self.defaultArgs['normed'] = False
    self.defaultArgs['percentage'] = False
    self.defaultArgs['histtype'] = 'bar' #  [ 'bar' | 'barstacked' | 'step' | 'stepfilled' ]
    self.defaultArgs['transparent'] = False
    #self.defaultArgs['histtype'] = 'step' #  [ 'bar' | 'barstacked' | 'step' | 'stepfilled' ]
    #self.defaultArgs['linewidth'] = 2 #  [ 'bar' | 'barstacked' | 'step' | 'stepfilled' ]
    self.defaultArgs['shadedLim'] = []
    # override with user input
    self.defaultArgs.update( defArgs )
    self.unit = self.defaultArgs.pop('unit')
    self.xlabel = self.defaultArgs.pop('xlabel')
    self.ylabel = self.defaultArgs.pop('ylabel')
    self.ax = None
    self.percentage = self.defaultArgs.pop('percentage')
    if self.percentage :
      self.defaultArgs['normed']=False
    # default color sequence for samples (typ 'r' for observation)
    self.defColors = ['r', 'b','g','k','m','c']
    self.colors = []
    
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
    errStr = self.xlabel
    if self.unit :
      errStr += ' ['+self.unit+']'
    ax.set_xlabel(errStr)
    ax.set_ylabel(self.ylabel)
  
  def setupAxes(self, fig, rect=111) :
    """Creates new axes and set them.
    A shorthand for calling createAxes and setAxes.
    Returns the axes."""
    ax = self.createAxes(fig, rect)
    self.setAxes(ax)
    return ax
  
  def addSample(self, error, label, color=None, weights=None, stats=None) :
    """Add an error signal to the diagram, identified with the given label string."""
    self.data.append( error )
    self.stats.append( stats )
    self.labels.append( label )
    if weights != None :
      self.weights.append( weights )
    elif self.percentage :
      self.weights.append( np.ones_like(error)/len(error) )
    if not color :
      color = self.defColors[len(self.data)-1]
    self.colors.append( color )
    self.dataRange = [ min(self.dataRange[0],error.min()) , max(self.dataRange[1],error.max()) ]
    
  def makePlot(self, **kwargs) :
    """Create histogram for all the error data. Same bins are used for all the data sets.
    kwargs is passed to matplotlib hist command."""
    if not self.ax :
      # create default figure
      fig = plt.figure(figsize=(9,4))
      self.setupAxes(fig)
    if not self.data :
      return # in case no error could be computed
    kw = dict(self.defaultArgs) # copy default
    kw.update(kwargs)           # override
    plotMinMax = kw.pop('plotMinMax',False)
    shadedLim = kw.pop('shadedLim',None)
    if len(self.weights) == len(self.data) :
      kw['weights'] = self.weights
    kw['color'] = self.colors
    kw['edgecolor'] = 'none'
    if 'title' in kw :
      self.addTitle( kw.pop('title') )
    if not 'range' in kw :
      kw['range'] = self.dataRange
    # set symmetric x axis
    xmax = 1.02 * max( abs(self.dataRange[0]), abs(self.dataRange[1]) )
    xlim = kw.pop('xlim', [-xmax, xmax])
    if shadedLim :
      self.ax.axvspan(shadedLim[0], shadedLim[1], facecolor=[0.8,1.0,0.85], edgecolor='none')
    transparent = kw.pop('transparent', False)
    if not transparent:
      kw['label'] = self.labels
      n, bins, patches = self.ax.hist(self.data, **kw)
    else:
      for i in range(len(self.data)):
        kw['color'] = self.colors[i]
        kw.setdefault('alpha', 0.4)
        kw['label'] = self.labels[i]
        if len(self.weights) == len(self.data):
          kw['weights'] = self.weights[i]
        n, bins, patches = self.ax.hist(self.data[i], **kw)
    # add space on top
    axheight_px = self.ax.bbox.extents[3]-self.ax.bbox.extents[1] # in px
    yspan = self.ax.get_ylim()[1]-self.ax.get_ylim()[0]
    ydataspan = self.ax.dataLim.ymax-self.ax.dataLim.ymin
    clearance_px = min( 20, 0.15*axheight_px )
    new_yspan = ydataspan/(1.0-(clearance_px/axheight_px))
    new_ymax = new_yspan-self.ax.get_ylim()[0]
    ymax = max( self.ax.get_ylim()[1], new_ymax )
    self.ax.set_ylim([self.ax.get_ylim()[0],ymax])
    if self.percentage :
      percent_formatter = FuncFormatter(tickLabelToPercent)
      self.ax.yaxis.set_major_formatter(percent_formatter)

    if plotMinMax :
      # add min/mean/max vertical lines
      ymax = self.ax.dataLim.ymax * 1.3
      meanStr = ''
      rmseStr = ''
      nmseStr = ''
      for i in range(len(self.data)) :
        d = self.data[i]
        x = np.percentile(d,5) #d.min()
        self.ax.axvspan(x, x, color='w', linewidth=1.0 )
        self.ax.axvspan(x, x, color=self.colors[i], linewidth=1.0, linestyle='solid' )
        x = np.percentile(d,95) #d.max()
        self.ax.axvspan(x, x, color='w', linewidth=1.0 )
        self.ax.axvspan(x, x, color=self.colors[i], linewidth=1.0, linestyle='solid' )
        x = d.mean()
        self.ax.axvspan(x, x, color='w', linewidth=1.5 )
        self.ax.axvspan(x, x, color=self.colors[i], linewidth=2.5, linestyle='solid' )

        rmse = np.sqrt(np.mean(d**2))
        if self.stats[i] != None :
          murphy = self.stats[i]['murphy']
          ioa = self.stats[i]['ioa']
          nmse = self.stats[i]['nmse']
          meanStr += '{0:4.2f} '.format(d.mean())
          rmseStr += '{0:4.2f} '.format(rmse)
          try:
            nmseStr += '{0:4.2f} '.format(nmse)
          except Exception as e:
            import pdb; pdb.set_trace()
        else :
          meanStr += '{0:4.2f} '.format(d.mean())
          rmseStr += '{0:4.2f} '.format(rmse)
    if nmseStr != '':
      statsText = u'BIAS: {0:s}{unit:s}\n RMSE: {1:s}{unit:s}\n NMSE: {2:s}'.format(meanStr, rmseStr, nmseStr, unit=self.unit)
    else:
      statsText = u'BIAS: {0:s}{unit:s}\n RMSE: {1:s}{unit:s}'.format(meanStr, rmseStr, unit=self.unit)
    self.ax.text(0.98, 0.9, statsText, horizontalalignment='right',
                 verticalalignment='top', transform=self.ax.transAxes)
    from matplotlib.ticker import MaxNLocator
    self.ax.yaxis.set_major_locator(MaxNLocator(4))
    self.ax.set_xlim(xlim)
    return bins, patches

  def addTitle(self, titleStr, **kwargs) :
    """Adds a title in the axes. Optional arguments are passed to set_title routine"""
    if self.ax :
      self.ax.set_title(titleStr, **kwargs)

  def showLegend(self, *args, **kwargs) :
    """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
    if self.ax :
      if not 'loc' in kwargs :
        kwargs['loc'] = 'best'
      self.ax.legend(*args, **kwargs)

# class that takes the error signal as a dataContainer
class errorHistogramDC(errorHistogramBase) :
  """A class where the signals are dataContainer objects."""
  def __init__(self, reference, **defArgs) :
    """Create new instance with the given reference signal"""
    errorHistogramDC.checkData( reference )
    self.referenceDC = reference
    errorHistogramBase.__init__(self, **defArgs)
    
  def addSample(self, sample, label, color=None) :
    """Add a sample in the error histogram plot. Error is computed versus the reference signal.
    If the time stamps are different, sample is interpolated on the reference time steps."""
    errorHistogramDC.checkData( sample )
    try :
      err = self.referenceDC.computeError( sample )
      r,o = self.referenceDC.alignTimes( sample )
      stats = statMod.computeSkillStatistics( r.data.ravel(), o.data.ravel() )
      errorHistogramBase.addSample(self,np.squeeze(err.data),label, color, stats=stats)
    except Exception as e:
      print 'Error could not be computed, skipping:', label, sample.description
      traceback.print_exc(file=sys.stdout)
      #print e

  @staticmethod
  def checkData( sample ) :
    if not isinstance( sample, dataContainer ) :
      raise Exception( 'sample must be a dataContainer object' )
    if sample.data.shape[0] > 1 :
      raise Exception( 'spatially varying data is not supported' )
    if sample.data.shape[1] > 1 :
      raise Exception( 'multiple fields are not supported' )
    return True

class stackHistogram(object) :
  def __init__(self, **defArgs) :
    self.defArgs = defArgs
    figsize = defArgs.pop( 'figsize' , (10,10) )
    self.fig = plt.figure( figsize=figsize )
    self.plots = dict() # stores all histograms
    self.tags = [] # to keep ordering
    self.legends = dict() # keep track of all legends
    self.axGrid = None

  def addPlot(self, tag, **kwargs) :
    """Adds a new subplot to the diagram"""
    self.tags.append( tag )
    kw = dict(self.defArgs)
    kw.update(kwargs)
    self.plots[tag] = errorHistogramBase(**kw)

  def addSample( self, tag, sample, label, color=None, stats=None ) :
    if not tag in self.tags :
      self.addPlot( tag )
    self.plots[tag].addSample( sample, label, color, stats=stats )

  def isEmpty( self ) :
    for tag in self.plots :
      for err in self.plots[tag].data :
        return False
    return True

  def makePlot(self,**kwargs) :
    nPlots = len(self.tags)
    if nPlots == 0 :
      return # no data
    from mpl_toolkits.axes_grid1 import AxesGrid
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
    # unify data ranges
    rmin =  10e99
    rmax = -10e99
    for i,tag in enumerate(self.tags) :
      r = self.plots[tag].dataRange
      rmin,rmax = min(rmin,r[0]) , max(rmax, r[1])
    for i,tag in enumerate(self.tags) :
      self.plots[tag].dataRange = [rmin,rmax]
    # set symmetric x axis
    xmax = 1.02* max( abs(rmin), abs(rmax) )
    xlim = kwargs.pop('xlim',[-xmax,xmax])
    # make plots
    from matplotlib.ticker import MaxNLocator
    for i,tag in enumerate(self.tags) :
      ax = self.axGrid[i]
      self.plots[tag].setAxes( ax )
      self.plots[tag].makePlot(**kwargs)
      if i != nPlots-1 :
        ax.set_xlabel('')
        plt.setp( ax.get_xticklabels(), visible=False)
      ax.yaxis.set_major_locator(MaxNLocator(3))
    self.axGrid[0].set_xlim(xlim)
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

class stackHistogramDC(stackHistogram) :
  def __init__(self, **defArgs) :
    """Create new instance with the given reference signal"""
    self.references = dict()
    stackHistogram.__init__(self, **defArgs)
    
  def addPlot(self, tag, reference, **kwargs) :
    errorHistogramDC.checkData( reference )
    self.references[tag] = reference
    stackHistogram.addPlot( self, tag, **kwargs )
    
  def addSample( self, tag, sample, label, color=None) :
    if not tag in self.references :
      raise Exception('Given tag not found')
    errorHistogramDC.checkData( sample )
    try :
      err = self.references[tag].computeError( sample )
      r,o = self.references[tag].alignTimes( sample )
      stats = statMod.computeSkillStatistics( r.data.ravel(), o.data.ravel() )
      stackHistogram.addSample(self,tag,np.squeeze(err.data),label, color, stats=stats)
    except Exception as e:
      print 'Error could not be computed, skipping:', tag, label, sample.description
      print e

if __name__=='__main__':

  from crane.data import timeArray
  from datetime import datetime

  ### examples with numpy array inputs
  # generate data
  startTime = datetime(2010,1,12,0,0,0)
  endTime = datetime(2010,2,13,3,30,0)
  dt = 900.0
  ta = timeArray.generateSimulationTimeArray(startTime,endTime,dt).asEpoch()
  t = ta.array
  
  T = 44714
  ref = np.sin(t/T) + 0.95*np.sin(0.95*t/T)
  m1 = 0.90*np.sin(t/T) + 0.7*np.sin(0.85*t/T+0.3) - 0.08
  m2 = 0.80*np.sin(t/T+0.8) + 0.9*np.sin(0.95*t/T+0.5) - 0.12
  m3 = 0.78*np.sin(t/T) + 0.82*np.sin(0.90*t/T) + 0.03
  
  dia = errorHistogramBase(unit='m',shadedLim=[-0.3,0.3])
  dia.addSample(m1-ref,label='model 1')
  dia.addSample(m2-ref,label='model 2')
  dia.addSample(m3-ref,label='model 3')
  dia.makePlot(bins=20)
  dia.addTitle('errorHistogramBase example')
  dia.showLegend()
  plt.show()

  dia = stackHistogram(unit='m',shadedLim=[-0.3,0.3])
  dia.addPlot('top',title='station 1')
  dia.addSample('top',m1,'model 1')
  dia.addSample('top',m2,'model 2')
  dia.addPlot('two',title='station 2')
  dia.addSample('two',m3,'model 3',color='c')
  dia.addPlot('drei',title='station 3')
  dia.addSample('drei',m2,'model 2')
  dia.addSample('drei',m3,'model 3')
  dia.makePlot()
  dia.addTitle('Stack histogram example')
  dia.showLegend()
  plt.show()
  
  ### example with dataContainer objects
  # generate data
  np.random.seed(int(3411))
  startTime = datetime(2010,1,12,0,0,0)
  startCorie = datetimeToCorieTime(startTime)
  #t0 = np.hstack( ( np.linspace(0,12,20), np.linspace(15.33,30,15) ) ) + startCorie
  t0 = np.linspace(0,12,20) + startCorie
  m0 = np.sin(t0)
  ta0 = timeArray.timeArray(t0,'corie')
  
  t1 = np.linspace(-10,34,100) + startCorie
  m1 = 0.8*np.sin(t1) + 0.03*np.random.randn(len(t1))
  ta1 = timeArray.timeArray(t1,'corie')
  
  t2 = np.linspace(-9,31.7,65) + startCorie
  m2 = 0.95*np.sin(t2) + 0.12*np.random.randn(len(t2))
  ta2 = timeArray.timeArray(t2,'corie')
  
  t3 = np.linspace(-9,32.2,100) + startCorie
  m3 = np.sin(t3-0.12) - 0.05*np.random.randn(len(t3))-3.0
  ta3 = timeArray.timeArray(t3,'corie')
  
  d0 = dataContainer.fromTimeSeries( 'Observation', ta0, m0, ['elev'] )
  d1 = dataContainer.fromTimeSeries( 'model Eins', ta1, m1, ['elev'] )
  d2 = dataContainer.fromTimeSeries( 'model Zwei', ta2, m2, ['elev'] )
  d3 = dataContainer.fromTimeSeries( 'model Drei', ta3, m3, ['elev'] )
 
  dia = errorHistogramDC(d0, unit='m',shadedLim=[-0.3,0.3],plotMinMax=True)
  dia.addSample(d1,label='model 1', color='b')
  dia.addSample(d2,label='model 2', color='y')
  dia.addSample(d3,label='model 3', color='m')
  dia.makePlot(bins=30)
  dia.addTitle('errorHistogramDC example')
  dia.showLegend()
  plt.show()

  dia = stackHistogramDC(unit='m',figsize=(10,10),normed=True,shadedLim=[-0.3,0.3],plotMinMax=True)
  dia.addPlot('top',d0,title='station 1',ylabel='freq')
  dia.addSample('top',d1,'model 1')
  dia.addSample('top',d2,'model 2')
  dia.addPlot('two',d0, title='station 2',ylabel='freq')
  dia.addSample('two',d2,'model 2')
  dia.addSample('two',d2,'model 3')
  dia.addPlot('drei',d0, title='station 3',ylabel='freq')
  dia.addSample('drei',d2,'model 2')
  dia.addSample('drei',d3,'model 3')
  dia.addPlot('4',d0, title='station 4',ylabel='freq')
  dia.addSample('4',d2,'model 2')
  dia.addSample('4',d3,'model 3')
  dia.addPlot('5',d0, title='station 5',ylabel='freq')
  dia.addSample('5',d2,'model 2')
  dia.addSample('5',d3,'model 3')
  dia.addPlot('6',d0, title='station 6',ylabel='freq')
  dia.addSample('6',d2,'model 2')
  dia.addSample('6',d3,'model 3')
  dia.addPlot('7',d0, title='station 7',ylabel='freq')
  dia.addSample('7',d2,'model 2')
  dia.addSample('7',d3,'model 3')
  dia.addPlot('8',d0, title='station 8',ylabel='freq')
  dia.addSample('8',d2,'model 2')
  dia.addSample('8',d3,'model 3')
  dia.makePlot(bins=20)
  dia.addTitle('stackHistogramDC example')
  dia.showLegend()
  plt.show()
  