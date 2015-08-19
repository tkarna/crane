#!/usr/bin/python

"""
Implementation of Taylor and bias vs. index-of-agreement diagrams.

For information on Taylor diagram, see:
http://www-pcmdi.llnl.gov/about/staff/Taylor/CV/Taylor_diagram_primer.htm

Tuomas Karna 2012-08-15
"""

from scipy import *
import sys

# TODO import only modules
from crane.data import statistics
from crane.data.dataContainer import dataContainer
from crane.plotting.plotBase import *

class diagramBase(object) :
  """Base class for statistics diagrams"""
  def __init__(self, refStd, refLabel='Reference', normalized=False, **defArgs) :
    """Creates a new diagram.
    refStd is the standard deviation of the reference (data) sample to be compared to.
    If normalized==True, dimensional measures (e.g. std and bias) are
    normalized by dividing by std of the reference signal."""
    self.refLabel = refLabel
    self.normalized = normalized
    self.refStd = refStd
    # default plot options for all diagrams
    self.defaultArgs = {}
    self.defaultArgs['marker'] = '^'
    self.defaultArgs['markersize'] = 7
    self.defaultArgs['linestyle'] = 'None'
    self.defaultArgs['unit'] = ''
    # override with user input
    self.defaultArgs.update( defArgs )
    # remove unit from the arg list (doesn't work with plotting commands)
    self.unit = self.defaultArgs.pop('unit')

    # default color sequence for samples
    self.colors = ['b','g','k','m','c']
    # number of plotted datasets
    self.nColors = 0
    
  def createAxes(self, fig, rect=111) :
    """Creates new axes for the plot in the given figure. Must be defined for each derived type.
    Returns the axes and parent axes (None if parent axes do not exist)."""
    raise NotImplementedError, "This method must be implemented in the derived class"

  def setAxes(self, ax, parent_ax=None) :
    """Set axes for the diagram. All data will be plotted in these axes. Must be defined for each derived type."""
    raise NotImplementedError, "This method must be implemented in the derived class"
  
  def plotReference(self) :
    """Plots the reference signal and associated metrics."""
    raise NotImplementedError, "This method must be implemented in the derived class"

  def setupAxes(self, fig, rect=111) :
    """Creates new axes and set them, and plots the reference.
    A shorthand for calling createAxes, setAxes and plotReference.
    Returns the axes and parent axes (None if parent axes do not exist)."""
    ax, pax = self.createAxes(fig, rect)
    self.setAxes(ax, pax)
    self.plotReference()
    return ax, pax
  
  def plotSample(self, x, y, **kwargs):
    """Add a sample to the diagram. Optional arguments are passed to plot routine."""
    kw = dict(self.defaultArgs) # copy default
    kw.update(kwargs)           # override
    if 'color' not in kw :
      kw['color'] = self.colors[self.nColors]
      self.nColors += 1
    kw.setdefault('zorder',5)
    l, = self.ax.plot(x, y, **kw)
    return l

  def addTitle(self, titleStr, **kwargs) :
    """Adds a title in the axes. Optional arguments are passed to set_title routine"""
    self.ax.set_title(titleStr, **kwargs)

  def showLegend(self, *args, **kwargs) :
    """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
    self.ax.legend(numpoints=1, *args, **kwargs)

class TaylorDiagram(diagramBase) :
  """Taylor diagram: plot model standard deviation and correlation
  versus reference data in a single-quadrant polar plot, with
  r=stddev and theta=arccos(correlation).
  """
  def __init__(self, refStd, refLabel, normalized, **defArgs):
    diagramBase.__init__(self, refStd, refLabel, normalized, **defArgs)
    # std axis range
    maxStd = 1.4*self.refStd if not self.normalized else 1.25
    self.maxStd = self.defaultArgs.pop('taylorMaxStd', maxStd)

  def createAxes(self, fig, rect=111) :
    """Create Taylor diagram axes, i.e. single quadrant polar
    plot, using mpl_toolkits.axisartist.floating_axes.
    """
    from matplotlib.projections import PolarAxes
    import mpl_toolkits.axisartist.floating_axes as FA
    import mpl_toolkits.axisartist.grid_finder as GF

    tr = PolarAxes.PolarTransform()

    # Correlation labels
    self.minCorr = 0.0
    rlocs = concatenate((arange(self.minCorr,1.0,0.1),[0.95,0.99]))
    tlocs = arccos(rlocs)        # Conversion to polar angles
    gl1 = GF.FixedLocator(tlocs) # Correlation tick/grid positions
    tf1 = GF.DictFormatter(dict(zip(tlocs, map(str,rlocs))))
    
    # TODO fix std ticks into more meaningful positions (now maxStd/10)
    # TODO BUG OverflowError: cannot convert float infinity to integer
    decimal = -int(floor(log10(self.maxStd)))
    stdTickMax = round(self.maxStd, decimal )
    digit = stdTickMax/10**int(floor(log10(stdTickMax)))
    nBaseInterv = max( 4/digit, 1.0 )
    nTicks = int(digit*nBaseInterv)+1
    tickInterval = stdTickMax/(nTicks-1)
    stdTicks = arange(0,stdTickMax+3*tickInterval,tickInterval)
    
    gl2 = GF.FixedLocator(stdTicks) # Std tick/grid positions
    #GridHelperCurveLinear(aux_trans, extreme_finder=None, grid_locator1=None,
    #    grid_locator2=None, tick_formatter1=None, tick_formatter2=None)
    ghelper = FA.GridHelperCurveLinear(tr,
                                        extremes=(0,arccos(self.minCorr),0,self.maxStd),
                                        grid_locator1=gl1,
                                        tick_formatter1=tf1,
                                        grid_locator2=gl2
                                        )
    
    # parent axes
    pax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
    fig.add_subplot(pax)
    # Adjust axes
    pax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
    pax.axis["top"].toggle(ticklabels=True, label=True)
    pax.axis["top"].major_ticklabels.set_axis_direction("top")
    pax.axis["top"].label.set_axis_direction("top")
    pax.axis["top"].label.set_text("Correlation")

    stdStr = 'Standard deviation'
    if self.unit and not self.normalized :
      stdStr += ' ['+self.unit+']'
    pax.axis["left"].set_axis_direction("bottom") # "X axis"
    pax.axis["left"].label.set_text(stdStr)
    #pax.axis["left"].toggle(ticklabels=False)

    pax.axis["right"].set_axis_direction("top")   # "Y axis"
    pax.axis["right"].toggle(ticklabels=True)
    #pax.axis["right"].label.set_text("Standard deviation")
    #pax.axis["right"].toggle(label=True)
    pax.axis["right"].major_ticklabels.set_axis_direction("left")

    pax.axis["bottom"].set_visible(False)         # Useless
    
    # Grid
    pax.grid()
    # polar coordinate axes
    ax = pax.get_aux_axes(tr)
    
    return ax, pax 

  def setAxes(self, ax, parent_ax) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    #self._ax = ax                   # Graphical axes
    #self.ax = ax.get_aux_axes(tr)   # Polar coordinates
    self.ax = ax                     # Polar coordinates
    self.parent_ax = parent_ax       # Graphical axes

  def plotReference(self) :
    """Adds reference point in the diagram"""
    corr = 1.0
    stddev = self.refStd
    self.plotSample(corr,stddev,color='r',marker='o',markersize=10,label=self.refLabel)
    # draw STD contour
    t = linspace(0,arccos(self.minCorr))
    r = zeros_like(t)
    if self.normalized :
      r += 1.0
    else :
      r += stddev
    self.ax.plot(t,r,color='r',linestyle='solid', label='_nolegend_')
    # plot RMSE circles
    decimal = -int(floor(log10(self.maxStd)))
    roundedMaxStd = round(self.maxStd, decimal )
    rmsContours = array([0.1,0.25,0.5,0.75])*roundedMaxStd
    for rms in rmsContours :
      self.addRmsCircle(rms,linestyle='dashed',color='green',text='')
    
  def plotSample(self, corr, stddev, **kwargs) :
    """Adds a sample in the diagram.
    corr -- correlation cofficient versus the reference signal
    stddev -- standard deviation of the sample
    """
    stddev = stddev
    if self.normalized and self.refStd > 1e-6:
        stddev = stddev/self.refStd
    l = diagramBase.plotSample(self,arccos(corr),stddev, **kwargs)
    return l 

  def addRmsCircle(self,rms, text=None, **kwargs) :
    """Add constant RMS curve in the Taylor Diagram"""
    refStd = self.refStd if not self.normalized else 1.0
    N = 20
    TOL = 1e-4
    rMin = refStd - rms + TOL
    if rms**2 - refStd**2 + self.minCorr**2*refStd**2 > 0 :
      # touches left boundary
      rMin = self.minCorr*refStd + sqrt(rms**2 - refStd**2 + self.minCorr**2*refStd**2)+TOL
    maxStd = self.parent_ax.get_xlim()[1]
    rMax = min(refStd + rms, maxStd-TOL) # right boundary
    r = linspace( rMin, rMax, N )
    t = ( r**2 + refStd**2 - rms**2 ) / (2*r*refStd)
    t = arccos(t)
    self.ax.plot(t,r, **kwargs)
    c = 'k'
    if 'color' in kwargs :
      c = kwargs['color']
    labelstr = str(rms)
    if text and len(text)>0 :
      labelstr = text+' '+str(rms)
    self.ax.annotate(labelstr,xy=(t[10],r[10]),xytext=(t[10],r[10]+0.1*rms),rotation=0,color=c)
        
  def addTitle(self, titleStr, **kwargs) :
    """Adds a title in the axes. Arguments are passed to set_title routine"""
    self.parent_ax.set_title(titleStr, **kwargs)



class BiasDiagram(diagramBase) :
  """Bias diagram: plot model index-of-agreement versus bias.
  """
  def __init__(self, refStd, refLabel='Reference', normalized=False,
               normalizeBias=False, **defArgs) :
    self.normalizeBias = normalizeBias
    diagramBase.__init__(self,refStd, refLabel, normalized, **defArgs)
    self.defaultArgs.pop('taylorMaxStd')

  def createAxes(self, fig, rect=111) :
    """Create rectangular diagram axes for the diagram.
    """
    ax = fig.add_subplot(rect)
    ax.grid()
    biasStr = 'Bias'
    if self.unit and not self.normalizeBias :
      biasStr += ' ['+self.unit+']'
    ax.set_xlabel(biasStr)
    ax.set_ylabel('Index of agreement')
    ax.set_ylim([0,1.01])
    return ax, None
    
  def setAxes(self, ax, parent_ax=None) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    self.ax = ax

  def plotReference(self) :
    # Add reference point
    bias = 0.0
    ioa = 1.0
    self.plotSample(bias,ioa,color='r',marker='o',linestyle='None',markersize=10,label=self.refLabel)
    # add bias=0 contour
    self.ax.plot([0,0],[0,1],color='r',linestyle='solid',label='_nolegend_')

  def plotSample(self, bias, ioa, **kwargs):
    """Adds a sample in the diagram.
    bias -- bias versus the reference signal
    ioa  -- index of agreement versus the reference
    """
    l = diagramBase.plotSample(self,bias,ioa, **kwargs)
    # set good axis bounds
    xMax = max( abs(self.ax.dataLim.xmin), abs(self.ax.dataLim.xmax) )
    if xMax != 0 :
      # round to next full decimal
      decimal = int(floor(log10(xMax)))
      xTickMax = math.ceil(xMax * 10**-decimal) / 10**-decimal
      xTicks = linspace(-xTickMax,xTickMax,7)
      self.ax.set_xticks( xTicks )
      self.ax.set_xticklabels( [ str( round(xt, -decimal+1) ) for xt in xTicks ] )
      self.ax.set_aspect(2*xTickMax)
    return l

class BiasMurphyDiagram(diagramBase) :
  """Bias diagram: plot model Murphy score versus bias.
  """
  def __init__(self, refStd, refLabel='Reference', normalized=False,
               normalizeBias=False, **defArgs) :
    self.normalizeBias = normalizeBias
    diagramBase.__init__(self,refStd, refLabel, normalized, **defArgs)
    self.defaultArgs.pop('taylorMaxStd', None)
    self.dataLim = [-1.01, 1.01]

  def createAxes(self, fig, rect=111) :
    """Create rectangular diagram axes for the diagram.
    """
    ax = fig.add_subplot(rect)
    ax.grid()
    biasStr = 'Bias'
    if self.unit and not self.normalizeBias :
      biasStr += ' ['+self.unit+']'
    ax.set_xlabel(biasStr)
    ax.set_ylabel('Murphy score')
    ax.set_ylim([-1.01,1.01])
    return ax, None

  def setAxes(self, ax, parent_ax=None) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    self.ax = ax

  def plotReference(self) :
    # Add reference point
    bias = 0.0
    mur = 1.0
    self.plotSample(bias,mur,color='r',marker='o',linestyle='None',markersize=10,label=self.refLabel)
    # add bias=0 contour
    self.ax.axvline(0.0, color='r', linestyle='solid', label='_nolegend_')
    self.shading = self.ax.axhspan(-15.0, 0.0, edgecolor='none', facecolor='Gainsboro')

  def plotSample(self, bias, murphy, **kwargs):
    """Adds a sample in the diagram.
    bias -- bias versus the reference signal
    murphy  -- murphy skill score versus the reference
    """
    l = diagramBase.plotSample(self, bias, murphy, **kwargs)
    # set good axis bounds
    xMax = max( abs(self.ax.dataLim.xmin), abs(self.ax.dataLim.xmax) )
    if xMax != 0 :
      # round to next full decimal
      decimal = int(floor(log10(xMax)))
      xTickMax = math.ceil(xMax * 10**-decimal) / 10**-decimal
      xTicks = linspace(-xTickMax,xTickMax,7)
      self.ax.set_xticks( xTicks )
      self.ax.set_xticklabels( [ str( round(xt, -decimal+1) ) for xt in xTicks ] )
      self.dataLim[0] = max(min(self.dataLim[0], l.get_ydata().min()), -10.0)
      ylim = self.dataLim
      if ylim[0] < -1.1:
          ylim[0] *= 1.05
      self.ax.set_ylim(ylim)
      yrange = self.dataLim[1]-self.dataLim[0]
      self.ax.set_aspect(2*xTickMax/yrange)
    return l

class BiasNMSEDiagram(diagramBase) :
  """Bias diagram: plot model Normalized Mean Square Error versus bias.
  """
  def __init__(self, refStd, refLabel='Reference', normalized=False,
               normalizeBias=False, **defArgs) :
    self.normalizeBias = normalizeBias
    diagramBase.__init__(self,refStd, refLabel, normalized, **defArgs)
    self.defaultArgs.pop('taylorMaxStd', None)
    self.dataLim = [-1.01, 1.01]

  def createAxes(self, fig, rect=111) :
    """Create rectangular diagram axes for the diagram.
    """
    ax = fig.add_subplot(rect)
    ax.grid()
    biasStr = 'Bias'
    if self.unit and not self.normalizeBias :
      biasStr += ' ['+self.unit+']'
    ax.set_xlabel(biasStr)
    ax.set_ylabel('NMSE')
    ax.set_ylim([0.0,1.51])
    return ax, None

  def setAxes(self, ax, parent_ax=None) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    self.ax = ax

  def plotReference(self) :
    # Add reference point
    bias = 0.0
    mur = 0.0
    self.plotSample(bias,mur,color='r',marker='o',linestyle='None',markersize=10,label=self.refLabel)
    # add bias=0 contour
    self.ax.axvline(0.0, color='r', linestyle='solid', label='_nolegend_')
    self.shading = self.ax.axhspan(1.0, 15.0, edgecolor='none', facecolor='Gainsboro')

  def plotSample(self, bias, nmse, **kwargs):
    """Adds a sample in the diagram.
    bias -- bias versus the reference signal
    nmse  -- normalized mean square error against the reference
    """
    l = diagramBase.plotSample(self, bias, nmse, **kwargs)
    # set good axis bounds
    xMax = max( abs(self.ax.dataLim.xmin), abs(self.ax.dataLim.xmax) )
    if xMax != 0 :
      # round to next full decimal
      decimal = int(floor(log10(xMax)))
      xTickMax = math.ceil(xMax * 10**-decimal) / 10**-decimal
      xTicks = linspace(-xTickMax,xTickMax,7)
      self.ax.set_xticks( xTicks )
      self.ax.set_xticklabels( [ str( round(xt, -decimal+1) ) for xt in xTicks ] )
      self.dataLim[1] = min(max(self.dataLim[1], l.get_ydata().max()), 10.0)
      ylim = self.dataLim
      if ylim[1] > 1.0:
          ylim[1] *= 1.05
      self.ax.set_ylim(ylim)
      yrange = self.dataLim[1]-self.dataLim[0]
      self.ax.set_aspect(2*xTickMax/yrange)
    return l

class statisticsDiagram(object) :
  """A class that plots both the Taylor and bias diagram in a same figure."""
  def __init__(self, refStd, refLabel='Reference', **defArgs) :
    """refStd is the standard deviation of the reference to be compared to."""
    self.fig = plt.figure(figsize=(9,5))
    self.fig.subplots_adjust(top=1.00,bottom=0.20) # for (9,5) fig
    self.taylor = TaylorDiagram(refStd, refLabel, **defArgs)
    #self.bias = BiasMurphyDiagram(refStd, refLabel, **defArgs)
    self.bias = BiasNMSEDiagram(refStd, refLabel, **defArgs)
    
    self.taylor.setupAxes(self.fig, 121)
    self.bias.setupAxes(self.fig, 122)
    
  @classmethod
  def fromStats(cls, stats, **defArgs) :
    """Creates a statistics diagram from given skillAssessmentStats object."""
    # create diagram object
    refStd = stats[(stats.refLabel,'stddev')]
    obj = cls(refStd, stats.refLabel, **defArgs)
    # add all samples
    labels = set()
    for key in stats.keys() :
      if key[0] != stats.refLabel :
        labels.add( key[0] )
        
    for l in labels :
      obj.plotSample(stats[(l,'corr')], stats[(l,'stddev')], stats[(l,'bias')], stats[(l,'murphy')],label=l)
      
    return obj
    
  def plotSample(self,corr,stddev,bias,murphy, **kwargs) :
    """Adds a sample in the diagram.
    corr -- correlation cofficient versus the reference signal
    stddev -- standard deviation of the sample
    bias -- bias versus the reference signal
    murphy  -- murphy skill score versus the reference
    """
    l1 = self.taylor.plotSample(corr,stddev, **kwargs)
    l2 = self.bias.plotSample(bias,murphy, **kwargs)
    return l1,l2
  
  def showLegend(self, **kwargs) :
    kw = dict()
    kw['prop'] = {'size':8}
    kw['numpoints'] = 1
    if not 'loc' in kwargs :
      kw['loc'] = 'upper center'
      kw['bbox_to_anchor'] = (-0.1, -0.15)
      kw['ncol'] = 4
    kw.update(kwargs)
    self.bias.ax.legend(**kw)
    
  def addTitle(self, titleStr, **kwargs) :
    self.fig.suptitle(titleStr, **kwargs)

class statisticsDiagramArray(statisticsDiagram, statistics.skillAssessmentStats) :
  """A version that takes signals as arguments and computes statistics automatically."""
  def __init__(self, reference, refLabel='Reference', **defArgs) :
    statistics.skillAssessmentStats.__init__(self, reference, refLabel)
    refStd = self.stats[(refLabel,'stddev')]
    statisticsDiagram.__init__(self, refStd, refLabel, **defArgs)
  
  def plotSample( self, sample, reference=None, **kwargs ) :
    if 'label' not in kwargs.keys() :
      raise Exception( 'keyword argument label must be provided' )
    l = kwargs['label']
    corr, stddev, bias, ioa, mur, crmse = self.addSample( sample, l, reference=reference)
    l1,l2 = statisticsDiagram.plotSample(self, corr, stddev, bias, mur, **kwargs)
    return l1,l2

class statisticsDiagramDC(statisticsDiagramArray) :
  """A wrapper class that uses dataContainer objects as inputs.
  Input data is automatically cast in same time units and
  interpolated to the reference time instances."""
  def __init__(self, refsample, refLabel='Reference', **defArgs) :
    statisticsDiagramDC.checkData( refsample )
    self.dataContainer = refsample
    statisticsDiagramArray.__init__( self, squeeze(refsample.data), refLabel, **defArgs )
  
  @staticmethod
  def checkData( sample ) :
    if not isinstance( sample, dataContainer ) :
      raise Exception( 'sample must be a dataContainer object' )
    if sample.data.shape[0] > 1 :
      raise Exception( 'spatially varying data is not supported' )
    if sample.data.shape[1] > 1 :
      raise Exception( 'multiple fields are not supported' )
    return True
    
  def plotSample(self, sample, **kwargs) :
    """Add a sample to the  diagram. args and kwargs are
    directly propagated to the plot command."""
    statisticsDiagramDC.checkData( sample )
    if not 'label' in kwargs :
      kwargs['label'] = sample.description
    s = sample
    r = self.dataContainer
    if self.dataContainer.time != sample.time :
      try :
        r,s = self.dataContainer.alignTimes( sample )
      except Exception as e:
        print 'Time series could not be aligned, skipping:', kwargs['label'], sample.description
        print e
        return None,None
    l1,l2 = statisticsDiagramArray.plotSample( self, squeeze(s.data), reference=squeeze(r.data), **kwargs )
    return l1,l2

class normalizedStatisticsDiagram(object) :
  """Taylor and bias diagrams where the dimensional measures (std, crmse and bias) have been normalized."""
  def __init__(self, figsize=None, normalizeBias=False) :
    """Creates a new empty normalized diagram"""
    if figsize == None:
      figsize=(9,5)
    self.fig = plt.figure(figsize=figsize)
    self.fig.subplots_adjust(top=1.00,bottom=0.20) # for (9,5) fig
    self.taylorAx = None
    self.taylorParentAx = None
    self.biasAx = None
    self.normalizeBias = normalizeBias
    
    self.taylorDiags = {}
    self.biasDiags = {}
    self.stats = {}
    self.murphyLim = [0.0, 1.51]

  def addDataSet(self, tag, refsample, refLabel='Reference', **defArgs) :
    """Adds a new reference sample to the diagram, identified by a string tag.
    Use the same tag to plot samples against this reference."""
    
    self.stats[tag] = statistics.skillAssessmentStats(refsample, refLabel)
    refStd = self.stats[tag][(refLabel,'stddev')]
    
    self.taylorDiags[tag] = TaylorDiagram(refStd, refLabel, normalized=True, **defArgs)
    #self.biasDiags[tag] = BiasMurphyDiagram(refStd, refLabel, normalized=True,
                                      #normalizeBias=self.normalizeBias, **defArgs)
    self.biasDiags[tag] = BiasNMSEDiagram(refStd, refLabel, normalized=True,
                                      normalizeBias=self.normalizeBias, **defArgs)
    
    if self.taylorAx == None :
      # first dataSet, create axes
      ax, pax = self.taylorDiags[tag].setupAxes(self.fig, 121)
      self.taylorAx = ax
      self.taylorParentAx = pax
      self.biasAx,foo = self.biasDiags[tag].setupAxes(self.fig, 122)
    else :
      # set axes
      self.taylorDiags[tag].setAxes(self.taylorAx,self.taylorParentAx)
      self.biasDiags[tag].setAxes(self.biasAx)
    
  def plotSample(self, tag, sample, reference=None, **kwargs) :
    """Plots a new sample in the diagram.
    The corresponding reference is identified by the tag."""

    if 'label' not in kwargs.keys() :
      raise Exception( 'keyword argument label must be provided' )
    l = kwargs['label']
    stats = self.stats[tag].addSample(sample, l, reference=reference)
    bias = stats['bias']
    nmse = stats['nmse']
    if self.normalizeBias :
      # divide bias by std of the reference
      bias = bias/self.taylorDiags[tag].refStd
    l1 = self.taylorDiags[tag].plotSample(stats['corr'], stats['stddev'],
                                          **kwargs)
    l2 = self.biasDiags[tag].plotSample(bias, nmse, **kwargs)
    self.murphyLim[0] = 0.0 #min(self.murphyLim[0], nmse*1.05)
    self.murphyLim[1] = max(self.murphyLim[1], nmse*1.05)
    self.biasAx.set_ylim(self.murphyLim)
    yrange = self.murphyLim[1]-self.murphyLim[0]
    xTickMax = max(self.biasAx.get_xticks())
    self.biasAx.set_aspect(2*xTickMax/yrange)
    return l1,l2
  
  def getStatistics(self,tag) :
    """Returns statistics dictionary for the given tag.
    In the dictionary, keys are (label, statStr) tuples, where statStr is one of the following:
    'bias'   -- bias
    'ioa'    -- index of agreement
    'mur'    -- murphy skill score
    'stddev' -- standard deviation 
    'crmse'   -- centered root mean square error
    'corr'   -- correlation coefficient
    """
    return self.stats[tag].getStatistics()
    
  def getStatisticsArray(self, tag) :
    """Returns an array of all the statistics for the given tag,
    a list of all labels and a list of variable names.
    """
    return self.stats[tag].getStatisticsArray()
    
  def showLegend(self, **kwargs) :
    if not self.biasAx :
      return
    kw = dict()
    kw['prop'] = {'size':8}
    kw['numpoints'] = 1
    if not 'loc' in kwargs :
      kw['loc'] = 'upper center'
      kw['bbox_to_anchor'] = (-0.1, -0.15)
      kw['ncol'] = 4
    kw.update(kwargs)
    self.biasAx.legend(**kw)
    
  def addTitle(self, titleStr, **kwargs) :
    self.fig.suptitle(titleStr, **kwargs)

  def isEmpty( self ) :
    for td in self.taylorDiags :
      return False
    return True

class normalizedStatisticsDiagramDC(normalizedStatisticsDiagram) :
  """A wrapper class that uses dataContainer objects as inputs.
  Input data is automatically cast in same time units and
  interpolated to the reference time instances."""
  
  def __init__(self, figsize=None, normalizeBias=False) :
    self.dataContainers = {}
    normalizedStatisticsDiagram.__init__(self, figsize, normalizeBias)

  # TODO add a version that takes skillAssessmentStats object
  def addDataSet(self, tag, refsample, refLabel='Reference', **defArgs) :
    """Adds a new reference sample to the diagram, identified by a string tag.
    Use the same tag to plot samples against this reference."""
    statisticsDiagramDC.checkData( refsample )
    self.dataContainers[tag] = refsample
    normalizedStatisticsDiagram.addDataSet( self, tag, squeeze(refsample.data), refLabel, **defArgs )
  
  def plotSample(self, tag, sample, **kwargs) :
    """Plots a new sample in the diagram.
    The corresponding reference is identified by the tag."""
    statisticsDiagramDC.checkData( sample )
    if not 'label' in kwargs :
      kwargs['label'] = sample.description
    r = self.dataContainers[tag]
    s = sample
    if self.dataContainers[tag].time != sample.time :
      try :
        r,s = self.dataContainers[tag].alignTimes( sample )
      except Exception as e:
        print 'Time series could not be aligned, skipping:', kwargs['label'], sample.description
        print e
        return None,None
    l1,l2 = normalizedStatisticsDiagram.plotSample(self, tag, squeeze(s.data), reference=squeeze(r.data), **kwargs)
    return l1,l2
    
if __name__=='__main__':

  # TODO add options for grid spacing, rms curves etc
  
  # -- examples with concurrent data arrays --
  # generate data
  random.seed(int(3411))
  x = linspace(0,4*pi,100)
  ref = sin(x) # Data
  m1 = ref + 0.2*random.randn(len(x)) # Model 1
  m2 = 0.8*ref - 0.1*random.randn(len(x)) # Model 2
  m3 = sin(x-pi/10)+ 0.03*random.randn(len(x))# Model 3
  
  # -- statistics object --
  sas = statistics.skillAssessmentStats(ref)
  sas.addSample(m1,'m1')
  sas.addSample(m2,'m2')
  sas.addSample(m3,'m3')
  st = sas.getStatistics()
  refStd = st[('reference','stddev')]

  # -- taylor --
  fig = plt.figure(figsize=(5,4))
  dia = TaylorDiagram(refStd, unit='m')
  ax1 = dia.setupAxes(fig, 111)
  dia.plotSample(st[('m1','corr')],st[('m1','stddev')], label='m1')
  dia.plotSample(st[('m2','corr')],st[('m2','stddev')], label='m2')
  dia.plotSample(st[('m3','corr')],st[('m3','stddev')], label='m3')
  dia.addTitle('Taylor example')
  dia.showLegend(prop=dict(size='small'),loc='lower left',bbox_to_anchor = (0.82, 0.75))
  plt.show()
  
  # -- bias --
  fig = plt.figure(figsize=(5,4))
  dia = BiasDiagram(refStd, marker='o', unit='m')
  ax2 = dia.setupAxes(fig, 111)
  dia.plotSample(st[('m1','bias')],st[('m1','ioa')], label='m1', color='m')
  dia.plotSample(st[('m2','bias')],st[('m2','ioa')], label='m2', color='y')
  dia.plotSample(st[('m3','bias')],st[('m3','ioa')], label='m3', color=[0.1,0.9,0.4])
  dia.addTitle('Bias example')
  dia.showLegend(('obs','m1','m2','m3'),prop=dict(size='small'),loc='lower left')
  plt.show()

  # -- statisticsDiagram --
  dia = statisticsDiagram(refStd, marker='*', unit='m')
  dia.plotSample(st[('m1','corr')],st[('m1','stddev')],st[('m1','bias')],st[('m1','ioa')],
                 label='model 1')
  dia.plotSample(st[('m2','corr')],st[('m2','stddev')],st[('m2','bias')],st[('m2','ioa')],
                 label='model 2',color=[0.1,0.9,0.4],marker='^',markersize=10)
  dia.plotSample(st[('m3','corr')],st[('m3','stddev')],st[('m3','bias')],st[('m3','ioa')],
                 label='model 3',markersize=8)
  dia.showLegend()
  dia.addTitle('statisticsDiagram example',size=14)
  plt.show()
    
  dia.showLegend(prop=dict(size='small'),loc='lower left')
  dia.addTitle('normalizedStatisticsDiagram example',size=14)
  ## -- statisticsDiagram directly from skillAssessmentStats object --
  #dia = statisticsDiagram.fromStats(sas, marker='*', unit='m')
  #dia.showLegend(prop=dict(size='small'),loc='lower left')
  #dia.addTitle('statisticsDiagram.fromStats example',size=14)
  #plt.show()
  
  # -- statisticsDiagram with arrays --
  dia = statisticsDiagramArray(ref, marker='*', unit='m')
  dia.plotSample(m1,label='model Eins')
  dia.plotSample(m2,label='model Zwei',color=[0.1,0.9,0.4],marker='^',markersize=10)
  dia.plotSample(m3,label='model Drei',markersize=8)
  dia.showLegend()
  dia.addTitle('statisticsDiagramArray example',size=14)
  plt.show()
  s = dia.getStatistics() # a dictionary with all the statistics
  print s
  arr,labels,statNames = dia.getStatisticsArray() # array with all the stats
  print labels, statNames
  print arr

  # -- normalizedStatisticsDiagram --
  dia = normalizedStatisticsDiagram(figsize=(9,8))
  dia.addDataSet('one',ref, marker='o')
  dia.plotSample('one',m1, label='station1 m1')
  dia.addDataSet('two',ref, marker='v')
  dia.plotSample('two',m2, label='station2 m1')
  dia.plotSample('two',m3, label='station2 m2')
  dia.addTitle('normalizedStatisticsDiagram example',size=14)
  dia.showLegend()
  plt.show()
  s = dia.getStatistics('two')
  print s

  # -- examples with dataContainer objects --
  # generate data
  from crane.data import timeArray
  random.seed(int(3411))
  t0 = hstack( ( linspace(0,2,20), linspace(2.33,6,15) ) ) + 2193
  m0 = sin(t0)
  ta0 = timeArray.timeArray(t0,'corie')
  
  t1 = linspace(-10,10,100) + 2193
  m1 = 0.8*sin(t1) + 0.03*random.randn(len(t1))
  ta1 = timeArray.timeArray(t1,'corie')
  
  t2 = linspace(-9,17.7,65) + 2193
  m2 = 0.95*sin(t2) + 0.12*random.randn(len(t2))
  ta2 = timeArray.timeArray(t2,'corie')
  
  t3 = linspace(-9,12.2,100) + 2193
  m3 = sin(t3-0.12) - 0.05*random.randn(len(t3))
  ta3 = timeArray.timeArray(t3,'corie')
  
  d0 = dataContainer.fromTimeSeries( 'Observation', ta0, m0, ['elev'] )
  d1 = dataContainer.fromTimeSeries( 'model Eins', ta1, m1, ['elev'] )
  d2 = dataContainer.fromTimeSeries( 'model Zwei', ta2, m2, ['elev'] )
  d3 = dataContainer.fromTimeSeries( 'model Drei', ta3, m3, ['elev'] )
  
  # plot
  dia = statisticsDiagramDC(d0, 'Observation', unit='m')
  dia.plotSample(d1)
  dia.plotSample(d2)
  dia.plotSample(d3)
  dia.showLegend()
  dia.addTitle('statisticsDiagramDC example',size=14)
  plt.show()
  arr,labels,statNames = dia.getStatisticsArray() # array with all the stats
  print labels, statNames
  print arr

  dia = normalizedStatisticsDiagramDC()
  dia.addDataSet('foo', d0, marker='^')
  dia.plotSample('foo', d1, label='foo '+d1.description)
  dia.plotSample('foo', d2, label='foo '+d2.description)
  dia.addDataSet('bar', d0, marker='*')
  dia.plotSample('bar', d3, label='bar '+d3.description)
  dia.showLegend()
  dia.addTitle('normalizedStatisticsDiagramDC',size=14)
  plt.show()
  arr,labels,statNames = dia.getStatisticsArray('foo') # array with all the stats
  print labels, statNames
  print arr
