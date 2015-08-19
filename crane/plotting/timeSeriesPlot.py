#!/usr/bin/python

"""
A class for plotting time series data.
Uses stackplot routine for plotting.

Tuomas Karna 2012-08-20
"""
import sys

import numpy as np
import matplotlib.pyplot as plt

from crane.data import dataContainer
from crane.data import timeArray
from crane.plotting import plotBase
# for machines where cmop pylib is not present
from crane.plotting import stackplot as wp

#try:
  #import cmop.webproduct as wp
#except ImportError:
  #import stackplot as wp # for machines where cmop pylib is not present

class timeSeriesPlot2(plotBase.plotBase) :

  def __init__(self, **defaultArgs) :
    """Constructor. defaultArgs are the default plot options."""
    self.nSamples = 0
    self.ylabel=defaultArgs.pop('ylabel','')
    self.xlabel=defaultArgs.pop('xlabel','Date')
    self.yunit=defaultArgs.pop('unit','')
    self.xIsTime=True
    super(timeSeriesPlot2, self).__init__(**defaultArgs)

  def setAxes(self, ax) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    self.ax = ax
    ax.grid()
    yStr = self.ylabel
    if self.yunit :
      yStr += ' ['+self.yunit+']'
    ax.set_ylabel(yStr,multialignment='center')
    xStr = self.xlabel
    ax.set_xlabel(xStr,multialignment='center')

  def addShadedRange(self, start, end, **kwargs) :
    """Adds a shaded time range in the background.

    start,end can be a datetime object or float representing epoch time.
    kwargs are passed to matplotlib axvspan routine."""
    if isinstance( start, datetime.datetime ) :
      start = timeArray.datetimeToEpochTime(start)
    if isinstance( end, datetime.datetime ) :
      end = timeArray.datetimeToEpochTime(end)
    start = plotBase.convertEpochToPlotTime(start)
    end = plotBase.convertEpochToPlotTime(end)
    kwargs.setdefault( 'facecolor', [0.8,1.0,0.85] )
    kwargs.setdefault( 'edgecolor', 'none' )
    self.ax.axvspan(start, end, **kwargs)
    
  def addSample(self, label, t, y, **kwargs):
    """Add time series to the plot. kwargs are
    passed to the matplotlib plot command."""
    kw = dict(self.defaultArgs)
    kw.update(kwargs)
    
    if 'title' in kw :
      self.ax.set_title( kw.pop('title') )
    xlim = kw.pop('xlim',None)
    ylim = kw.pop('ylim',None)
    if 'color' not in kw :
      kw['color'] = self.colorSequence[ self.nSamples ]
    # detect gaps
    te = plotBase.convertPlotToEpochTime(t)
    ta = timeArray.timeArray(te,'epoch')
    gapDt = kw.pop('gapDt', None)
    gaps,ranges,_ = ta.detectGaps(dt=gapDt)
    for i in range(ranges.shape[0]) :
      lab = label if i==0 else '_none'
      self.ax.plot( t[ranges[i,0]:ranges[i,1]], y[ranges[i,0]:ranges[i,1]],
                    label=lab, **kw )
    self.nSamples += 1
    # set xlim to dates
    self.updateXAxis( xlim=xlim )
    if ylim :
      self.ax.set_ylim(ylim)

class timeSeriesPlotDC2(timeSeriesPlot2) :

  def addSample(self, dc, **kwargs) :
    t = plotBase.convertEpochToPlotTime( dc.time.asEpoch().array )
    y = np.squeeze( dc.data )
    label = kwargs.pop('label')
    timeSeriesPlot2.addSample(self, label, t, y, **kwargs)

class stackTimeSeriesPlot(plotBase.stackPlotBase) :
  """A class for stacking multiple time series axes in the same plot."""
  def addPlot(self, tag, **kwargs) :
    """Adds a new subplot to the diagram"""
    kw = self.defArgs
    kw.update(kwargs)
    plot = timeSeriesPlot2(**kw)
    #plot.updateXAxis( [ pb.convertEpochToPlotTime( 0.0 ), pb.convertEpochToPlotTime( 365*86400 )] )
    super(stackTimeSeriesPlot, self).addPlot(plot, tag)

  def addSample( self, tag, *args, **kwargs ) :
    if not tag in self.tags :
      raise Exception('No plot with tag: '+tag)
      #self.addPlot( trackTimeSeriesPlot(**self.defArgs), tag=tag)
    self.plots[tag].addSample( *args, **kwargs )
    #self.plots[tag].showColorBar()
    #xlim = [+1e30,-1e30]
    #for k in self.plots :
      #xlim[0] = min( self.plots[k].ax.dataLim.xmin, xlim[0] )
      #xlim[1] = max( self.plots[k].ax.dataLim.xmax, xlim[1] )
    #self.plots[tag].ax.set_xlim(xlim)

class stackTimeSeriesPlotDC(stackTimeSeriesPlot) :
  """stackTimeSeriesPlot class that uses dataContainer as an input"""
  def __init__(self, **defArgs) :
    super(stackTimeSeriesPlotDC, self).__init__(**defArgs)

  def addSample( self, tag, dc, **kwargs ) :
    """Add data from dataContainer to subplot identified with tag.
    If given tag does not exist, new subplot is appended.
    """
    t = plotBase.convertEpochToPlotTime( dc.time.asEpoch().array )
    data = dc.data.flatten()
    label = kwargs.pop('label',dc.getMetaData('tag'))
    super(stackTimeSeriesPlotDC, self).addSample(tag, label,t,data, **kwargs )

  
class timeSeriesPlot(object) :
  """A wrapper class to call stackplot. """
  def __init__(self, varLabel, tref=None, yref=None, refLabel='Reference', unit='', ylim=None, **kwargs) :
    """refsample is the reference (data) sample to be compared to."""
    self.varLabel = varLabel
    self.unit = unit
    self.tref = tref
    self.yref = yref
    self.data = []
    
    # default stackplot arguments
    self.defArgs = {}
    self.defArgs['x_is_time'] = True
    self.defArgs['timetick'] = 'new'
    self.defArgs['xlabel'] = 'Time'
    if tref != None :
      self.defArgs['xlim'] = [ tref.min(), tref.max() ]
    self.defArgs['marker'] = 'None'
    self.defArgs['linewidth'] = 1.0
    self.defArgs['linestyle'] = 'solid'
    self.defArgs['width'] = 7.5
    self.defArgs['handlegaps'] = True
    self.defArgs.update( kwargs )
    self.ylim = ylim

    # default color sequence
    self.colors = ['r','b','g','k','m','c']
    # add reference
    if tref != None and yref != None :
      timeSeriesPlot.addSample( self, refLabel, tref, yref )

  def addSample(self, label, t, y, **kwargs):
    """Add time series to the plot. kwargs are
    passed to the stackplot command."""
    if 'color' not in kwargs :
      kwargs['color'] = self.colors[ len(self.data) ]
    if not kwargs.has_key('ylim') and self.ylim :
      kwargs['ylim'] = self.ylim
    if not 'timestep' in kwargs : # handlegaps workaround for sat03
      kwargs['timestep'] = np.mean( np.diff( t ) )
    sample = wp.dataset( label, X=t, Y=y, Z=None, **kwargs ) 
    self.data.append( sample )

  def getVarLabelWithUnit(self) :
    varLabel = str(self.varLabel)
    if self.unit :
      varLabel += ' ['+self.unit+']'
    return varLabel
    
  def getDataTuple(self) :
    """Returns a (varLabel, dataList) tuple for stackplot routine"""
    return (self.getVarLabelWithUnit(), self.data)
    
  def makePlot( self, **kwargs ) :
    """Creates the plot. Should be called after all the data has been added.
    kwargs are passed to the stackplot command."""
    kw = dict(self.defArgs)
    kw.update(kwargs)
    wp.stackplot( [self.getDataTuple()], **kw )

class timeSeriesPlotDC(timeSeriesPlot) :
  """timeSeriesPlot that uses dataContainer objects as inputs"""
  def __init__(self, varLabel, ref=None, refLabel='Reference', unit='', ylim=None, **kwargs) :
    if 'xlim' in kwargs :
      kwargs['xlim'] = [ plotBase.convertEpochToPlotTime( timeArray.datetimeToEpochTime( dt ) ) for dt in kwargs['xlim'] ]
    if ref == None :
      tref = yref = None
    else :
      tref = plotBase.convertEpochToPlotTime( ref.time.asEpoch().array )
      yref = np.squeeze( ref.data )
    timeSeriesPlot.__init__(self, varLabel, tref, yref, refLabel, unit, ylim, **kwargs)

  def addSample(self, sampleData, **kwargs) :
    """Add time series to the plot."""
    timeSeriesComboPlotDC.checkData( sampleData )
    if not 'label' in kwargs :
      label = sampleData.description
    else :
      label = kwargs.pop('label')
    t = plotBase.convertEpochToPlotTime( sampleData.time.asEpoch().array )
    y = np.squeeze( sampleData.data )
    timeSeriesPlot.addSample(self,label, t, y, **kwargs)

class timeSeriesErrorPlot(timeSeriesPlot) :
  """Plots deviation from the reference signal versus time.
  All signals must be defined at the same time steps, so that error
  can be computer directly."""
  
  def addSample(self, label, t, error, **kwargs) :
    """Add time series to the plot. Time array must agree with the reference data set.
    kwargs are passed to the stackplot command."""
    if 'color' not in kwargs :
      kwargs['color'] = self.colors[ len(self.data) ]
    if not kwargs.has_key('ylim') and self.ylim :
      kwargs['ylim'] = self.ylim
    if not 'timestep' in kwargs : # handlegaps workaround for sat03
      kwargs['timestep'] = np.mean( np.diff( t ) )
    sample = wp.dataset( label, X=t, Y=error, Z=None, **kwargs ) 
    self.data.append( sample )
  
  def getDataTuple(self) :
    """Returns a (varLabel, dataList) tuple for stackplot routine"""
    return (self.getVarLabelWithUnit(), self.data[1:])

class timeSeriesComboPlot(object) :
  """Plots time series, and deviation from the reference signal versus time.
  All the signals must be defined on the same time instance, so that error
  can be computer directly."""
  
  def __init__(self, varLabel, tref, yref, refLabel='Reference', errLabel='Error', unit='', ylim=None, err_ylim=None, **kwargs) :
    figsize = kwargs.pop('figsize',(12,5))
    fig = plt.figure(figsize=figsize)
    width = figsize[0]
    self.tsPlot = timeSeriesPlot(varLabel, tref, yref, refLabel, unit, width=width, ylim=ylim, **kwargs)
    self.errPlot = timeSeriesErrorPlot(errLabel, tref, yref, refLabel, unit, width=width, ylim=err_ylim, **kwargs)
    
  def addSample(self, label, t, y, error, **kwargs) :
    """Add time series to the plot. kwargs are
    passed to the stackplot command."""
    self.tsPlot.addSample(label, t, y, **kwargs)
    self.errPlot.addSample(label, t, error, **kwargs)
  
  def makePlot( self, **kwargs ) :
    """Creates the plot. Should be called after all the data has been added.
    kwargs are passed to the stackplot command."""
    tsTuple = self.tsPlot.getDataTuple()
    errTuple = self.errPlot.getDataTuple()
    # if ylim not set, force a minimum ylim
    errData = errTuple[1]
    if errData :
      errKw = errData[0].kwargs
      if not errKw.has_key('ylim') :
        dataLim = np.array([1e40,-1e40])
        for d in errData :
          dataLim[0] = min( d.Y.min(), dataLim[0] )
          dataLim[1] = max( d.Y.max(), dataLim[1] )
        dataRange = dataLim[1]-dataLim[0]
        ylimThreshold = .2
        if dataRange < ylimThreshold :
          dataLim = [ dataLim.mean()-ylimThreshold/2,dataLim.mean()+ylimThreshold/2 ]
          errKw['ylim'] = dataLim
    data = [ tsTuple, errTuple ]
    self.tsPlot.defArgs.update(kwargs)
    wp.stackplot( data, **self.tsPlot.defArgs )



class timeSeriesComboPlotDC(timeSeriesComboPlot) :
  """timeSeriesComboPlot that uses dataContainer objects as inputs"""
  def __init__(self, varLabel, refData, refLabel=None, errLabel='Error', unit='', ylim=None, err_ylim=None, **kwargs) :
    self.refData = refData
    tref = plotBase.convertEpochToPlotTime( refData.time.asEpoch().array )
    if 'xlim' in kwargs :
      kwargs['xlim'] = [ plotBase.convertEpochToPlotTime( timeArray.datetimeToEpochTime( dt ) ) for dt in kwargs['xlim'] ]
    yref = np.squeeze( refData.data )
    if not refLabel :
      refLabel = refData.description
    timeSeriesComboPlot.__init__(self,varLabel, tref, yref, refLabel, errLabel, unit, ylim=ylim, err_ylim=err_ylim, **kwargs)
    
  @staticmethod
  def checkData( sample ) :
    if not isinstance( sample, dataContainer.dataContainer ) :
      raise Exception( 'sample must be a dataContainer object' )
    if sample.data.shape[0] > 1 :
      raise Exception( 'spatially varying data is not supported' )
    if sample.data.shape[1] > 1 :
      raise Exception( 'multiple fields are not supported' )
    return True

  def addSample(self, sampleData, **kwargs) :
    """Adds a sample to the plot. For the error plot the data is first interpolated to
    the reference time steps. kwargs are passed to the stackplot command."""
    timeSeriesComboPlotDC.checkData( sampleData )
    if not 'label' in kwargs :
      label = sampleData.description
    else :
      label = kwargs.pop('label')
    t = plotBase.convertEpochToPlotTime( sampleData.time.asEpoch().array )
    y = np.squeeze( sampleData.data )
    self.tsPlot.addSample(label, t, y, **kwargs)
    try :
      err = self.refData.computeError( sampleData )
      t = plotBase.convertEpochToPlotTime( err.time.asEpoch().array )
      self.errPlot.addSample(label, t, np.squeeze( err.data ), **kwargs)
    except Exception as e :
      print 'Error could not be computed, skipping:', label, sampleData.description
      print e

class timeSeriesStackPlot(object) :
  def __init__( self, **kwargs) :
    figsize = kwargs.get('figsize', (12,5) )
    fig = plt.figure(figsize=figsize)
    if not 'width' in kwargs :
      kwargs['width'] = figsize[0]
    if 'xlim' in kwargs :
      kwargs['xlim'] = [ plotBase.convertEpochToPlotTime( timeArray.datetimeToEpochTime( dt ) ) for dt in kwargs['xlim'] ]
    self.tags = [] # one for each subplot
    self.data = dict() # list of time series for each subplot
    self.varLabel = dict() # one for each subplot
    self.unit = dict()
    # default stackplot options
    self.defArgs = {}
    self.defArgs['x_is_time'] = True
    self.defArgs['timetick'] = 'new'
    self.defArgs['xlabel'] = 'Time'
    self.defArgs['marker'] = 'None'
    self.defArgs['linewidth'] = 1.0
    self.defArgs['linestyle'] = 'solid'
    self.defArgs['width'] = 17.5
    self.defArgs['handlegaps'] = True
    self.defArgs.update( kwargs )

  def getVarLabelWithUnit(self, tag) :
    varLabel = str(self.varLabel[tag])
    if self.unit[tag] :
      varLabel += ' ['+self.unit[tag]+']'
    return varLabel

  def addSubplot( self, tag, varLabel, unit ) :
    if not tag in self.tags :
      self.tags.append( tag )
      self.data[tag] = []
      self.unit[tag] = unit
      self.varLabel[tag] = varLabel

  def addSample( self, tag, label, t,y, **kwargs ) :
    if not tag in self.tags :
      raise Exception( 'Given tag cannot be found')
    if not 'timestep' in kwargs : # handlegaps workaround for sat03
      kwargs['timestep'] = np.mean( np.diff( t ) )
    sample = wp.dataset(label, X=t, Y=y, Z=None, **kwargs)
    self.data[tag].append( sample )

  def makePlot( self, **kwargs ) : # TODO remove kwargs from here?
    data = []
    for tag in self.tags :
      data.append( ( self.getVarLabelWithUnit(tag), self.data[tag] ) )
    self.defArgs.update(kwargs)
    wp.stackplot( data, **self.defArgs )

  def isEmpty( self ) :
    for tag in self.data :
      for ts in self.data[tag] :
        return False
    return True

class timeSeriesStackPlotDC(timeSeriesStackPlot) :
  def addSample( self, tag, sampleData, **kwargs ) :
    timeSeriesComboPlotDC.checkData( sampleData )
    if not 'label' in kwargs :
      label = sampleData.description
    else :
      label = kwargs.pop('label')
    if 'xlim' in kwargs :
      kwargs['xlim'] = [ plotBase.convertEpochToPlotTime( timeArray.datetimeToEpochTime( dt ) ) for dt in kwargs['xlim'] ]
    t = plotBase.convertEpochToPlotTime( sampleData.time.asEpoch().array )
    y = np.squeeze( sampleData.data )
    timeSeriesStackPlot.addSample(self,tag,label, t, y, **kwargs)

if __name__=='__main__':

  import datetime

  ### examples with numpy array inputs
  # generate data
  startTime = datetime.datetime(2010,1,12,0,0,0)
  endTime = datetime.datetime(2010,2,13,3,30,0)
  dt = 900.0
  sta = timeArray.generateSimulationTimeArray(startTime,endTime,dt).asEpoch()
  t = sta.array
  
  T = 44714
  ref = np.sin(t/T) + 0.95*np.sin(0.95*t/T)
  m1 = 0.90*np.sin(t/T) + 0.7*np.sin(0.85*t/T+0.3) - 0.08
  m2 = 0.80*np.sin(t/T+0.8) + 0.9*np.sin(0.95*t/T+0.5) - 0.12
  m3 = 0.78*np.sin(t/T) + 0.82*np.sin(0.90*t/T) + 0.03
  
  N = len(t)
  ref[N/2:N/2+N/10] = np.nan # missing values
  tn = plotBase.convertEpochToPlotTime( t )
 
  # plot
  fig = plt.figure(figsize=(10,5))
  kwargs = {}
  kwargs['title'] = 'timeSeriesPlot example'
  kwargs['width'] = 15
  dia = timeSeriesPlot('Elevation', tn, ref, refLabel='observation', unit='m', ylim=[-2.2,2.5])
  dia.addSample('m1', tn, m1, linestyle='-', color='b')
  dia.addSample('m2', tn, m2, linestyle='-', color=[0.1,0.5,0.1])
  dia.addSample('m3', tn, m3, linestyle='dashed', color='k')
  dia.makePlot(**kwargs)
  plt.show()
  
  fig = plt.figure(figsize=(10,5))
  kwargs['title'] = 'timeSeriesErrorPlot example'
  kwargs['width'] = 15
  dia = timeSeriesErrorPlot('Error', tn, ref, refLabel='observation', unit='m',linewidth=1)
  dia.addSample('m1', tn, m1-ref, linestyle='-', color='b')
  dia.addSample('m2', tn, m2-ref, linestyle='-', color=[0.1,0.5,0.1])
  dia.addSample('m3', tn, m3-ref, linewidth=4, linestyle='dashed', color='k')
  dia.makePlot(**kwargs)
  plt.show()
  
  kwargs['title'] = 'timeSeriesComboPlot example'
  dia = timeSeriesComboPlot('Elevation', tn, ref, refLabel='observation', unit='m', figsize=(15,5), ylim=[-2,2.2], err_ylim=[-1.5,1.5])
  #dia = timeSeriesComboPlot('Elevation', tn, ref, refLabel='observation', unit='m', linewidth=3)
  dia.addSample('m1', tn, m1, m1-ref, linestyle='-')
  dia.addSample('m2', tn, m2, m2-ref, linestyle='-')
  dia.addSample('m3', tn, m3, m3-ref, linestyle='dashed')
  dia.makePlot(**kwargs)
  plt.show()

  dia = timeSeriesStackPlot(title='timeSeriesStackPlot example',figsize=(15,5))
  dia.addSubplot( 'top', 'Elevation', 'm')
  dia.addSubplot( 'mid', 'Elevation', 'm')
  dia.addSubplot( 'bot', 'Elevation', 'm')
  dia.addSample('top', 'm1', tn, m1)
  dia.addSample('top', 'm2', tn, m2,color='m')
  dia.addSample('mid', 'm1', tn, m1)
  dia.addSample('mid', 'm3', tn, m3)
  dia.addSample('bot', 'm1', tn, m1)
  dia.makePlot()
  plt.show()

  # example with dataContainer objects
  # generate data
  np.random.seed(int(3411))
  startTime = datetime.datetime(2010,1,12,0,0,0)
  startCorie = timeArray.datetimeToCorieTime(startTime)
  t0 = np.hstack( ( np.linspace(0,12,20), np.linspace(15.33,30,15) ) ) + startCorie
  m0 = np.sin(t0)
  ta0 = timeArray.timeArray(t0,'corie')
  
  t1 = np.hstack( ( np.linspace(-10,17,50), np.linspace(22,34,50) ) )  + startCorie
  m1 = 0.8*np.sin(t1) + 0.03*np.random.randn(len(t1))
  ta1 = timeArray.timeArray(t1,'corie')
  
  t2 = np.linspace(-9,31.7,65) + startCorie
  m2 = 0.95*np.sin(t2) + 0.12*np.random.randn(len(t2))
  ta2 = timeArray.timeArray(t2,'corie')
  
  t3 = np.linspace(-9,32.2,100) + startCorie
  m3 = np.sin(t3-0.12) - 0.05*np.random.randn(len(t3))
  ta3 = timeArray.timeArray(t3,'corie')
  
  d0 = dataContainer.dataContainer.fromTimeSeries( 'Observation', ta0, m0, ['elev'] )
  d1 = dataContainer.dataContainer.fromTimeSeries( 'model Eins', ta1, m1, ['elev'] )
  d2 = dataContainer.dataContainer.fromTimeSeries( 'model Zwei', ta2, m2, ['elev'] )
  d3 = dataContainer.dataContainer.fromTimeSeries( 'model Drei', ta3, m3, ['elev'] )

  kwargs = {}
  kwargs['title'] = 'timeSeriesPlotDC example'
  dia = timeSeriesPlotDC( 'Elevation', d0, refLabel='observation', unit='m')
  dia.addSample(d1, label='m1', linestyle='-', color='b')
  dia.addSample(d2, label='m2', linestyle='-', color=[0.1,0.5,0.1])
  dia.addSample(d3, label='m3', linestyle='dashed', color='k')
  dia.makePlot(**kwargs)
  plt.show()

  # plot
  kwargs = {}
  kwargs['title'] = 'timeSeriesComboPlotDC example'
  xlim = [datetime.datetime(2010,1,5,0,0,0), datetime.datetime(2010,2,10,0,0,0)]
  dia = timeSeriesComboPlotDC('Elevation', d0 ,refLabel='custom ref', unit='m', xlim=xlim,
                              ylim=[-1.2,1.1], err_ylim=[-0.4,0.5], linestyle='dashed',linewidth=1.0)
  dia.addSample(d1, color='b')
  dia.addSample(d2, color=[0.1,0.5,0.1], label='custom label')
  dia.addSample(d3, color='k', linestyle='dashed')
  dia.makePlot(**kwargs)
  plt.show()

  dia = timeSeriesStackPlotDC(title='timeSeriesStackPlotDC example',ylim=[-4.5,5],xlabel='custom x label')
  dia.addSubplot( 'top', 'tansy'+'\n''Elevation', 'm')
  dia.addSubplot( 'mid', 'Elevation', 'm')
  dia.addSubplot( 'bot', 'Elevation', 'm')
  dia.addSample('top', d0, color='r')
  dia.addSample('top', d1, color='b')
  dia.addSample('mid', d0, color='r')
  dia.addSample('mid', d2, color='g')
  dia.addSample('bot', d0, color='r')
  dia.addSample('bot', d3, color='m')
  dia.makePlot()
  plt.show()
