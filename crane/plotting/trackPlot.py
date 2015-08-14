#!/usr/bin/python

"""
A class for plotting time dependent tracks, where x = x(t) etc., e.g. AUV, glider and profiler data.

Tuomas Karna 2012-11-02
"""

from data.dataContainer import dataContainer
from data.timeArray import *

from plotting.plotBase import *

class trackTimeSeriesPlot(colorPlotBase) :
  """trackTimeSeriesPlot class"""
  def __init__(self, **defaultArgs) :
    # default plot options for all diagrams
    self.unit =   defaultArgs.pop('unit')
    self.xlabel = defaultArgs.pop('xlabel','Date (PST)')
    self.ylabel = defaultArgs.pop('ylabel','depth')
    self.clabel = defaultArgs.pop('clabel')
    self.clim =   defaultArgs.pop('clim',[])
    self.ylim =   defaultArgs.pop('ylim',[])
    self.xlim =   defaultArgs.pop('xlim',[])
    self.xIsTime =   defaultArgs.pop('xIsTime',True)
    self.yunit =  defaultArgs.pop('yunit','m')
    self.xunit =  defaultArgs.pop('xunit','')
    self.climIsLog = defaultArgs.pop('climIsLog',False)
    self.logScale = defaultArgs.pop('logScale',False)
    if self.logScale and self.unit :
      self.unit = 'log10('+self.unit+')'
    colorPlotBase.__init__(self,**defaultArgs)
    if self.xlim and self.xIsTime :
      self.xlim = [ convertEpochToPlotTime( datetimeToEpochTime( dt ) ) for dt in self.xlim ]

  def setAxes(self, ax) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    self.ax = ax
    ax.grid()
    yStr = self.ylabel
    if self.yunit :
      yStr += ' ['+self.yunit+']'
    ax.set_ylabel(yStr)
    xStr = self.xlabel
    if self.xunit :
      xStr += ' ['+self.xunit+']'
    ax.set_xlabel(xStr)

  def addSample(self, time, zCoord, variable, **kwargs) :
    """Add signal to the diagram. Time is assumed to be in epoch format."""
    if self.xIsTime :
      time = convertEpochToPlotTime( time )
    kw = dict(self.defaultArgs)
    kw.update(kwargs)
    
    xlim = kw.pop('xlim',None)
    
    data = variable if not self.logScale else np.log10( variable )

    if not 's' in kw :
      kw['s'] = 15
    if not 'edgecolors' in kw :
      kw['edgecolors'] = 'None'
    if self.clim and not 'vmin' in kw:
      cmin, cmax = self.clim
      if self.logScale and not self.climIsLog :
        cmin, cmax = np.log10(cmin), np.log10(cmax)
      kw['vmin'] = cmin
      kw['vmax'] = cmax
    self.cax = self.ax.scatter(time, zCoord, c=data, **kw)

    # tight y-axis
    pad = 0.02*(self.ax.dataLim.ymax - self.ax.dataLim.ymin)
    ylim = [ self.ax.dataLim.ymin-pad, self.ax.dataLim.ymax+pad]
    if self.ylim :
      ylim = self.ylim
    self.ax.set_ylim( ylim )
    self.updateXAxis( xlim )


class trackTimeSeriesPlotDC(trackTimeSeriesPlot) :
  def addSample( self, dc , **kwargs ) :
    # add sample from dataContainer
    if not dc.isTrack() :
      raise Exception( 'given dataContainer does not contain track information' )
    t = dc.time.asEpoch().array
    z = -dc.z.flatten()
    data = dc.data[:,0,:].flatten()
    trackTimeSeriesPlot.addSample( self, t, z, data, **kwargs )

class stackTrackPlot(stackPlotBase) :
  """A class for stacking multiple tracks in the same plot."""
  def addPlot(self, tag, **kwargs) :
    """Adds a new subplot to the diagram"""
    kw = self.defArgs
    kw.update(kwargs)
    plot = trackTimeSeriesPlot(**kw)
    stackPlotBase.addPlot(self,plot,tag)

  def addSample( self, tag, *args, **kwargs ) :
    if not tag in self.tags :
      raise Exception('No plot with tag: '+tag)
      #self.addPlot( trackTimeSeriesPlot(**self.defArgs), tag=tag)
    self.plots[tag].addSample( *args, **kwargs )
    self.plots[tag].showColorBar()
    xlim = [+1e30,-1e30]
    for k in self.plots :
      xlim[0] = min( self.plots[k].ax.dataLim.xmin, xlim[0] )
      xlim[1] = max( self.plots[k].ax.dataLim.xmax, xlim[1] )
    self.plots[tag].ax.set_xlim(xlim)

class stackTrackPlotDC(stackTrackPlot) :
  """stackTrackPlot class that uses dataContainer as an input"""
  def __init__(self, **defArgs) :
    stackTrackPlot.__init__(self,**defArgs)

  def addSample( self, tag, dc, **kwargs ) :
    """Add track from dataContainer to subplot identified with tag.
    If given tag does not exist, new subplot is appended.
    """

    t = dc.time.asEpoch().array
    z = -dc.z.flatten()
    data = dc.data.flatten()
    stackTrackPlot.addSample( self, tag, t,z,data, **kwargs )


if __name__ == '__main__' :
  #dataDir = '/home/tuomas/workspace/cmop/projects/AUV_comp/f22/f22/data/track/'
  #d0 = dataContainer.loadFromNetCDF(dataDir+'AUV71_salt_trck_2012-11-01_2012-11-01.nc')
  dataDir = '/home/tuomas/workspace/cmop/projects/AUV_comp/AUV/data/'
  d0 = dataContainer.loadFromNetCDF(dataDir+'AUV71_salt_auv_2012-11-01_2012-11-01.nc')

  t = d0.time.asEpoch().array
  z = -d0.z.flatten()
  data = d0.data.flatten()

  # low level example
  fig = plt.figure(figsize=(10,5))
  dia = trackTimeSeriesPlot(clabel='Salinity',unit='psu',clim=[0,32], ylim=[-20,0],ylabel='depth below surface')
  dia.setupAxes( fig )
  dia.addSample( t, z, data, s=30 )
  dia.showColorBar( )
  dia.addTitle( 'trackTimeSeriesPlot example' )
  plt.show()

  # example with dataContainer
  fig = plt.figure(figsize=(10,5))
  dia = trackTimeSeriesPlotDC(clabel='Salinity',unit='psu',clim=[0,32], ylim=[-20,0],ylabel='depth below surface')
  dia.setupAxes( fig )
  dia.addSample( d0 )
  dia.showColorBar( )
  dia.addTitle( 'trackTimeSeriesPlotDC example' )
  plt.show()

  # low level stack example
  dia = stackTrackPlot(clabel='Salinity',unit='psu',clim=[0,32], ylim=[-20,0],ylabel='depth below surface')
  dia.addSample( 'top',t, z, data, s=30 )
  dia.addTitle( 'stackTrackPlot example',tag='top' )
  dia.addSample( 'next',t, z, data, s=30 )
  dia.addTitle( 'bottom title',tag='next' )
  plt.show()

  # stack example with dataContainer
  dia = stackTrackPlotDC(clabel='Salinity',unit='psu',clim=[0,32], ylim=[-20,0],ylabel='depth below surface')
  dia.addSample( 'top',d0, s=30 )
  dia.addTitle( 'stackTrackPlotDC example',tag='top' )
  dia.addSample( 'next',d0, s=30 )
  dia.addTitle( 'bottom title',tag='next' )
  plt.show()
  