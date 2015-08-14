#!/usr/bin/python
"""
A class for plotting time series of vertical profiles.

Tuomas Karna 2012-11-02
"""

import numpy as np
from data.timeArray import *
from plotting.plotBase import *

class verticalProfilePlot(plotBase) :
  """Class that plots vertical profiles in (variable,z) line graphs."""
  def __init__(self, **defaultArgs) :
    """Constructor. defaultArgs are the default plot options.
    
    Parameters
    ----------
    ylabel : str, optional
        Label for vertical axis (default 'Depth')
    yunit : str, optional
        Unit for vertical axis (default 'm')
    xlabel : str, optional
        Label for horizontal axis (default '')
    xunit : str, optional
        Unit for horizontal axis (default '')
    invert_yaxis : bool, optional
        If set inverts vertical axis so y decreases upwards (default False)
    
    """
    self.nSamples = 0
    self.ylabel=defaultArgs.pop('ylabel','Depth')
    self.yunit=defaultArgs.pop('yunit','m')
    self.xlabel=defaultArgs.pop('xlabel','')
    self.xunit=defaultArgs.pop('xunit','')
    self.invert_yaxis=defaultArgs.pop('invert_yaxis',False)
    plotBase.__init__(self, **defaultArgs)

  def setAxes(self, ax) :
    """Set axes for the diagram where data will be plotted.
    
    Parameters
    ----------
    ax : matplotlib.axes.AxesSubplot
    
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
    loc = matplotlib.ticker.MaxNLocator(nbins=6)
    ax.xaxis.set_major_locator(loc)
    if self.invert_yaxis :
      self.ax.invert_yaxis()

  def addSample(self, z, x, label, **kwargs):
    """Add profile to the plot.
    
    Parameters
    ----------
    z, x : numpy.ndarray
            Singleton arrays for the (x,z) plot. Must be of equal length.
    label : str
            Label identifying the data set. Will be printed in legend.
    Additional kwargs are passed to pyplot.plot command.

    """
    kw = dict(self.defaultArgs)
    kw.update(kwargs)

    if 'title' in kw :
      self.ax.set_title( kw.pop('title') )
    xlim = kw.pop('xlim',None)
    ylim = kw.pop('ylim',None)
    zerolinecolor = kw.pop('zerolinecolor','gray')
    if 'color' not in kw :
      kw['color'] = self.colorSequence[ self.nSamples ]
    lab = label
    self.ax.plot( x, z, label=lab, **kw )
    self.nSamples += 1
    if xlim :
      self.ax.set_xlim(xlim)
    if ylim :
      self.ax.set_ylim(ylim)
    self.ax.axvspan(0.0, 0.0, facecolor='none', edgecolor=zerolinecolor, lw=0.5)

class verticalProfilePlotDC(verticalProfilePlot) :
  """Class that uses dataContainer as inputs"""
  def addSample(self, dc, label, iTime=0, iField=0, **kw) :
    """Add profile to the plot.

    Parameters
    ----------
    dc : dataContainer
        vertical profile data. dc.z must have shape (nZ,) or (nZ,nTime)
        dc.data must have shape (nZ,nFields,nTime)
    label : str
         Label identifying the data set. Will be printed in legend.
    iTime : int, optional
         Specifies which time instance to plot (default 0)
    iField : int, optional
         Specifies which field to plot (default 0)

    Additional kwargs are passed to pyplot.plot command.

    """
    z = dc.z[:,iTime] if dc.zDependsOnTime else dc.z
    x = dc.data[:,iField,iTime]
    verticalProfilePlot.addSample(self,z,x,label,**kw)

class profileTimeSeries(colorPlotBase) :
  """Class that plots time-dependent vertical profiles in (t,z) plot,
  colored by variable."""
  def __init__(self, **defaultArgs) :
    # default plot options for all diagrams
    self.unit =   defaultArgs.pop('unit')
    self.xlabel = defaultArgs.pop('xlabel','Date')
    self.ylabel = defaultArgs.pop('ylabel','depth')
    self.clabel = defaultArgs.pop('clabel')
    self.clim =   defaultArgs.pop('clim',[])
    self.climIsLog = defaultArgs.pop('climIsLog',False)
    self.ylim =   defaultArgs.pop('ylim',[])
    self.xlim =   defaultArgs.pop('xlim',[])
    self.xIsTime =   defaultArgs.pop('xIsTime',True)
    self.yunit =  defaultArgs.pop('yunit','m')
    self.logScale = defaultArgs.pop('logScale',False)
    self.invert_yaxis=defaultArgs.pop('invert_yaxis',False)
    defaultArgs.setdefault('plotType','contourf')
    defaultArgs.setdefault('extend','both')
    defaultArgs.setdefault('N',20)
    if self.logScale and self.unit :
      self.unit = r'$\log10($'+self.unit+'$)$'
    colorPlotBase.__init__(self,**defaultArgs)
    if self.xlim and self.xIsTime :
      self.xlim = [ convertEpochToPlotTime( datetimeToEpochTime( dt ) )
                 for dt in self.xlim ]

  def setAxes(self, ax) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    self.ax = ax
    ax.grid()
    yStr = self.ylabel
    if self.yunit :
      yStr += ' ['+self.yunit+']'
    ax.set_ylabel(yStr)
    ax.set_xlabel(self.xlabel)
    if self.invert_yaxis :
      self.ax.invert_yaxis()

  def addSample(self, time, zCoord, variable, **kwargs) :
    """Add signal to the diagram. Time is assumed to be in epoch format."""
    t = convertEpochToPlotTime( time )
    z = zCoord.copy()
    if len(t.shape) == 1 and z.shape[1] == len(t) :
      # z time dependent, copy singleton time to correct shape
      t = np.tile( t, (z.shape[0],1) )
    else :
      raise Exception('arrays of given shape are not supported',
                    t.shape,z.shape )
    kw = dict(self.defaultArgs)
    kw.update(kwargs)
    plotType = kw.pop('plotType')

    data = variable.copy() if not self.logScale else np.log10( variable )

    N = kw.pop('N')
    if not 'levels' in kw :
      if self.clim :
        cmin, cmax = self.clim
        if self.logScale and not self.climIsLog :
          cmin, cmax = np.log10(cmin), np.log10(cmax)
        kw['levels'] = np.linspace( cmin, cmax, N )
        # saturate limits
        #data[ data < cmin] = cmin
        #data[ data > cmax] = cmax
    if plotType == 'contour' :
      self.cax = self.ax.contour(t,z,data,**kw)
    else :
      self.cax = self.ax.contourf(t,z,data,**kw)

    # tight y-axis
    ylim = [ self.ax.dataLim.ymin, self.ax.dataLim.ymax]
    if self.ylim :
      ylim = self.ylim
    self.ax.set_ylim( ylim )
    if self.invert_yaxis :
      self.ax.invert_yaxis()
    self.updateXAxis()

  def addOverlay(self, time, zCoord, variable, **kwargs) :
    """Adds time series on top of contour image. Time is assumed to be in epoch format."""
    t = convertEpochToPlotTime( time )
    z = zCoord.copy()
    if len(z) == 1 :
      z = np.ones_like(t)*z

    kw = dict(kwargs)

    data = variable.copy() if not self.logScale else np.log10( variable )

    if self.clim and not 'vmin' in kw:
      cmin, cmax = self.clim
      if self.logScale and not self.climIsLog :
       cmin, cmax = np.log10(cmin), np.log10(cmax)
      kw['vmin'] = cmin
      kw['vmax'] = cmax
    kw.setdefault('s',12)
    kw.setdefault('edgecolors','none')
    kw['s'] = 3*kw['s']
    kw.setdefault('marker','s')
    self.ax.scatter(t, z, c='w', zorder=10, **kw)
    kw['s'] = kw['s']/3
    self.ax.scatter(t, z, c=data, zorder=11, **kw)

  def updateXAxis( self, xlim=None, numticks=12 ) :
    # dateticks
    if not xlim :
      if self.xlim :
        xlim = self.xlim
      else :
        xlim = [ self.ax.dataLim.xmin, self.ax.dataLim.xmax ]
    self.ax.set_xlim( xlim )
    if self.xIsTime :
      loc, fmt = date_ticker_factory( xlim[1] - xlim[0], numticks=numticks)
      self.ax.xaxis_date()
      self.ax.xaxis.set_major_formatter(fmt)
      self.ax.xaxis.set_major_locator(loc)

class profileTimeSeriesDC(profileTimeSeries) :
  def addSample( self, dc , iField=0, **kwargs ) :
    """Add profile time series"""
    gaps, ranges, time = dc.detectGaps(gapFactor=10)
    for r in ranges :
      # plot each continuous data range separately
      t = time[r[0]:r[1]+1]
      zCoord = dc.z
      var = np.squeeze( dc.data[:,iField,r[0]:r[1]+1] )
      profileTimeSeries.addSample( self, t, zCoord, var, **kwargs )

  def addOverlay( self, dc , iField=0, **kwargs ) :
    """Overlay time series on plot"""
    gaps, ranges, time = dc.detectGaps(gapFactor=10)
    for r in ranges :
      # plot each continuous data range separately
      t = time[r[0]:r[1]+1]
      zCoord = dc.z
      var = np.squeeze( dc.data[:,iField,r[0]:r[1]+1] )
      profileTimeSeries.addOverlay( self, t, zCoord, var, **kwargs )

class stackProfileTimeSeries(stackPlotBase) :
  """A class for stacking multiple profiles in the same plot."""
  def addPlot(self, tag, **kwargs) :
    """Adds a new subplot to the diagram"""
    kw = dict(self.defArgs)
    kw.update(kwargs)
    plot = profileTimeSeries(**kw)
    stackPlotBase.addPlot(self,plot,tag)

  def addSample( self, tag, *args, **kwargs ) :
    if not tag in self.tags :
      self.addPlot( tag, **kwargs )
    self.plots[tag].addSample( *args, **kwargs )

  def addOverlay( self, tag, *args, **kwargs ) :
    if not tag in self.tags :
      self.addPlot( tag, **kwargs )
    self.plots[tag].addOverlay( *args, **kwargs )

class stackProfileTimeSeriesDC(stackProfileTimeSeries) :
  """stackProfileTimeSeries class that uses dataContainer as an input"""

  def addSample( self, tag, dc, **kwargs ) :
    """Add transect from dataContainer wih appropriate data restructuring.
       Args:
       tag       -- (string) tag for identifiying the subplot
       dc        -- (dataContainer) data
       kwargs    -- (**) passed to profileTimeSeries.addSample
    """
    gaps, ranges, time = dc.detectGaps(gapFactor=10)
    for r in ranges :
      # plot each continuous data range separately
      t = time[r[0]:r[1]+1]
      if len(t) < 2 :
        continue
      if dc.zDependsOnTime :
        zCoord = dc.z[:,r[0]:r[1]+1]
      else :
        zCoord = dc.z
      var = np.squeeze( dc.data[:,0,r[0]:r[1]+1] )
      stackProfileTimeSeries.addSample( self, tag, t, zCoord, var, **kwargs )

  def addOverlay( self, tag, dc, **kwargs ) :
    """Add transect from dataContainer wih appropriate data restructuring.
       Args:
       tag       -- (string) tag for identifiying the subplot
       dc        -- (dataContainer) data
       kwargs    -- (**) passed to profileTimeSeries.addSample
    """
    gaps, ranges, time = dc.detectGaps(gapFactor=10)
    for r in ranges :
      # plot each continuous data range separately
      t = time[r[0]:r[1]+1]
      if len(t) < 2 :
        continue
      if dc.zDependsOnTime :
        zCoord = dc.z[:,r[0]:r[1]+1]
      else :
        zCoord = dc.z
      var = np.squeeze( dc.data[:,0,r[0]:r[1]+1] )
      stackProfileTimeSeries.addOverlay( self, tag, t, zCoord, var, **kwargs )
