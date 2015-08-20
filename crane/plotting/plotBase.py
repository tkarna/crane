#!/usr/bin/python

"""
Abstract base class(es) for plotting routines, to enforce standard interface,
and plotting helper routines.

Should only depend on matplotlib.

Tuomas Karna 2012-11-29
"""

import os
import numpy as np
import datetime
import matplotlib
if 'TACC_SYSTEM' in os.environ and os.environ['TACC_SYSTEM'] == 'stampede':
    # change backend if no display present (e.g. on supercomputers)
    matplotlib.use("Agg", warn=False)
elif 'stccmop' in os.environ.get('HOSTNAME', []):
    matplotlib.use("Agg", warn=False)
elif 'edison' in os.environ.get('HOSTNAME', []):
    matplotlib.use("Agg", warn=False)
elif 'hopper' in os.environ.get('HOSTNAME', []):
    matplotlib.use("Agg", warn=False)

import matplotlib.pyplot as plt
import matplotlib.dates
from mpl_toolkits.axes_grid1 import make_axes_locatable

from crane.data import timeArray

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

# sequence of colors for line plots etc
defaultColorSequence = ['r','b','g','k','#D100D1','#00C2C2','Chocolate','MediumPurple','GoldenRod','#348ABD', '#7A68A6', '#A60628', '#467821', '#CF4457', '#188487', '#E24A33']

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------


def updateXAxis(ax, xIsTime=False, xlim=None, minticks=3, maxticks=12, prune=None):
    """Updates x axis ticks.

    xlim sets the range for the axis, if xlim==None the range is not changed.
    If xlim='tight' x limits are set to data range.
    prune='lower'|'upper'|'both'|None can be used to omit low/high tick labels
    """
    if xlim == 'tight' :
      xlim = [ ax.dataLim.xmin, ax.dataLim.xmax ]
    if xlim == None :
      xlim = ax.get_xlim()
    if xIsTime :
      if xlim[0]==xlim[1] :
        offset = 1./24 # days
        xlim[0] -= offset/2
        xlim[1] += offset/2
      ax.set_xlim(xlim)

      hourfmt = '%H:%M'
      monthfmt = '%b'
      timerange = [ax.dataLim.xmin, ax.dataLim.xmax]
      span = timerange[1]-timerange[0]
      if span <= 1.0:
        # all data from same day
        meanTime = timeArray.epochToDatetime(convertPlotToEpochTime(np.mean(timerange)))
        ax.set_xlabel(meanTime.strftime('%Y-%m-%d'))
      elif span < 25.0:
        # all data from same month
        meanTime = timeArray.epochToDatetime(convertPlotToEpochTime(np.mean(timerange)))
        ax.set_xlabel(meanTime.strftime('%Y'))
        hourfmt = '%H:%M\n%b %d'
      elif span < 366.1:
        # all data from same year
        meanTime = timeArray.epochToDatetime(convertPlotToEpochTime(np.mean(timerange)))
        ax.set_xlabel(meanTime.strftime('%Y'))
      else:
        monthfmt = '%b\n%Y'
        ax.set_xlabel('')

      loc = matplotlib.dates.AutoDateLocator(minticks=minticks,
                                             maxticks=maxticks,
                                             interval_multiples=True)
      fmt = matplotlib.dates.AutoDateFormatter(loc)
      # tick format for tick scale (difference between major ticks)
      fmt.scaled = {365.0  : '%Y',
                    30.    : monthfmt,
                    1.0    : '%b %d',
                    0.5    : '%b-%d %H:%M',
                    1./24. : hourfmt,
                    1. / (24.*60.): '%H:%M:%S',
                    }
      ax.xaxis_date()
      ax.xaxis.set_major_formatter(fmt)
      ax.xaxis.set_major_locator(loc)
      ax.set_xlim(xlim)
    else :
      if xlim :
        ax.set_xlim(xlim)
      loc = matplotlib.ticker.MaxNLocator(nbins=maxticks, prune=prune)
      ax.xaxis.set_major_locator(loc)


def updateYAxis(ax, yIsTime=False, ylim=None, minticks=3, maxticks=12, prune=None):
    """Updates y axis ticks.

    ylim sets the range for the axis, if ylim==None the range is not changed.
    If ylim='tight' y limits are set to data range.
    prune='lower'|'upper'|'both'|None can be used to omit low/high tick labels
    """
    if ylim == 'tight' :
      ylim = [ ax.dataLim.ymin, ax.dataLim.ymax ]
    if ylim == None :
      ylim = ax.get_ylim()
    if yIsTime :
      if ylim[0]==ylim[1] :
        offset = 1./24 # days
        ylim[0] -= offset/2
        ylim[1] += offset/2
      ax.set_ylim(ylim)

      hourfmt = '%H:%M'
      monthfmt = '%b'
      timerange = [ax.dataLim.ymin, ax.dataLim.ymax]
      span = timerange[1]-timerange[0]
      if span <= 1.0:
        # all data from same day
        meanTime = timeArray.epochToDatetime(convertPlotToEpochTime(np.mean(timerange)))
        ax.set_ylabel(meanTime.strftime('%Y-%m-%d'))
      elif span < 25.0:
        # all data from same month
        meanTime = timeArray.epochToDatetime(convertPlotToEpochTime(np.mean(timerange)))
        ax.set_ylabel(meanTime.strftime('%Y'))
        hourfmt = '%H:%M\n%b %d'
      elif span < 366.1:
        # all data from same year
        meanTime = timeArray.epochToDatetime(convertPlotToEpochTime(np.mean(timerange)))
        ax.set_ylabel(meanTime.strftime('%Y'))
      else:
        monthfmt = '%b\n%Y'
        ax.set_ylabel('')

      loc = matplotlib.dates.AutoDateLocator(minticks=minticks,
                                             maxticks=maxticks,
                                             interval_multiples=True)
      fmt = matplotlib.dates.AutoDateFormatter(loc)
      # tick format for tick scale (difference between major ticks)
      fmt.scaled = {365.0  : '%Y',
                    30.    : monthfmt,
                    1.0    : '%b %d',
                    0.5    : '%b-%d %H:%M',
                    1./24. : hourfmt,
                    1. / (24.*60.): '%H:%M:%S',
                    }
      ax.yaxis_date()
      ax.yaxis.set_major_formatter(fmt)
      ax.yaxis.set_major_locator(loc)
    else :
      if ylim :
        ax.set_ylim(ylim)
      loc = matplotlib.ticker.MaxNLocator(nbins=maxticks, prune=prune)
      ax.yaxis.set_major_locator(loc)


def parseTimeStr( timeStr ):
  if len(timeStr.split('-')) == 3 :
    dt = datetime.datetime.strptime( timeStr ,'%Y-%m-%d')
  elif len(timeStr.split('-')) == 4 :
    dt = datetime.datetime.strptime( timeStr ,'%Y-%m-%d-%H')
  elif len(timeStr.split('-')) == 5 :
    dt = datetime.datetime.strptime( timeStr ,'%Y-%m-%d-%H-%M')
  elif len(timeStr.split('-')) == 6 :
    dt = datetime.datetime.strptime( timeStr ,'%Y-%m-%d-%H-%M-%S')
  else :
    raise Exception( 'Could not parse date string:',timeStr )
  return dt

def createDirectory(path) :
  if os.path.exists(path) :
    if not os.path.isdir(path) :
      raise Exception( 'file with same name exists',path )
  else :
    os.makedirs(path)
  return path

def saveFigure( path, filename, extensions, verbose=False, dpi=200, bbox_tight=False ) :
  """Saves currently open figure in the given format.
  If extensions is a list, figure is saved in multiple formats."""
  kw = {}
  if bbox_tight :
    kw['bbox_inches'] = 'tight'
  if not isinstance( extensions, list ) :
    extensions = [ extensions ]
  for ext in extensions :
    f = os.path.join(path,filename+'.'+ext )
    if verbose :
      print 'saving to', f
    plt.savefig(f, dpi=dpi, **kw)

def convertEpochToPlotTime( t ) :
  """Converts python datetime epoch value to num used in matplotlib.
  Epoch time is always in UTC (see data.timeArray)"""
  tzCorrection = -8.0*3600 # time zone correction for PST
  return matplotlib.dates.epoch2num( t + tzCorrection )

def convertPlotToEpochTime( t ) :
  """Converts python matplotlib num to datetime epoch time format.
  Epoch time is always in UTC (see data.timeArray)"""
  tzCorrection = -8.0*3600 # time zone correction for PST
  return matplotlib.dates.num2epoch( t ) - tzCorrection

#-------------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------

class plotBase(object) :

  """The mother of all plots"""
  def __init__(self, **defaultArgs) :
    """Constructor. defaultArgs are the default plot options."""
    if not hasattr(self,'xIsTime') :
      self.xIsTime = defaultArgs.pop('xIsTime',False)
    if not hasattr(self,'yIsTime') :
      self.yIsTime = defaultArgs.pop('yIsTime',False)
    self.defaultArgs = defaultArgs
    self.ax = None
    self.titleObj = None
    # default color sequence for line plots etc
    self.colorSequence = defaultColorSequence

  def createAxes(self, fig=None, rect=111) :
    """Create rectangular diagram axes for the diagram. Returns the axes.
    """
    if fig==None :
      fig = plt.figure()
    ax = fig.add_subplot(rect)
    return ax

  def setAxes(self, ax) :
    """Set axes for the diagram. All data will be plotted in these axes.
    """
    # this is a good spot for setting x/ylabel etc properties
    self.ax = ax

  def setupAxes(self, fig=None, rect=111) :
    """Creates new axes and set them.
    A shorthand for calling createAxes and setAxes.
    Returns the axes."""
    ax = self.createAxes(fig, rect)
    self.setAxes(ax)
    return ax

  def hideXTicks(self):
    """Removes xlabel and makes x axis tick labels invisible"""
    self.ax.set_xlabel('')
    plt.setp( self.ax.get_xticklabels(), visible=False)

  def hideYTicks(self):
    """Removes ylabel and makes y axis tick labels invisible"""
    self.ax.set_ylabel('')
    plt.setp( self.ax.get_yticklabels(), visible=False)

  def updateXAxis( self, xlim=None, minticks=3, maxticks=12, prune=None ) :
    """Updates x axis ticks.

    xlim sets the range for the axis, if xlim==None the range is not changed.
    If xlim='tight' x limits are set to data range.
    prune='lower'|'upper'|'both'|None can be used to omit low/high tick labels
    """
    updateXAxis(self.ax, self.xIsTime, xlim, minticks, maxticks, prune)

  def updateYAxis( self, ylim=None, minticks=3, maxticks=12, prune=None ) :
    """Updates y axis ticks.

    ylim sets the range for the axis, if ylim==None the range is not changed.
    If ylim='tight' y limits are set to data range.
    prune='lower'|'upper'|'both'|None can be used to omit low/high tick labels
    """
    updateYAxis(self.ax, self.yIsTime, ylim, minticks, maxticks, prune)

  def addSample(self, sample, **kwargs) :
    """Add signal to the diagram. Time is assumed to be in epoch format."""
    raise NotImplementedError( 'This method must be implemented in the derived class' )

  def addTitle(self, titleStr, **kwargs) :
    """Adds a title in the axes. Optional arguments are passed to set_title routine"""
    if self.ax is not None:
      self.titleObj = self.ax.set_title(titleStr, **kwargs)

  def updateTitle(self, titleStr):
    """Re-renders title string in animation mode."""
    if self.ax is not None and self.titleObj is not None:
      self.titleObj.set_text(titleStr)

  def showColorBar(self, **kwargs) :
    """Adds a colorbar in the axes."""
    if self.ax :
      # create an axes on the right side of ax. The width of cax will be 5%
      # of ax and the padding between cax and ax will be fixed at 0.05 inch.
      divider = make_axes_locatable(self.ax)
      cax = divider.append_axes("right", size=0.3, pad=0.1, add_to_figure=False)
      #cax.get_xaxis().set_visible(False)

  def addShadedRange(self, start, end, **kwargs) :
    """Adds a shaded time range in the background.

    If x-axis is time start,end can be a datetime object or float representing epoch time.
    kwargs are passed to matplotlib axvspan routine."""
    if self.xIsTime :
      if isinstance( start, datetime.datetime ) :
        start = timeArray.datetimeToEpochTime(start)
      if isinstance( end, datetime.datetime ) :
        end = timeArray.datetimeToEpochTime(end)
      start = convertEpochToPlotTime(start)
      end = convertEpochToPlotTime(end)
    kwargs.setdefault( 'facecolor', [0.8,1.0,0.85] )
    kwargs.setdefault( 'edgecolor', 'none' )
    self.ax.axvspan(start, end, **kwargs)

  def showLegend(self, *args, **kwargs) :
    """Adds a legendbox in the axes. Optional arguments are passed to legend routine"""
    if self.ax :
      if not 'loc' in kwargs :
        # place outside axis on left
        kwargs['loc'] = 'upper left'
        kwargs['bbox_to_anchor'] =(1.01, 1.)
      kwargs.setdefault('prop',{}).setdefault('size','small')
      self.ax.legend(*args, **kwargs)

class colorPlotBase(plotBase) :
  """Base class for all plots that need a colorbar
  (slab, track, transect, profile etc.)"""
  def __init__( self, **defaultArgs ) :
    self.cax = None # cax is the return value of scatter/mesh/contourf/tripcolor
    plotBase.__init__(self,**defaultArgs)

  def showColorBar(self, location='right', size=0.3, pad=0.1, fontsize=None,
                   maxticks=None, prune=None, integer=False, symmetric=False,
                   **kwargs) :
    """Adds a colorbar in the axes."""
    if self.ax and self.cax:
      #title = r'TKE [$\log_{10} (m^2 s^{-2}) $]'
      title = self.clabel
      if self.unit :
        title += r'$\ [$'+self.unit+'$]$'
      title = kwargs.pop( 'title',title  )
      # try to guess good format
      clim = self.cax.get_clim()
      magnitude = max([abs(v) for v in clim])
      default_format = '%.1f'
      if magnitude > 100 : default_format = None
      if magnitude < 0.5 : default_format = '%.2f'
      if magnitude < 0.05 : default_format = '%.1e'
      kwargs.setdefault('format',default_format)
      # create an axes on the right side of ax. The width of cax will be 5%
      # of ax and the padding between cax and ax will be fixed at 0.05 inch.
      divider = make_axes_locatable(self.ax)
      cax = divider.append_axes(location, size=size, pad=pad)
      orientation = 'vertical' if location in ['right','left'] else 'horizontal'
      self.cb = plt.colorbar(self.cax,cax=cax, orientation=orientation,
                             **kwargs)
      if maxticks is not None:
        loc = matplotlib.ticker.MaxNLocator(nbins=maxticks, prune=prune,
                                            integer=integer,
                                            symmetric=symmetric)
        self.cb.locator = loc
        self.cb.update_ticks()
      self.cb.set_label(title, fontsize=fontsize)

  def updateColorBar(self):
      """Updates the color bar. This is meant to be called when the associated
      plot has been updated trough updatePlot routine."""
      if self.cax:
        self.cax.autoscale()

class stackPlotBase(object) :
  """A class for stacking multiple transects in the same plot."""
  def __init__(self, **defArgs) :
    self.plotheight = defArgs.pop( 'plotheight' , 2 ) # ~height of each subplot
    self.figwidth = defArgs.pop( 'figwidth' , 10 ) # width of figure
    self.topPad = defArgs.pop( 'topPad' , 0.3 )
    self.bottomPad = defArgs.pop( 'bottomPad' , 0.25 )
    self.rightPad = defArgs.pop( 'rightPad' , 0.8 )
    self.sharex = defArgs.pop( 'sharex' , True )
    self.fig = plt.figure( figsize=(self.figwidth,self.plotheight) )
    self.plots = dict() # stores all plots
    self.tags = [] # to keep ordering
    self.axGrid = []
    self.defArgs = defArgs

  def addPlot(self, plot, tag=None) :
    """Adds a new subplot to the diagram"""
    if tag==None :
      tag = str(len(self.plots.keys())+1)
    if tag in self.plots.keys() :
      print 'Warning: subplot with tag {tag:s} already exists'.format(tag=tag)
    self.plots[tag]=plot
    self.tags.append( tag )
    # update figure height
    newHeight = len(self.tags)*self.plotheight +self.topPad+self.bottomPad
    self.fig.set_size_inches( (self.figwidth,newHeight), forward=True )

    # create new subplot
    n = len(self.axGrid)
    if n == 0 :
      newAx = self.fig.add_subplot(n+1, 1, n+1)
    else :
      # fix geometry of previous plots
      for i in range(n):
        self.axGrid[i].change_geometry(n+1, 1, i+1)
      # remove ticks and xlabel
      self.axGrid[-1].set_xlabel('')
      plt.setp( self.axGrid[-1].get_xticklabels(), visible=False)
      # new axes
      if self.sharex:
        newAx = self.fig.add_subplot(n+1, 1, n+1, sharex=self.axGrid[0])
      else:
        newAx = self.fig.add_subplot(n+1, 1, n+1)
    # add new axes
    self.axGrid.append( newAx )
    # assign axes
    self.plots[tag].setAxes( newAx )
    self.fig.tight_layout()
    w,h = self.fig.get_size_inches()
    topMargin = 1.0 - self.topPad/h
    rightMargin = 1.0 - self.rightPad/w
    self.fig.subplots_adjust(top=topMargin,right=rightMargin)
    return newAx

  #def isEmpty( self ) :
    #for tag in self.plots :
      #for err in self.plots[tag].data :
        #return False
    #return True

  def addTitle(self, titleStr, tag=None, **kwargs) :
    """Adds a title in the axes. Optional arguments are passed to set_title routine"""
    if tag == None :
      tag = self.tags[0]
    self.plots[tag].addTitle(titleStr, **kwargs)

  def updateTitle(self, titleStr, tag=None, **kwargs) :
    """Updates the title string"""
    if tag == None :
      tag = self.tags[0]
    self.plots[tag].updateTitle(titleStr)

  def showColorBar(self, **kwargs) :
    if 'tag' in kwargs :
      tag = kwargs.pop('tag')
      self.plots[tag].showColorBar(**kwargs)
    else :
      for tag in self.plots :
        #if hasattr(self.plots[tag],'showColorBar') :
        self.plots[tag].showColorBar(**kwargs)

  def updateColorBar(self, **kwargs) :
    if 'tag' in kwargs :
      tag = kwargs.pop('tag')
      self.plots[tag].updateColorBar()
    else :
      for tag in self.plots :
        #if hasattr(self.plots[tag],'showColorBar') :
        self.plots[tag].updateColorBar()

  def showLegend(self, *args, **kwargs) :
    for tag in self.plots :
      if hasattr(self.plots[tag],'showLegend') :
        self.plots[tag].showLegend(*args, **kwargs)
