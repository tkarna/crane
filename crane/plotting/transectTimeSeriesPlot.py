#!/usr/bin/python
"""
A class for plotting transect timeseries derived from transectPlot.py

Jesse Lopez 2014-06-04
"""
#------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------
from plotting.plotBase import *
from data.dataContainer import dataContainer
from plotting.transectPlot import projectXYOnTransect


#------------------------------------------------------------------------------
# Classes 
#------------------------------------------------------------------------------
class transectTimeSeriesSnapShot(colorPlotBase) :
  """transectTimeSeriesSnapShot class"""
  def __init__(self, **defaultArgs) :
    # More legible for plotting
    plt.rcParams['figure.figsize'] = 10, 8
    # default plot options for all diagrams
    self.unit =   defaultArgs.pop('unit')
    self.xlabel = defaultArgs.pop('xlabel','Date')
    self.xunit =  defaultArgs.pop('xunit', None)
    self.ylabel = defaultArgs.pop('ylabel','Distance along transect')
    self.yunit =  defaultArgs.pop('yunit','km')
    self.clabel = defaultArgs.pop('clabel')
    self.clim =   defaultArgs.pop('clim',[])
    self.climIsLog = defaultArgs.pop('climIsLog',False)
    self.ylim =   defaultArgs.pop('ylim',[])
    self.logScale = defaultArgs.pop('logScale',False)
    defaultArgs.setdefault('xIsTime',True)
    defaultArgs.setdefault('yIsTime',False)
    if self.logScale :
      self.unit = 'log10('+self.unit+')'
    colorPlotBase.__init__(self,**defaultArgs)
    self.ax = None
    self.pax = None
    self.cax = None

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

  def addSample(self, time, x, variable, **kwargs):
    """Add transect to the diagram"""
    kw = dict(self.defaultArgs)
    kw.update(kwargs)
    
    plotType = kw.pop('plotType', 'contourf')

    # convert inputs for plot here!
    t = convertEpochToPlotTime(time)
    x_conv = x.copy()
    if not self.xIsTime and self.xunit == 'km':
        x_conv /= 1000.0
    if not self.yIsTime and self.yunit == 'km':
        x_conv /= 1000.0
    data = variable.copy() if not self.logScale else np.log10( variable )

    clim = kw.pop('clim', self.clim)
    if clim:
        cmin, cmax = clim
        if self.logScale and not self.climIsLog:
            cmin, cmax = np.log10(cmin), np.log10(cmax)
        # saturate limits
        data[ data < cmin] = cmin
        data[ data > cmax] = cmax
    else:
        cmin, cmax = data.min(), data.max()+1e-9
    N = kw.pop('N', 20)
    kw.setdefault('levels', np.linspace( cmin, cmax, N ) )

    # Handle colormap here
    if 'cmap' in kw and isinstance(kw['cmap'],str) :
      kw['cmap'] = plt.get_cmap(kw['cmap'])
    if plotType == 'contour' and 'colors' in kw:
        kw.pop('cmap', None)
    #self.pax = self.ax.pcolormesh(data, shading='gouraud', **kw)
    if self.xIsTime:
        if plotType == 'contourf':
            self.cax = self.ax.contourf(t, x_conv, data, shading='gouraud', **kw)
        else:
            self.cax = self.ax.contour(t, x_conv, data, **kw)
        self.updateXAxis()
    else:
        if plotType == 'contourf':
            self.cax = self.ax.contourf(x_conv, t, data, shading='gouraud', **kw)
        else:
            self.cax = self.ax.contour(x_conv, t, data, **kw)
        self.updateYAxis()

  def addContourSample(self, variable, **kwargs) :
    """Adds contours of another varible to plot"""
    kw = dict(self.defaultArgs)
    kw.update(kwargs)

    data = variable.copy()
    self.cax = self.ax.contour(data, colors='k', linewidth=1, **kw)
    self.cax.clabel(inline=1, fontsize=12, fmt='%d')

  def addStationMarker(self, x, label, textkw={}, **kwargs):
    """Adds a vertical line in the given position with a label text.
    kwargs are passed to ax.plot routine."""
    x0 = x
    if self.xunit == 'km':
        x0 /= 1000
    xmin = self.ax.dataLim.xmin
    xmax = self.ax.dataLim.xmax
    ymin = self.ax.dataLim.ymin
    ymax = self.ax.dataLim.ymax
    if x0 >= xmin and x0 <= xmax:
        self.ax.axvspan(x0, x0, **kwargs)
        xRange = (xmax-xmin)
        hAlign = 'right' if (xmax-x0) < 0.05*xRange else 'left'
        tkw = {}
        tkw['color'] = kwargs['color']
        tkw['verticalalignment'] = 'bottom'
        tkw['horizontalalignment'] = hAlign
        tkw.update( textkw )
        self.ax.text(x0+0.005*xRange, ymin+0.01*(ymax-ymin), label, **tkw)


# TODO: Eliminate assumptions in the method
# XXX: Hacks
  def makePlotPretty(self, nNodes, dist, ndays, dates, **kwargs) :
    """Adds changes to plots specific to this class"""
    # Handle ticks on pcolormesh
    # This number looks "good", put at module level as PCOLORMESH_YTICK_STRIDE
    N = 10
    # Default levels here are for salt.  Put at module level
    SALT_LEVELS = np.arange(5, 31, 5)
    # Day tick factory?
    DAILY_TICKS = np.arange(24, 24*(ndays+1), 24) 
    XLABEL_FONTSIZE = 14
    YLABEL_FONTSIZE = 14
    CMIN, CMAX = 0.0, 0.1

    kw = dict(self.defaultArgs)
    kw.update(kwargs)
    kw.setdefault('levels', SALT_LEVELS) 
    self.pax.axes.set_yticklabels(np.floor(dist[::N]))
    self.pax.axes.set_yticks(np.arange(nNodes)[::N])

    # XXX: Assuming days along x-axis, place 
    self.pax.axes.set_xticks(DAILY_TICKS)
    self.pax.axes.set_xticklabels(dates[1:])

    self.ax.set_ylabel(self.ylabel, fontsize=YLABEL_FONTSIZE)
    self.ax.set_xlabel(self.xlabel, fontsize=XLABEL_FONTSIZE)

    if self.clim :
      cmin, cmax = self.clim
      if self.logScale and not self.climIsLog :
        cmin, cmax = np.log10(cmin), np.log10(cmax)
      # saturate limits
      data[ data < cmin] = cmin
      data[ data > cmax] = cmax
    else :
      cmin, cmax = CMIN, CMAX 
    self.pax.set_clim(vmin=cmin, vmax=cmax)
    N = kw.pop('N', 10)
    kw.setdefault('levels', SALT_LEVELS)


#------------------------------------------------------------------------------
# Functions 
#------------------------------------------------------------------------------
def generateTransectTimeSeriesFromDC( dc, level, iComp, average=4):
  """ Unravel dataContainer arrays for transect time series plots.
  
  params: 
  -------
  dc : dataContainer
    dataContainer of transect data
  level : int
    Vertical level of transect to plot (bottom=, ..., surface=-1) 
  iComp : int
    dataContainer field to process (scalar=0, u=0, v=1)
  average : int
    Number of timesteps to average over (noAveraging=1)

  returns:
  --------
  C : np.array
    Reshaped array of size (n_transect_nodes, times)
  T : timeArray
    Times of the data
  """

  if dc.xDependsOnTime or dc.yDependsOnTime:
    raise Exception('This method cannot handle time varying x or y')
  if np.isnan(dc.x).any() or np.isnan(dc.x).any():
    raise Exception('This method cannot handle nans in x or y')
  x = dc.x
  y = dc.y
  # find unique x,y pairs
  xyHashed = x + 1e6*y
  uniqueXYHash, uniqueIx = np.unique(xyHashed,return_index=True)
  uniqueIx.sort() # keep ordering
  uniqueXYHash = xyHashed[uniqueIx]
  breakIx = []
  uniXYIter = iter(uniqueXYHash)
  nextXY = uniXYIter.next()
  for i in range(len(x)) :
    if xyHashed[i] == nextXY :
      breakIx.append( i )
      try :
        nextXY = uniXYIter.next()
      except StopIteration :
        break
  breakIx.append( len(x) )
  breakIx = np.array(breakIx)

  nTime = dc.z.shape[1]
  nPoints = len(breakIx)-1
  NVert = [breakIx[i+1] - breakIx[i] for i in range(nPoints)]
  NVert.sort()
  if NVert[0] != NVert [-1]:
    raise Exception('This method cannot handle varying numbers of vertical levels')
  T = dc.time.array
  ix = np.arange(level, dc.data.shape[0], NVert[0]) 
  
  # Reshape data to (nPoints, time) with averaging over 'average' time steps 
  C = dc.data[ix, iComp, :].reshape(nPoints, average, -1).mean(axis=1)
  T = dc.time[::average]

  return C, T


def calcDistanceAlongTransect(dc) :
  """Calculate the distances along the transect """
  xpos, ix = np.unique(dc.x, return_index=True)
  ypos = dc.y[ix] 
  dpos = np.hypot(np.diff(xpos), np.diff(ypos))
  # In KM
  dist = np.cumsum(dpos)/1000.0

  return dist


def calcDateStrings(dc) :
  """Calculate datestrings and number of days for xticklabels """
  baseDate = timeArray.epochToDatetime(dc.time[0])
  endDate = timeArray.epochToDatetime(dc.time[-1])
  ndays = (endDate - baseDate).days

  dates = [baseDate + datetime.timedelta(days=d) for d in range(ndays+1)]
  dates = [d.strftime('%b %d') for d in dates]
  
  return dates, ndays


def calcNumberTransectPoints(dc) :
  """Calculate number of points in transect""" 
  x = dc.x
  y = dc.y
  # find unique x,y pairs
  xyHashed = x + 1e6*y
  uniqueXYHash, uniqueIx = np.unique(xyHashed,return_index=True)
  uniqueIx.sort() # keep ordering
  uniqueXYHash = xyHashed[uniqueIx]
  breakIx = []
  uniXYIter = iter(uniqueXYHash)
  nextXY = uniXYIter.next()
  for i in range(len(x)) :
    if xyHashed[i] == nextXY :
      breakIx.append( i )
      try :
        nextXY = uniXYIter.next()
      except StopIteration :
        break
  breakIx.append( len(x) )
  breakIx = np.array(breakIx)

  nPoints = len(breakIx)-1

  return nPoints

#------------------------------------------------------------------------------
# Main 
#------------------------------------------------------------------------------
if __name__ == '__main__':
  import os
  basePath = '/disk/ambfs21/0/users/lopezj/sediment/etm_cruise_runs/db22/fall/run08/run08/data/transect'
  sed_path = os.path.join(basePath, 'north_channel_ex_turbidity_0_2012-10-22_2012-11-01.nc')
  salt_path = os.path.join(basePath, 'north_channel_ex_salt_0_2012-10-22_2012-11-01.nc')

  bottom_level = 0
  field = 0
  # Average over 4 timesteps (hourly)
  avg = 4
  sed = dataContainer.loadFromNetCDF(sed_path)
  sed_trans, sed_time = generateTransectTimeSeriesFromDC(sed, bottom_level, field,
                                                         average=avg)
  dist = calcDistanceAlongTransect(sed)
  dates, ndays = calcDateStrings(sed)
  nPoints = calcNumberTransectPoints(sed)

  salt = dataContainer.loadFromNetCDF(salt_path)
  salt_trans, salt_time = generateTransectTimeSeriesFromDC(salt, bottom_level, field,
                                                           average=avg)

  fig = plt.figure(figsize=(10,8))
  dia = transectTimeSeriesSnapShot(clabel='SSC', unit='kg/m3')
  dia.setupAxes( fig )
  dia.addSample(sed_trans)
  dia.showColorBar()
  dia.addContourSample(salt_trans)
  dia.addTitle('North Channel\nNear-bottom suspended sediment')
  dia.makePlotPretty(nPoints, dist, ndays, dates)
  plt.show()
