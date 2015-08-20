#!/usr/bin/python

"""
A class for plotting transects.

Tuomas Karna 2012-11-16
"""

import numpy as np
import datetime

import matplotlib.pyplot as plt
from crane.data import timeArray
from crane.plotting import plotBase
from crane.plotting.plotBase import saveFigure, createDirectory

STA_MIN_DIST = 800

# class that takes the data (transect time series) and plots a snaphot (for certain time)
class transectSnapshot(plotBase.colorPlotBase) :
  """transectSnapshot class"""
  def __init__(self, **defaultArgs) :
    # default plot options for all diagrams
    self.unit =   defaultArgs.pop('unit')
    self.xlabel = defaultArgs.pop('xlabel','Distance')
    self.xunit =  defaultArgs.pop('xunit','km')
    self.ylabel = defaultArgs.pop('ylabel','Depth')
    self.yunit =  defaultArgs.pop('yunit','m')
    self.clabel = defaultArgs.pop('clabel')
    self.clim =   defaultArgs.pop('clim',[])
    self.climIsLog = defaultArgs.pop('climIsLog',False)
    self.ylim =   defaultArgs.pop('ylim',[])
    self.logScale = defaultArgs.pop('logScale',False)
    if self.logScale :
      self.unit = 'log10('+self.unit+')'
    super(transectSnapshot, self).__init__(**defaultArgs)
    self.ax = None
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

  def prepDataForPlotting(self, xCoord0, zCoord0, variable, **kwargs):
    xCoord = xCoord0.copy()
    zCoord = zCoord0.copy()
    kw = dict(self.defaultArgs)
    kw.update(kwargs)

    data = variable.copy() if not self.logScale else np.log10( variable )
    if self.clim :
      cmin, cmax = self.clim
      if self.logScale and not self.climIsLog :
        cmin, cmax = np.log10(cmin), np.log10(cmax)
      # saturate limits
      data[ data < cmin] = cmin
      data[ data > cmax] = cmax
    else :
      cmin, cmax = data.min(), data.max()+1e-9
    N = kw.pop('N', 20)
    kw.setdefault('levels', np.linspace( cmin, cmax, N ) )

    if self.xunit == 'km' :
      xCoord /= 1000
    if 'cmap' in kw and isinstance(kw['cmap'],str) :
      kw['cmap'] = plt.get_cmap(kw['cmap'])

    return xCoord, zCoord, data, kw

  def addSample(self, xCoord0, zCoord0, variable, **kwargs) :
    """Add transect to the diagram"""
    xCoord, zCoord, data, kw = self.prepDataForPlotting(xCoord0, zCoord0,
                                                        variable, **kwargs)
    plotType = kw.pop('plotType', 'contourf')
    if plotType == 'contourf' :
      self.cax = self.ax.contourf(xCoord, zCoord, data, **kw)
    else :
      self.cax = self.ax.contour(xCoord, zCoord, data, **kw)
    # tight y-axis
    ylim = [ self.ax.dataLim.ymin, self.ax.dataLim.ymax]
    if self.ylim :
      ylim = self.ylim
    self.ax.set_ylim( ylim )
    xlim = kw.pop('xlim',None)
    if xlim :
      self.ax.set_xlim( xlim )

  def addStationMarker(self, x, label, **kwargs) :
    """Adds a vertical line in the given position with a label text.
    kwargs are passed to ax.plot routine."""
    x0 = x
    if self.xunit == 'km' :
      x0 /= 1000
    xmin = self.ax.dataLim.xmin
    xmax = self.ax.dataLim.xmax
    ymin = self.ax.dataLim.ymin
    ymax = self.ax.dataLim.ymax
    if x0 >= xmin and x0 <= xmax :
      self.ax.axvline(x0,**kwargs)
      kwargs.pop('linewidth',None)
      kwargs.pop('linestyle',None)
      xRange = (xmax-xmin)
      hAlign = 'right' if (xmax-x0) < 0.05*xRange else 'left'
      self.ax.text( x0+0.005*xRange, ymin+0.01*(ymax-ymin), label, verticalalignment='bottom',horizontalalignment=hAlign, **kwargs )

def dcToTransect( dc2, iComp=0 ) :
  """
  Unravel dataContainer arrays for transect plots.
  
  Args:
   dc2    (dataContainer) representing transect
   iComp  (int) dataContainer field id to process (scalar=0, u=0, v=1)
  Outputs:
   T.shape (nTime,) epoch time stamps
   X.shape (nVertical,nPoints,) transect x coordinates
   Y.shape (nVertical,nPoints,) transect y coordinates
   Z.shape (nVertical,nPoints,nTime) transect z coordinates
   C.shape (nVertical,nPoints,nTime) values from iComp component
  
  Here nVertical is the max nb of vertical levels in the dataset.
  If not uniform, the arrays are padded with nans.
  """
  x = dc2.x[:,0] if dc2.xDependsOnTime else dc2.x
  y = dc2.y[:,0] if dc2.yDependsOnTime else dc2.y
  if np.isnan(dc2.x).any() or np.isnan(dc2.x).any():
    raise Exception('This method cannot handle nans in x or y')
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

  nTime = dc2.z.shape[1]
  nPoints = len(breakIx)-1
  maxNVert = max( breakIx[i+1]-breakIx[i] for i in range(nPoints) )
  T = dc2.time.array
  dtype = dc2.data.dtype
  X = (np.ones( (maxNVert,nPoints) )*np.nan).astype(dtype)
  Y = (np.ones( (maxNVert,nPoints) )*np.nan).astype(dtype)
  Z = (np.ones( (maxNVert,nPoints,nTime) )*np.nan).astype(dtype)
  C = (np.ones( (maxNVert,nPoints,nTime) )*np.nan).astype(dtype)
  for i in range(nPoints) :
    n = breakIx[i+1]-breakIx[i]
    X[:n,i] = dc2.x[breakIx[i]:breakIx[i+1]]
    Y[:n,i] = dc2.y[breakIx[i]:breakIx[i+1]]
    Z[:n,i,:] = dc2.z[breakIx[i]:breakIx[i+1],:]
    C[:n,i,:] = dc2.data[breakIx[i]:breakIx[i+1],iComp,:]
  return T,X,Y,Z,C

def generateTransectFromDataContainer( dc, timeStamp, flipDirection=False ) :
  """Reshapes and massages dataContainer data into transect grids.
  """
  # TODO save transects as masked array in vector
  # TODO optimize
  if isinstance( timeStamp, datetime.datetime ) :
    # interpolate to correct time
    timeStamp = timeArray.datetimeToEpochTime( timeStamp )
    newTime = timeArray.timeArray( np.array( [timeStamp] ), 'epoch' )
    dc = dc.interpolateInTime( newTime, acceptNaNs=True )
    x = dc.x.flatten()
    y = dc.y.flatten()
    z = dc.z.flatten()
    data = dc.data[:,:,0] # shape (npts,nComp)
    time = dc.time.getDatetime(0)
  else :
    x = dc.x[:,timeStamp] if dc.xDependsOnTime else dc.x
    y = dc.y[:,timeStamp] if dc.yDependsOnTime else dc.y
    z = dc.z[:,timeStamp] if dc.zDependsOnTime else dc.z
    data = dc.data[:,:,timeStamp] # shape (npts,nComp)
    time = dc.time.getDatetime(timeStamp)
  # find unique x,y pairs
  xyHashed = x + 1e6*y
  uniqueXYHash, uniqueIx = np.unique(xyHashed,return_index=True)
  uniqueIx.sort() # keep ordering
  uniqueXYHash = xyHashed[uniqueIx]
  uniqueXYCoords = np.vstack( (x[uniqueIx],y[uniqueIx]) ).T
  # covert to grids
  # in dataContainer transects are stored as vectors
  nx = uniqueXYCoords.shape[0]
  nVert = np.array([ len(np.nonzero( xyHashed == uniqueXYHash[i] )[0]) for i in range(nx) ])
  maxNVert = nVert.max()
  nComp = data.shape[1]
  X = np.ma.ones( (maxNVert,nx) )*np.ma.masked
  Y = np.ma.ones( (maxNVert,nx) )*np.ma.masked
  Z = np.ma.ones( (maxNVert,nx) )*np.ma.masked
  C = np.ma.ones( (maxNVert,nx,nComp) )*np.ma.masked
  for i in range(nx) :
    ix = np.nonzero( xyHashed == uniqueXYHash[i] )[0]
    X[:,i] = uniqueXYCoords[i,0] # fill all levels
    Y[:,i] = uniqueXYCoords[i,1] # fill all levels
    Z[-len(ix):,i] = z[ix]
    C[-len(ix):,i,:] = data[ix,:]
  # flip transect direction
  if flipDirection :
    X = X[:,::-1]
    Y = Y[:,::-1]
    Z = Z[:,::-1]
    C = C[:,::-1,:]
  Z = np.ma.masked_invalid( Z )
  Xdiff = np.diff(X,axis=1)
  Ydiff = np.diff(Y,axis=1)
  Xdelta = np.hstack( ( np.zeros((Xdiff.shape[0],1)), Xdiff ) )
  Ydelta = np.hstack( ( np.zeros((Ydiff.shape[0],1)), Ydiff ) )
  if nComp == 2 :
    # compute tangential velocity
    # mean deltaX for each node
    Xdelta2 = np.hstack( ( Xdiff[:,1,None], Xdiff ) )/2
    Ydelta2 = np.hstack( ( Ydiff[:,1,None], Ydiff ) )/2
    Xdelta2 += np.hstack( ( Xdiff, Xdiff[:,-1,None] ) )/2
    Ydelta2 += np.hstack( ( Ydiff, Ydiff[:,-1,None] ) )/2
    # compute along transect unit vectors
    XYmag = np.sqrt( Xdelta2**2 + Ydelta2**2 ) # shape (maxNVert,nx )
    # shape (maxNVert,nx,2)
    XYunit = np.ma.zeros((maxNVert,nx,2))
    XYunit[:,:,0] = Xdelta2/XYmag # preserves mask
    XYunit[:,:,1] = Ydelta2/XYmag
    # tangential velocity
    C = np.sum(C * XYunit, axis=2)
  elif nComp == 1 :
    C = C[:,:,0]
  else :
    raise Exception( 'nFields > 2 is not supported' )
  C = np.ma.masked_invalid( C )
  # compute along-channel coordinates
  Xalong = np.sqrt( Xdelta**2 + Ydelta**2 )
  Xalong = np.cumsum( Xalong, axis=1 )

  return Xalong, Z, C, time, uniqueXYCoords

def projectXYOnTransect( x,y,xTran,yTran,xAlong ) :
  """Computes an approximation for x,y in the along-transect coordinates.
  Transect is defined by xTran,yTran coordinates, the along-transect
   coordinates of these points are xAlong."""
  dist = np.sqrt( (xTran-x)**2 + (yTran-y)**2 )
  nearestPt = np.argmin(dist)
  minDist = dist[nearestPt]
  xDiff = xTran[nearestPt] - x
  yDiff = yTran[nearestPt] - y
  if nearestPt < len(xTran) -1 :
    xVec = xTran[nearestPt+1] - xTran[nearestPt]
    yVec = yTran[nearestPt+1] - yTran[nearestPt]
  else :
    xVec = xTran[nearestPt] - xTran[nearestPt-1]
    yVec = yTran[nearestPt] - yTran[nearestPt-1]
  diffAlong = (xDiff*xVec + yDiff*yVec)/np.sqrt( xVec**2+yVec**2 )
  xNew = xAlong[nearestPt]
  return xNew, minDist
  
class transectSnapshotDC(transectSnapshot) :
  """Transect plot class that uses dataContainer as input."""
  def __init__(self, **defaultArgs) :
    self.xyPoints = None
    self.xAlong = None
    transectSnapshot.__init__(self, **defaultArgs)

  def addSample( self, dc, timeStamp, **kwargs ) :
    """Add transect from dataContainer wih appropriate data restructuring.
       Args:
       dc        -- (dataContainer) data
       timeStamp -- (int or datetime) time index/stamp to plot.
                    If datetime object, temporal interpolation is performed.
       kwargs    -- (**) passed to transectSnapshot.addSample
    """
    flipDir = kwargs.get('flipDirection',False)
    Xalong, Z, C, time, uniqueXYCoords = generateTransectFromDataContainer( dc, timeStamp, flipDir )
    # replace mask with nans for plotting
    Z = Z.filled( -1e4 )
    # for mapping coordinates on transect
    self.xyPoints = uniqueXYCoords
    self.xAlong = Xalong[0,:]
    
    if not self.ylim :
      self.ylim = [np.nanmin(dc.z),np.nanmax(dc.z)] 
    if not 'clim' in kwargs or not kwargs['clim'] :
      if dc.fieldNames == ['u','v'] :
        maxAbs = max( abs(dc.data.min()),abs(dc.data.max()) )
        kwargs['clim'] = [-maxAbs,maxAbs]
      else :
        kwargs['clim'] = [dc.data.min(),dc.data.max()]
    transectSnapshot.addSample( self, Xalong, Z, C, **kwargs )
    self.addTitle( dc.description.split('.')[0]+' '+time.strftime( '%Y-%m-%d %H:%M:%S' )+' (PST)' )

  def addStationMarker(self, x, y, label, staMinDist=STA_MIN_DIST, **kwargs) :
    """Adds a vertical line in the given position with a label text.
    If station is more than staMinDist from any transect point, the marker is not drawn.
    kwargs are passed to ax.plot routine."""
    xAlong, minDist = projectXYOnTransect( x,y,self.xyPoints[:,0],self.xyPoints[:,1],self.xAlong )
    if minDist > staMinDist :
      print 'point too far from the transect, omitting',label,minDist,staMinDist
      return
    # if close enough, just use x coordinate
    transectSnapshot.addStationMarker( self, xAlong, label, **kwargs )

class stackTransectPlot(plotBase.stackPlotBase) :
  """A class for stacking multiple transects in the same plot."""

  def addPlot(self, tag, **kwargs) :
    """Adds a new subplot to the diagram"""
    kw = dict(self.defArgs)
    kw.update(kwargs)
    plot = transectSnapshot(**kw)
    super(stackTransectPlot, self).addPlot(plot,tag)

  def addSample( self, tag, xCoord, zCoord, variable, **kwargs ) :
    if not tag in self.tags :
      self.addPlot(tag, **kwargs)
    self.plots[tag].addSample(xCoord, zCoord, variable, **kwargs)

  def addStationMarker(self, tag, x, label, **kwargs) :
    """Adds a vertical line in the given position with a label text.
    If tag is 'all', marker is added to all subplots.
    kwargs are passed to ax.plot routine."""
    if tag != 'all' :
      self.plots[tag].addStationMarker(x, label, **kwargs)
    else :
      for t in self.plots : # recursion
        self.addStationMarker(t, x, label, **kwargs)

  ## TODO is this needed ???
  #def isEmpty( self ) :
    #for tag in self.plots :
      #for err in self.plots[tag].data :
        #return False
    #return True

class stackTransectPlotDC(stackTransectPlot) :
  """stackTransectPlot class that uses dataContainer as an input"""
  def __init__(self, **defArgs) :
    self.xyPoints = {} # for each tag
    self.xAlong = {}
    stackTransectPlot.__init__(self,**defArgs)
    
  def addSample( self, tag, dc, timeStamp, **kwargs ) :
    """Add transect from dataContainer wih appropriate data restructuring.
       Args:
       tag       -- (string) tag for identifiying the subplot
       dc        -- (dataContainer) data
       timeStamp -- (int or datetime) time index/stamp to plot.
                    If datetime object, temporal interpolation is performed.
       kwargs    -- (**) passed to transectSnapshot.addSample
    """
    flipDir = kwargs.get('flipDirection',False)
    Xalong,Z,C,time,uniqueXYCoords = generateTransectFromDataContainer( dc, timeStamp, flipDir )
    # replace mask with nans for plotting
    C = C.filled(np.nan)
    Z = Z.filled(-1e4)
    self.xyPoints[tag] = uniqueXYCoords
    self.xAlong[tag] = Xalong[0,:]
    if not 'ylim' in kwargs or not kwargs['ylim'] :
      kwargs['ylim'] = [np.nanmin(dc.z),np.nanmax(dc.z)]
    if not 'clim' in kwargs or not kwargs['clim'] :
      if dc.fieldNames == ['u','v'] :
        maxAbs = max( abs(dc.data.min()),abs(dc.data.max()) )
        kwargs['clim'] = [-maxAbs,maxAbs]
      else :
        kwargs['clim'] = [dc.data.min(),dc.data.max()]
    stackTransectPlot.addSample( self, tag, Xalong, Z, C, **kwargs )

  def addStationMarker(self, tag, x, y, label, staMinDist=STA_MIN_DIST, **kwargs) :
    """Adds a vertical line in the given position with a label text.
    If station is more than staMinDist from any transect point, the marker is not drawn.
    If tag is 'all', marker is added to all subplots.
    kwargs are passed to ax.plot routine."""
    if tag != 'all' :
      xPts = self.xyPoints[tag][:,0]
      yPts = self.xyPoints[tag][:,1]
      xAlo = self.xAlong[tag]
      xAlong, minDist = projectXYOnTransect( x,y,xPts,yPts,xAlo )
      if minDist > staMinDist :
        print 'point too far from the transect, omitting',tag,label,minDist,staMinDist
        return
      # if close enough, just use x coordinate
      stackTransectPlot.addStationMarker( self, tag, xAlong, label, **kwargs )
    else :
      for t in self.plots : # recursion
        self.addStationMarker( t,x,y,label,staMinDist, **kwargs )

# TODO add acrossVel, magnitude etc computation
# TODO add station etc position markers
# TODO convert m to km if suitable

if __name__ == '__main__' :
  # generate dummy data
  x = np.linspace(1e3,5e3,100)
  z = np.linspace(-20,0,40)
  [X,Z]=np.meshgrid(x,z)
  C = ( np.sin(X/400)*(-Z) - z.min() ) / 2

  fig = plt.figure()
  dia = transectSnapshot(clabel='Salinity',unit='psu')
  dia.setupAxes( fig )
  dia.addSample( X, Z, C, s=30 )
  dia.addStationMarker( 3555, 'sta1', color='k', linewidth=2.0, linestyle='dashed' )
  dia.showColorBar( )
  dia.addTitle('transectSnapshot example')
  plt.show()

  dia = stackTransectPlot()
  dia.addSample( 'first', X, Z, C, clabel='Salinity',unit='psu' )
  dia.addSample( 'second', X, Z, C, clabel='Salinity',unit='psu' )
  dia.addStationMarker( 'first', 3555, 'sta1', color='k', linewidth=2.0, linestyle='dashed' )
  dia.addStationMarker( 'all', 4555, 'sta2', color='k', linewidth=1.5, linestyle='dashed' )
  dia.addTitle('stackTransectPlot example')
  plt.show()

  from crane.data import dataContainer
  d0 = dataContainer.dataContainer.loadFromNetCDF('~/workspace/cmop/projects/db29dev/transect/run76/data/transect/nchannel_salt_tran_2012-05-01_2012-05-19.nc')
  print d0

  fig = plt.figure(figsize=(15,5))
  dia = transectSnapshotDC(clabel='Salinity',unit='psu', clim=[0,34])
  dia.setupAxes( fig )
  dia.addSample( d0,354, N=35 )
  #dia.addSample( d0,datetime.datetime(2012,5,12,11,12,34), N=35 )
  dia.addStationMarker( x=349653.0,y=290854.0, label='sat01', color='k', linewidth=2.0, linestyle='dashed' )
  dia.addStationMarker( x=344755, y=291695, label='ncbn1', color='k', linewidth=2.0, linestyle='dashed' )
  dia.showColorBar( )
  plt.show()

  dia = stackTransectPlotDC()
  dia.addSample( 'first', d0,datetime.datetime(2012,5,12,11,12,34), N=35, clim=[0,34], clabel='Salinity',unit='psu' )
  dia.addSample( 'second', d0, 500, N=35, clim=[], clabel='Salinity',unit='psu', logScale=True )
  dia.addSample( 'third', d0, 510, N=35, clim=[], clabel='Salinity',unit='psu', logScale=True )
  dia.addSample( '4', d0, 520, N=35, clim=[], clabel='Salinity',unit='psu', logScale=True )
  dia.addSample( '5', d0, 530, N=35, clim=[], clabel='Salinity',unit='psu', logScale=True )
  dia.addStationMarker( 'first', x=349653.0,y=290854.0, label='sat01', color='k', linewidth=2.0, linestyle='dashed' )
  dia.addStationMarker( 'all', x=344755, y=291695, label='ncbn1', color='k', linewidth=1.5, linestyle='dashed' )
  dia.addTitle('stackTransectPlotDC example')
  plt.show()
