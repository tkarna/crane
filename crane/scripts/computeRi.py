"""
A script for computing gradient Richardson number from Winched Profiler track.
"""

import numpy as np
from crane.data import dataContainer
from crane.data import timeArray
from data.stationCollection import *
import data.dirTreeManager as dtm
from data.pca import *
from data.eqState import *
from scipy.interpolate import interp1d, griddata
from scipy.signal import convolve
from scipy.spatial import cKDTree
from plotting.plotBase import convertEpochToPlotTime

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------

oldRule = dtm.oldTreeRule()
defRule = dtm.defaultTreeRule()

tra = dtm.netcdfTreeTraverser( rule=oldRule )
tra2 = dtm.netcdfTreeTraverser( rule=defRule )

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def doPCA( dc ) :
  A = dc.data[0,:2,:].T.copy()
  #print A.mean(axis=0)
  #print A.shape
  p = PCA(A,fraction=0.9999)
  print "npc:", p.npc
  print "% variance:", p.sumvariance * 100
  B = p.pc()
  #print B.shape
  dc2 = dc.copy()
  dc2.data = B.T[None,:,:]
  dc2.fieldNames = dc.fieldNames[:dc2.data.shape[1]]
  dc2.setMetaData('variable','pcavel')
  return dc2

def interpolateOnGrid( dc, T, Z ) :
  timeRange = (dc.time.array.max() - dc.time.array.min())/1000
  X = np.vstack((dc.time.array/timeRange,dc.z[0,:])).T
  X2 = np.vstack((T.ravel()/timeRange,Z.ravel())).T
  Y = dc.data[0,0,:][:,None]
  #print X.shape, Y.shape, X2.shape
  #print np.nanmin(Y), np.nanmax(Y)
  val_int = griddata( X,Y, X2 )
  val_int = np.reshape(val_int,T.shape)
  #print val_int.shape, np.nanmin(val_int), np.nanmax(val_int)
  return val_int

def binOnGrid( dc, T, Z, tScale, zScale ) :
  # construct a tree for finding nearest point in grid
  points = np.vstack( (T.ravel()/tScale,Z.ravel()/zScale) ).T
  tree = cKDTree( points )
  # find nearest point for each input point
  q = np.vstack((dc.time.array/tScale,dc.z[0,:]/zScale)).T
  vals = dc.data[0,0,:]
  k = 5
  nix = tree.query( q,k=k )[1]
  multip = np.zeros_like( T.ravel() )
  valsum = np.zeros_like( T.ravel() )
  for i in range(len(nix)) :
    for j in range(k) :
      multip[nix[i,j]] += 1.0/(j+1)
      valsum[nix[i,j]] += 1.0/(j+1)*vals[i]
    #multip[nix[i,1]] += 0.5
    #valsum[nix[i,1]] += 0.5*vals[i]
  goodIx = multip > 0
  valsum[goodIx] = valsum[goodIx]/multip[goodIx]
  valsum[~goodIx] = np.nan
  # fill gaps by interpolation
  q = np.vstack((dc.time.array,dc.z[0,:])).T
  valsum[~goodIx] = griddata( q[goodIx,:],valsum[goodIx], q[~goodIx,:] )
  return np.reshape(valsum,T.shape)

def interpolateOnGridDC( dc, T, Z ) :
  vals = interpolateOnGrid(dc,T,Z)
  ta = timeArray(T.ravel(),'epoch',acceptDuplicates=True)
  meta = dc.getMetaData()
  fn = dc.fieldNames[:1]
  dc2 = dataContainer('',ta, dc.x, dc.y, Z.ravel()[None,:], vals.ravel()[None,None,:],fn,metaData=meta, acceptNaNs=True)
  return dc2

# not needed
def runningX( dc, T=None, operator=computeRunningMean ) :
  if T == None :
    T = 44714. # M2 period in seconds
  gaps,ranges,t = dc.detectGaps()
  x = dc.data[0,0,:]
  z = dc.z[0,:]
  tRes = np.array([])
  xRes = np.array([])
  zRes = np.array([])
  for i in range( ranges.shape[0] ) :
    tt = t[ ranges[i,0]:ranges[i,1] ]
    xx = x[ ranges[i,0]:ranges[i,1] ]
    zz = z[ ranges[i,0]:ranges[i,1] ]
    if len(tt) == 0 :
      continue
    t0, xx = operator( tt,xx,T )
    tt, zz = computeRunningMean( tt,zz,T )
    tRes = np.append( tRes, tt )
    xRes = np.append( xRes, xx )
    zRes = np.append( zRes, zz )
  if len(tRes) == 0 :
    print 'Running mean could not be computed, skipping (time series too short?)'
    return
  ta = timeArray( tRes, 'epoch' )
  data = xRes.reshape( (1,1,-1) )
  z = zRes[None,:]
  meta = dc.getMetaData()
  dc2 = dataContainer( '', ta, dc.x,dc.y,z, data,
                            dc.fieldNames[:1], coordSys=dc.coordSys,metaData=meta)
  return dc2

def runningMean( dc, T=None ) :
  return runningX( dc, T, operator=computeRunningMean )

def runningMin( dc, T=None ) :
  return runningX( dc, T, operator=computeRunningMin )

def runningMax( dc, T=None ) :
  return runningX( dc, T, operator=computeRunningMax )

def runningRange( dc, T=None ) :
  return runningX( dc, T, operator=computeRunningRange )

def gauss_kern(size, sizey=None):
  """ Returns a normalized 2D gauss kernel array for convolutions """
  size = int(size)
  if not sizey:
    sizey = size
  else:
    sizey = int(sizey)
  x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
  g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
  return g / g.sum()

def blur_image(im, n, ny=None) :
  """ blurs the image by convolving with a gaussian kernel of typical
      size n. The optional keyword argument ny allows for a different
      size in the y direction.
  """
  g = gauss_kern(n, sizey=ny)
  #improc = convolve(im,g, mode='same')
  improc = convolve(im,g, mode='valid')
  return improc

def computeRi( tag, loc ) :
  st = datetime.datetime(2012,10,26)
  et = datetime.datetime(2012,11,3)

  hvel = dtm.getDataContainer( dataType='track', tag=tag, location=loc,variable='hvel', startTime=st, endTime=et )
  salt = dtm.getDataContainer( dataType='track', tag=tag, location=loc,variable='salt', startTime=st, endTime=et )
  temp = dtm.getDataContainer( dataType='track', tag=tag, location=loc,variable='temp', startTime=st, endTime=et )

  pc_u = doPCA( hvel )
  # remove outliers from second dir
  good_ix = np.abs(pc_u.data[0,1,:]) < 0.8
  pc_u.time.array = pc_u.time.array[good_ix]
  pc_u.data = pc_u.data[:,:,good_ix]
  pc_u.z = pc_u.z[:,good_ix]
  pc_u = doPCA( pc_u )
  # rotate velocity so that flood is positive
  ix = pc_u.z[0,:] < 10
  mvel = pc_u.data[0,0,ix].mean()
  print 'mean surf vel',mvel
  if mvel > 0 :
    print 'flipping'
    pc_u.data *= -1

  # generate mesh for interpolation
  # find top and bottom data boundaries
  tzmin, zmin = computeRunningMin( pc_u.time.asEpoch().array, pc_u.z.ravel(), T=45*60.0 )
  tzmax, zmax = computeRunningMax( pc_u.time.asEpoch().array, pc_u.z.ravel(), T=45*60.0 )
  tzmin, zmin = computeRunningMean( tzmin, zmin, T=2*3600.0 )
  tzmax, zmax = computeRunningMean( tzmax, zmax, T=2*3600.0 )

  #plt.plot(tzmin,zmin)
  #plt.plot(tzmax,zmax)
  #plt.show()
  # grid resolution
  dt = 18*60 #25*60
  t = np.arange( tzmin[0],tzmin[-1], dt )
  nZ = 25 #20
  T = np.tile(t,(nZ,1)).T
  zmin2 = interp1d( tzmin, zmin )(t)
  zmax2 = interp1d( tzmax, zmax )(t)
  lam = np.linspace(0,1,nZ)
  Z = (np.outer( lam, zmin2 ) + np.outer( 1-lam, zmax2 )).T
  #plt.pcolormesh(T,Z,np.ones_like(Z), facecolor=None, edgecolor='k')
  #plt.axis('equal')
  #plt.show()

  #uDC = interpolateOnGridDC( pc_u, T, Z )
  #tDC = interpolateOnGridDC( temp, T, Z )
  #sDC = interpolateOnGridDC( salt, T, Z )
  #u_int = interpolateOnGrid( pc_u, T, Z )
  #t_int = interpolateOnGrid( temp, T, Z )
  #s_int = interpolateOnGrid( salt, T, Z )

  #pc_u = runningMean( pc_u, T=60.0*2 )
  #temp = runningMean( temp, T=60.0*2 )
  #salt = runningMean( salt, T=60.0*2 )

  #print T.min(),T.max(),T.max()-T.min()
  #print Z.min(),Z.max(),Z.max()-Z.min()
  tScale = (T.max()-T.min())*1.0
  zScale = (Z.max()-Z.min())*2.0
  u_int = binOnGrid( pc_u, T, Z, tScale, zScale )
  t_int = binOnGrid( temp, T, Z, tScale, zScale )
  s_int = binOnGrid( salt, T, Z, tScale, zScale )

  # smooth the data before computing gradients
  radius = 3 #3
  u_int = blur_image( u_int, radius )
  t_int = blur_image( t_int, radius )
  s_int = blur_image( s_int, radius )
  # crop T,Z to match the smoothed shape
  off0 = (T.shape[0]-u_int.shape[0])/2
  off1 = (T.shape[1]-u_int.shape[1])/2
  T = T[off0:-off0,off1:-off1]
  Z = Z[off0:-off0,off1:-off1]

  ## plot interpolated u,t,s, rho data
  #fig = plt.figure()
  #ax = fig.add_subplot(4,1,1)
  #ax.pcolormesh(T,-Z,u_int,vmin=-2.0,vmax=2.0)
  #ax = fig.add_subplot(4,1,2)
  #ax.pcolormesh(T,-Z,t_int,vmin=12.0,vmax=13.5)
  #ax = fig.add_subplot(4,1,3)
  #ax.pcolormesh(T,-Z,s_int,vmin=1.0,vmax=31.0)
  #rho_int = equationOfStateJackett( s_int, t_int, np.zeros_like(s_int) )
  #print np.nanmin(rho_int), np.nanmax(rho_int)
  #ax = fig.add_subplot(4,1,4)
  #ax.pcolormesh(T,-Z,rho_int,vmin=1000,vmax=1024)
  #plt.show()

  # compute vertical gradient
  # cell height
  dz = np.diff(-Z,axis=1)
  # differences
  du = np.diff(u_int,axis=1)
  ds = np.diff(s_int,axis=1)
  dt = np.diff(t_int,axis=1)
  # vert gradient
  dudz = du/dz
  dsdz = ds/dz
  dtdz = dt/dz
  # values at cell center
  s_cnrt = s_int[:,:-1] + 0.5*ds
  t_cnrt = t_int[:,:-1] + 0.5*dt
  p_cntr = np.zeros_like(t_cnrt)
  Zcntr = -Z[:,:-1] + 0.5*dz
  Tcntr = T[:,:-1]
  # density gradient
  drhodz = equationOfStateJackettAlpha( s_cnrt, t_cnrt, p_cntr )*dtdz + equationOfStateJackettBeta( s_cnrt, t_cnrt, p_cntr )*dsdz
  # gradient Richardson number
  rho0 = 1000.0
  g = 9.81
  N2 = -g/rho0*drhodz
  S2 = dudz**2
  Ri = N2/S2

  # create profile dataContainers
  def makeProfDataContainer(trackDC,Z,T,C,var=None) :
    meta = trackDC.getMetaData()
    #meta['location'] = station
    #meta['instrument'] = 'model'
    #meta['variable'] = var
    meta['bracket'] = 'A'
    meta['dataType'] = 'profile'
    #meta.pop('instrument',None)
    meta.pop('institute',None)
    fieldNames = trackDC.fieldNames
    if var != None :
      meta['variable'] = var
      fieldNames = [ var ]
    nZ = C.shape[1]
    x = trackDC.x[0,0] if len(trackDC.x.shape) == 2 else trackDC.x
    y = trackDC.y[0,0] if len(trackDC.y.shape) == 2 else trackDC.y
    x = x*np.ones((nZ,))
    y = y*np.ones((nZ,))
    z = Z.T
    data = C.T[:,None,:]
    ta = timeArray(T[:,0],'epoch')
    dc = dataContainer('', ta, x,y,z, data, fieldNames,
                        coordSys='spcs',metaData=meta,acceptNaNs=True)
    return dc

  RiDC = makeProfDataContainer( salt, Zcntr, Tcntr, Ri, var='Ri' )
  S2DC = makeProfDataContainer( salt, Zcntr, Tcntr, S2, var='shearFreq' )
  N2DC = makeProfDataContainer( salt, Zcntr, Tcntr, N2, var='buoyFreq' )
  
  return salt, temp, hvel, pc_u, S2DC, N2DC, RiDC

def makeRiPlot( S2DC, N2DC, RiDC, xlim=None, ylim=None ) :
  # make stackPlot with all 3 panels
  from plotting.profilePlot import *
  dia = stackProfileTimeSeriesDC(xlabel=str(RiDC.time.getDatetime(0).year),ylim=ylim)
  dia.addPlot( 'S2', clabel='S2', unit='s-2',
              logScale=True,clim=[-4,-1.5],climIsLog=True)
  dia.addSample( 'S2', S2DC )
  dia.addPlot( 'N2', clabel='N2', unit='s-2',
              logScale=True,clim=[-6,-1.5],climIsLog=True)
  dia.addSample( 'N2', N2DC )
  #dia.addPlot( 'ri', clabel='Ri', unit='-',clim=[0,2.0] )#, clim=[0,1] ), clim=[0,0.5] )
  #dia.addSample( 'ri', RiDC, zorder=2, plotType='contour',
                #levels=[0.25],colors='w',linewidths=2,linestyles='solid' )
  #dia.addSample( 'ri', RiDC, zorder=3, plotType='contour',
                #levels=[0.25],colors='k',linewidths=2,linestyles='dashed' )
  #dia.addSample( 'ri', RiDC, zorder=0 )
  RiDC.data[RiDC.data<0] = 1e-12
  dia.addPlot( 'ri', clabel='Ri', unit='n/a',
              logScale=True,clim=[-2,2],climIsLog=True)#, clim=[0,1] ), clim=[0,0.5] )
  dia.addSample( 'ri', RiDC, zorder=2, plotType='contour',
                levels=[np.log10(0.25)],colors='w',linewidths=2.5,linestyles='solid' )
  dia.addSample( 'ri', RiDC, zorder=3, plotType='contour',
                levels=[np.log10(0.25)],colors='k',linewidths=2,linestyles='dashed' )
  dia.addSample( 'ri', RiDC, zorder=0, logScale=True )
  dia.showColorBar()
  dia.addTitle(' '.join([RiDC.getMetaData('tag'),RiDC.getMetaData('location').upper()]))
  if xlim :
    # force xaxis to the original data limits
    dia.plots['ri'].updateXAxis( xlim=xlim, numticks=8 )

  dateStr = RiDC.time.getDatetime(0).strftime('%Y-%m-%d')
  fn = '_'.join(['riNumber2',RiDC.getMetaData('tag'),RiDC.getMetaData('location'),dateStr])
  saveFigure( 'paperPlots', fn, 'png', verbose=True, dpi=100, bbox_tight=True )
  #plt.show()

  #print np.nanmin(dudz), np.nanmax(dudz)
  #fig = plt.figure()
  #ax = fig.add_subplot(3,1,1)
  #ax.pcolormesh(Tcntr,Zcntr,dudz**2,vmin=0.,vmax=0.05)
  #ax = fig.add_subplot(3,1,2)
  #print np.nanmin(-g/rho0*drhodz), np.nanmax(-g/rho0*drhodz)
  #ax.pcolormesh(Tcntr,Zcntr,-g/rho0*drhodz,vmin=-0.0,vmax=0.05)
  #ax = fig.add_subplot(3,1,3)
  #print np.nanmin(Ri), np.nanmax(Ri)
  #ax.pcolormesh(Tcntr,Zcntr,Ri,vmin=0.0,vmax=1)
  #plt.show()

#-------------------------------------------------------------------------------
# Main script
#-------------------------------------------------------------------------------

#tags = ['obs','db31-2012','run011']
#tags = ['run011']
tags = ['db31-2012']
#locs = ['oc1','oc2']
locs = ['oc1']

for tag in tags :
  for loc in locs :
    salt, temp, hvel, pc_u, S2, N2, Ri = computeRi( tag, loc )

    # save to disk
    #dtm.saveDataContainerInTree( pc_u, rule=oldRule, overwrite=True )
    #dtm.saveDataContainerInTree( Ri, rule=oldRule, overwrite=True )
    #dtm.saveDataContainerInTree( S2, rule=oldRule, overwrite=True )
    #dtm.saveDataContainerInTree( N2, rule=oldRule, overwrite=True )

    xlim=[convertEpochToPlotTime(salt.time[0]),
          convertEpochToPlotTime(salt.time[-1])]

    makeRiPlot( S2, N2, Ri, xlim=xlim, ylim=[-19,0] )