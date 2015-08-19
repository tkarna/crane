"""
A script for computing gradient Richardson number from model profiles.
"""

import numpy as np
from crane.data import dataContainer
from crane.data import timeArray
from data.stationCollection import *
import data.dirTreeManager as dtm
from data.pca import *
from data.eqState import *
from scipy.spatial import cKDTree
from plotting.plotBase import convertEpochToPlotTime

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------

st = datetime.datetime(2012,10,26)
et = datetime.datetime(2012,11,3)

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def computeRi( tag, loc ) :

  hvel = dtm.getDataContainer( dataType='profile', tag=tag, location=loc,variable='hvel', startTime=st, endTime=et )
  salt = dtm.getDataContainer( dataType='profile', tag=tag, location=loc,variable='salt', startTime=st, endTime=et )
  temp = dtm.getDataContainer( dataType='profile', tag=tag, location=loc,variable='temp', startTime=st, endTime=et )

  def doPCAprof( dc ) :
    A = np.swapaxes( dc.data[:,:2,:], 1,2 )
    A = np.reshape( A, (-1,dc.data.shape[1]) )
    print A.mean(axis=0)
    print A.shape
    p = PCA(A,fraction=0.9)
    print "npc:", p.npc
    print "% variance:", p.sumvariance * 100
    B = p.pc()
    print B.shape
    dc2 = dc.copy()
    dc2.data = np.reshape(B, (dc.data.shape[0],1,dc.data.shape[2]) )
    dc2.fieldNames = dc.fieldNames[:dc2.data.shape[1]]
    dc2.setMetaData('variable','pcavel')
    return dc2

  pc_u = doPCAprof( hvel )

  # rotate velocity so that flood is positive
  ix = pc_u.z[:,0] < 10
  mvel = pc_u.data[ix,0,:].mean()
  print 'mean surf vel',mvel
  if mvel > 0 :
    print 'flipping'
    pc_u.data *= -1

  # compute vertical gradients of u,temp,salt
  dz = np.diff(pc_u.z,axis=0)
  print pc_u.z.shape, dz.shape
  du = np.diff(pc_u.data[:,0,:],axis=0)
  ds = np.diff(salt.data[:,0,:],axis=0)
  dt = np.diff(temp.data[:,0,:],axis=0)
  # vert gradient
  dudz = du/dz
  dsdz = ds/dz
  dtdz = dt/dz
  # values at cell center
  s_cnrt = salt.data[:-1,0,:] + 0.5*ds
  t_cnrt = temp.data[:-1,0,:] + 0.5*dt
  p_cntr = np.zeros_like(t_cnrt)
  z_cntr = salt.z[:-1,:] + 0.5*dz
  T_cntr = np.tile(salt.time.array,(z_cntr.shape[0],1))
  print z_cntr.shape, T_cntr.shape
  # density gradient
  drhodz = equationOfStateJackettAlpha( s_cnrt, t_cnrt, p_cntr )*dtdz + equationOfStateJackettBeta( s_cnrt, t_cnrt, p_cntr )*dsdz
  # gradient Richardson number
  rho0 = 1000.0
  g = 9.81
  N2 = -g/rho0*drhodz
  S2 = dudz**2
  Ri = N2/S2

  #print np.nanmin(dudz), np.nanmax(dudz)
  #plt.figure()
  #plt.pcolormesh(T_cntr,z_cntr,dudz**2,vmin=0.,vmax=0.05)
  #plt.figure()
  #print np.nanmin(-g/rho0*drhodz), np.nanmax(-g/rho0*drhodz)
  #plt.pcolormesh(T_cntr,z_cntr,-g/rho0*drhodz,vmin=-0.0,vmax=0.05)
  #plt.figure()
  #print np.nanmin(Ri), np.nanmax(Ri)
  #plt.pcolormesh(T_cntr,z_cntr,Ri,vmin=0.0,vmax=1)
  #plt.colorbar()

  # create profile dataContainers
  def makeProfDataContainer(profDC,Z,T,C,var=None) :
    meta = profDC.getMetaData()
    #meta['location'] = station
    #meta['instrument'] = 'model'
    #meta['variable'] = var
    meta['bracket'] = 'A'
    meta['dataType'] = 'profile'
    #meta.pop('instrument',None)
    meta.pop('institute',None)
    fieldNames = profDC.fieldNames
    if var != None :
      meta['variable'] = var
      fieldNames = [ var ]
    nZ = C.shape[0]
    x = profDC.x[0]
    y = profDC.y[0]
    x = x*np.ones((nZ,))
    y = y*np.ones((nZ,))
    z = Z
    data = C[:,None,:]
    ta = profDC.time #timeArray(T[:,0],'epoch')
    dc = dataContainer('', ta, x,y,z, data, fieldNames,
                        coordSys='spcs',metaData=meta,acceptNaNs=True)
    return dc

  RiDC = makeProfDataContainer( salt, z_cntr, T_cntr, Ri, var='Ri' )
  S2DC = makeProfDataContainer( salt, z_cntr, T_cntr, S2, var='shearFreq' )
  N2DC = makeProfDataContainer( salt, z_cntr, T_cntr, N2, var='buoyFreq' )

  return salt, temp, hvel, pc_u, S2DC, N2DC, RiDC
  
def profileZToDepth( dc ) :
  dc2 = dc.copy()
  dc2.z = dc2.z-np.tile(dc2.z.max(axis=0),(dc2.z.shape[0],1))
  return dc2

def makeRiPlot( S2DC, N2DC, RiDC, xlim=None, ylim=None ) :
  # make stackPlot with all 3 panels
  from plotting.profilePlot import *
  dia = stackProfileTimeSeriesDC(xlabel=str(RiDC.time.getDatetime(0).year),ylim=ylim)
  dia.addPlot( 'S2', clabel='S2', unit='s-2',
              logScale=True,clim=[-4,-1.5],climIsLog=True)

  dia.addSample( 'S2', profileZToDepth(S2DC) )
  dia.addPlot( 'N2', clabel='N2', unit='s-2',
              logScale=True,clim=[-6,-1.5],climIsLog=True)
  dia.addSample( 'N2', profileZToDepth(N2DC) )
  #dia.addPlot( 'ri', clabel='Ri', unit='-',clim=[0,2.0] )#, clim=[0,1] ), clim=[0,0.5] )
  #dia.addSample( 'ri', RiDC, zorder=2, plotType='contour',
                #levels=[0.25],colors='w',linewidths=2,linestyles='solid' )
  #dia.addSample( 'ri', RiDC, zorder=3, plotType='contour',
                #levels=[0.25],colors='k',linewidths=2,linestyles='dashed' )
  #dia.addSample( 'ri', RiDC, zorder=0 )
  RiDC.data[RiDC.data<0] = 1e-12
  dia.addPlot( 'ri', clabel='Ri', unit='n/a',
              logScale=True,clim=[-2,2],climIsLog=True)#, clim=[0,1] ), clim=[0,0.5] )
  dia.addSample( 'ri', profileZToDepth(RiDC), zorder=2, plotType='contour',
                levels=[np.log10(0.25)],colors='w',linewidths=2.5,linestyles='solid' )
  dia.addSample( 'ri', profileZToDepth(RiDC), zorder=3, plotType='contour',
                levels=[np.log10(0.25)],colors='k',linewidths=2,linestyles='dashed' )
  dia.addSample( 'ri', profileZToDepth(RiDC), zorder=0, logScale=True )
  dia.showColorBar()
  dia.addTitle(' '.join([RiDC.getMetaData('tag'),'profile',
                         RiDC.getMetaData('location').upper()]))
  if xlim :
    # force xaxis to the original data limits
    dia.plots['ri'].updateXAxis( xlim=xlim, numticks=8 )

  dateStr = RiDC.time.getDatetime(0).strftime('%Y-%m-%d')
  fn = '_'.join(['riNumberProf',RiDC.getMetaData('tag'),RiDC.getMetaData('location'),dateStr])
  saveFigure( 'paperPlots', fn, 'png', verbose=True, dpi=100, bbox_tight=True )

#-------------------------------------------------------------------------------
# Main script
#-------------------------------------------------------------------------------

tags = ['db31-2012']
#locs = ['oc1','oc2']
locs = ['oc1']

for tag in tags :
  for loc in locs :
    salt, temp, hvel, pc_u, S2, N2, Ri = computeRi( tag, loc )

    ## save to disk
    #dtm.saveDataContainerInTree( pc_u, rule=oldRule, overwrite=True )
    #dtm.saveDataContainerInTree( Ri, rule=oldRule, overwrite=True )
    #dtm.saveDataContainerInTree( S2, rule=oldRule, overwrite=True )
    #dtm.saveDataContainerInTree( N2, rule=oldRule, overwrite=True )

    sTrack = dtm.getDataContainer( dataType='track', tag=tag, location=loc,variable='salt', startTime=st, endTime=et )
    xlim=[convertEpochToPlotTime(sTrack.time[0]),
          convertEpochToPlotTime(sTrack.time[-1])]

    makeRiPlot( S2, N2, Ri, xlim=xlim, ylim=[-19,0] )