"""
Methods for manipulating winched profiler data.

Tuomas Karna 2013-01-17
"""
import numpy as np
from scipy.interpolate import interp1d
from crane.data import dataContainer
from crane.data import timeArray

def generateSat01ProfilerModData( obsWProfilerDC, modProfileDC ) :

  # merge time ranges
  tmin = max( obsWProfilerDC.time.getDatetime(0), modProfileDC.time.getDatetime(0) )
  tmax = min( obsWProfilerDC.time.getDatetime(-1), modProfileDC.time.getDatetime(-1) )

  o = obsWProfilerDC.timeWindow( tmin,tmax,includeEnd=True )
  m = modProfileDC.timeWindow( tmin,tmax,includeEnd=True )

  # get all distinct time stamps in obs (it's a track)
  to = o.time.asEpoch().array
  tArr, ixTo = np.unique( to, return_index=True )
  ix = np.zeros((len(ixTo),2),dtype=int)
  ix[:,0] = ixTo
  ix[:-1,1] = ixTo[1:]
  ix[-1,1] = len(to)

  # align times
  tix = m.time.getAlignedTimeIndices( timeArray.timeArray(tArr,'epoch') )
  tArr = tArr[tix]
  ta = timeArray.timeArray( tArr, 'epoch' )
  m = m.interpolateInTime( ta, acceptNaNs=True )

  if not np.array_equal( tArr, m.time.array ) :
    print tArr.shape, m.time.array.shape
    print tArr.min(), tArr.max()
    print m.time.array.min(), m.time.array.max()
    raise Exception( 'the time stamps of observation and model do not agree: '+str(ta.getDatetime(0))+', '+str(m.time.getDatetime(0)) )
  if not m.zDependsOnTime or m.z.shape[1] != len(tArr) :
    raise Exception( 'model z coordinate does not have correct time dimension: '+str(m.z.shape))

  # for all time stamps, process each vertical
  znew = []
  vnew = []
  tnew = []
  for i,t in enumerate( tArr ) :
    # get observation z coordinates (depth)
    zobs = o.z[0,ix[i,0]:ix[i,1]]
    # get model z coords
    zmod = m.z[:,i]
    # convert to depth
    zmod = zmod.max()-zmod # surf_val - z_coord
    zsortix = np.argsort(zmod)
    zmod = zmod[zsortix]
    # get model values
    vmod = m.data[:,0,i][zsortix]
    # discard zobs values that are out of range
    goodIx = np.logical_and( zobs <= zmod.max(), zobs >= zmod.min() )
    zobs = zobs[goodIx]
    # do the interpolation
    vint = interp1d( zmod, vmod )( zobs )
    znew.append( zobs )
    vnew.append( vint )
    tnew.append( t*np.ones_like(zobs) )
  znew = np.concatenate( tuple(znew) )
  vnew = np.concatenate( tuple(vnew) )
  tnew = np.concatenate( tuple(tnew) )
  ta = timeArray.timeArray( tnew, 'epoch',acceptDuplicates=True )

  # create dataContainer
  z = znew[None,:]
  data = vnew[None,None,:]
  meta = m.getMetaData()
  meta['dataType'] = 'profiler'
  #meta.pop( 'msldepth', None ) # TODO
  dc = dataContainer.dataContainer('', ta, m.x[0],m.y[0],z,data,o.fieldNames,
                     m.coordSys,meta,acceptNaNs=True)
  return dc

def generateSat01ProfilerObsData( depDC, varDC ) :
  """Generate winched profiler data for Saturn01.
  Observation data is binned in time and depth to produce vertical columns of
  data at regular intervals.
  
  Args:
  depDC -- (dataContainer) time series of profiler depth [m below surf]
  varDC -- (dataContainer) time series of measured variable

  Returns:
  dc -- (dataContainer) track data with correct depth

  """
  tReso = 15*60 # target time step (15min to compare with model)
  zReso = 0.25 # target vertical resolution
  # interpolate pressure data on instrument time
  depDC,varDC = depDC.alignTimes( varDC )
  # bin data in bot time and depth
  # round start/end time to nearest sharp 15min
  t = depDC.time.asEpoch()
  z = depDC.data.squeeze()
  varArr = varDC.data.squeeze()
  ts = t[0]
  te = t[-1]
  binDt = tReso
  ts = round(ts/binDt)*binDt
  te = round(te/binDt)*binDt
  # divide t to tReso bins
  tBins = np.arange(ts-binDt/2,te+binDt*3/2,binDt) # bin bnds
  ixBin = np.digitize( t, tBins )-1 # bin index for each time stamp
  # generate list of time stamps for each bin TODO slow, optimize
  binMembers = [ [] for _ in tBins[1:] ]
  for i,ibin in enumerate(ixBin) :
    binMembers[ibin].append(i)
  #tmp, bix = np.unique( ixBin, return_index=True )
  #ix = np.zeros((len(bix),2),dtype=int)
  #ix[:,0] = bix
  #ix[:-1,1] = bix[1:]
  #ix[-1,1] = len(ixBin)
  
  tNew = []
  zNew = []
  varNew = []
  for i in range(ixBin[0],ixBin[-1]) :
    tCenter = (tBins[i]+tBins[i+1])/2
    # find min/max z
    iii = binMembers[i]
    if len(iii) == 0 :
      continue
    zBin = z[iii]
    tBin = t[iii]
    varBin = varArr[iii]
    # generate z grid at bin mean time
    zMin = zBin.min()
    zMax = zBin.max()
    zGrid = np.arange( zMin, zMax+zReso,zReso )
    tGrid = np.ones_like(zGrid)*tCenter
    zBins = np.arange( zMin-zReso/2, zMax+3/2*zReso,zReso )
    histW,b = np.histogram( zBin, zBins, weights=varBin )
    hist,b = np.histogram( zBin, zBins )
    goodIx = hist > 0
    varGrid = histW[goodIx]/hist[goodIx]
    tGrid = tGrid[goodIx]
    zGrid = zGrid[goodIx]
    # store data
    tNew.append( tGrid )
    zNew.append( zGrid )
    varNew.append( varGrid )

  if len(tNew) == 0 :
    print 'Could not generate obs profiler data'
    return None
  tNew = np.concatenate( tuple(tNew) )
  zNew = np.concatenate( tuple(zNew) )
  varNew = np.concatenate( tuple(varNew) )
  # create new dataContainer
  tNew = timeArray.timeArray(tNew,'epoch',acceptDuplicates=True)
  z = np.reshape(zNew,(1,-1))
  data = np.reshape(varNew,(1,1,-1))
  var = varDC.fieldNames[0]
  meta = varDC.getMetaData()
  #meta.pop( 'msldepth', None ) # TODO
  meta['dataType'] = 'profiler'
  dc = dataContainer.dataContainer('',tNew,
                        varDC.x,varDC.y,z,
                        data,[var],varDC.coordSys,metaData=meta)
  return dc