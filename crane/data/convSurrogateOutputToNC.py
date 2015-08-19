#!/usr/bin/env python
"""
A script for extracting data from Matlab surrogate model profile time series,
and storing it in netCDF format.

Tuomas Karna 2013-11-06

"""
import os
import sys
import traceback
import numpy as np
from crane.data import timeArray
from crane.data import dataContainer
import glob
from scipy.io.matlab import loadmat
import datetime
from data.coordSys import *
import data.dirTreeManager as dtm
from scipy.interpolate import interp1d

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def datetime2matlabdn(dt):
  """Convert datetime to matlab datenum format"""
  mdn = dt + datetime.timedelta(days = 366)
  frac_seconds = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
  frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
  return mdn.toordinal() + frac_seconds + frac_microseconds

datenum_epoch_offset = datetime2matlabdn(datetime.datetime(1970,1,1))
def datenumUTCToEpoch( dn ) :
  """Covert Matlab datenum (in UTC time zone) to epoch format (UTC)"""
  return (dn-datenum_epoch_offset)*86400.0

def datenumPSTToEpoch( dn ) :
  """Covert Matlab datenum (in PST time zone) to epoch format (UTC)"""
  return (dn-datenum_epoch_offset)*86400.0 + 8*3600.0

def read_mat_sil_files( path, fn ) :
  """Reads SIL time series from matlab file and converts data
  to python format"""
  f = os.path.join(path,fn)
  print 'reading',f
  d = loadmat(f)
  t = d['t'].flatten() # (1,nTime)
  data = d['data'].flatten() # (1,nTime)

  # convert time from Matlab datenum (in PST) to epoch (UTC)
  time = datenumPSTToEpoch(t)
  # round to nearest minute
  time = np.round(time/60.)*60.
  print 'Loaded data for range:\n ',str(timeArray.epochToDatetime(time[0])),' -> ',  str(timeArray.epochToDatetime(time[-1]))
  return time,data

def genSILDC( path,fn,runTag ) :
  """Reads matlab SIL time series and returns a dataContainer"""
  time,vals = read_mat_sil_files( path,fn )
  var = fn.split('.')[0]
  meta={}
  meta['description'] = 'surrogate model'
  meta['instrument'] = 'surrogate'
  #meta['msldepth'] = '0'
  meta['dataType'] = 'sil'
  meta['location'] = 'mainChannel'
  meta['bracket'] = 'A' # fixed z coordinate
  meta['variable'] = var
  meta['tag'] = runTag
  goodIx = np.logical_and( np.isfinite(time), np.any(np.isfinite(vals),axis=0) )
  ta = timeArray.timeArray(time[goodIx],'epoch')
  data = vals[goodIx][None,None,:] / 1000.0
  nReg = 1
  xx = np.zeros((nReg,))
  yy = np.zeros((nReg,))
  zz = np.zeros((nReg,))
  fieldNames = [var]
  
  dc = dataContainer.dataContainer('', ta, xx,yy,zz, data, fieldNames, coordSys='spcs',metaData=meta,acceptNaNs=True)
  return dc

def read_mat_sho_files( path, fn ) :
  """Reads SHO time series from matlab file and converts data
  to python format"""
  f = os.path.join(path,fn)
  print 'reading',f
  d = loadmat(f)
  t = d['t'].flatten() # (1,nTime)
  data = d['data'].T # (nTime,nReg)
  indData = d['shoIndv'].T # (nTime,nCriteria,nReg)

  # convert time from Matlab datenum (in PST) to epoch (UTC)
  time = datenumPSTToEpoch(t)
  # round to nearest minute
  time = np.round(time/60.)*60.
  print 'Loaded data for range:\n ',str(timeArray.epochToDatetime(time[0])),' -> ',  str(timeArray.epochToDatetime(time[-1]))
  
  return time,data,indData
  
def genSHODC( path,fn,runTag ) :
  """Reads matlab SHO time series and returns a dataContainer"""
  time,vals,indData = read_mat_sho_files( path,fn )
  meta={}
  meta['description'] = 'surrogate model'
  meta['instrument'] = 'surrogate'
  #meta['msldepth'] = '0'
  meta['dataType'] = 'sho'
  meta['location'] = 'cr'
  meta['bracket'] = 'A' # fixed z coordinate
  meta['variable'] = 'shov1'
  meta['tag'] = runTag
  goodIx = np.logical_and( np.isfinite(time), np.any(np.isfinite(vals),axis=0) )
  ta = timeArray.timeArray(time[goodIx],'epoch')
  data = vals[:,goodIx][:,None,:]
  indd = indData[:,:,goodIx]
  data = np.concatenate( (data,indd), axis=1 )
  nReg = vals.shape[0]
  xx = np.zeros((nReg,))
  yy = np.zeros((nReg,))
  zz = np.zeros((nReg,))
  fieldNames = ['sho','sho_t','sho_s','sho_d','sho_v']
  
  dc = dataContainer.dataContainer('', ta, xx,yy,zz, data, fieldNames, coordSys='spcs',metaData=meta,acceptNaNs=True)
  return dc

def read_mat_profile_files( path,loc,var,dataSetName='test',dataSetType='ms' ) :
  """Reads generic time series from matlab file and converts data
  to python format"""
  varToChar = {'salt':'s', 'elev':'e', 'temp':'t', 'u':'u', 'v':'v'}
  pattern = os.path.join(path,dataSetName+'.'+dataSetType+'.'+varToChar[var]+'.'+loc+'.mat')
  fList = sorted(glob.glob(pattern))
  if not fList :
    raise Exception('File not found: '+pattern)
  f = fList[0]
  print 'Reading',f
  d = loadmat(f)
  t = d['t'].flatten() # (1,nTime)
  z = d['z'] # (nVert,nTime)
  data = d['data'] # (nVert,nTime)
  # convert time from Matlab datenum (in PST) to epoch (UTC)
  time = datenumPSTToEpoch(t)
  # round to nearest minute
  time = np.round(time/60.)*60.
  print '  Loaded data range: ',str(timeArray.epochToDatetime(time[0])),' -> ',  str(timeArray.epochToDatetime(time[-1]))
  return time,z,data

def genProfileDC(dataDir,runTag,location,var,x,y,
                 dataSetName='test', dataSetType='ms') :
  """Reads matlab time series and returns a dataContainer with correct
  station location"""
  if var=='hvel':
    dc = genProfileDC(dataDir,runTag,location,'u',x,y,dataSetName,dataSetType)
    dc2 = genProfileDC(dataDir,runTag,location,'v',x,y,dataSetName,dataSetType)
    dc.mergeFields( dc2 )
    dc.fieldNames = ['u','v']
    dc.setMetaData( 'variable',var )
    return dc
  fieldNames = [var]
  time,z,vals = read_mat_profile_files(dataDir, location, var,
                                       dataSetName=dataSetName,
                                       dataSetType=dataSetType)
  meta={}
  meta['description'] = 'surrogate model'
  meta['instrument'] = 'surrogate'
  if var == 'elev' :
    meta['msldepth'] = '0'
    meta['dataType'] = 'timeseries'
  else :
    meta['dataType'] = 'profile'
  meta['location'] = location
  meta['bracket'] = 'A' # fixed z coordinate
  meta['variable'] = var
  meta['tag'] = runTag
  goodIx = np.logical_and( np.isfinite(time), np.any(np.isfinite(vals),axis=0) )
  if var != 'elev' :
    goodIx = np.logical_and( goodIx, np.any(np.isfinite(z),axis=0) )
  #print goodIx.shape
  ta = timeArray.timeArray(time[goodIx],'epoch')
  data = vals[:,goodIx][:,None,:]
  zz = z[:,goodIx]
  if var == 'elev' :
    zz = np.array([0])
  else :
    goodZ = np.isfinite( vals ).all(axis=1)
    zz = zz[goodZ,:]
    data = data[goodZ,:]
  xx = np.tile( x, (zz.shape[0],) )
  yy = np.tile( y, (zz.shape[0],) )
  dc = dataContainer.dataContainer('', ta, xx,yy,zz, data, fieldNames, coordSys='spcs',metaData=meta,acceptNaNs=True)
  return dc

def interpolateInVertical(Z,V,z=None,k=None,zRelToSurf=False):
  """
  Interpolates vertical profile at given depth
  
  Parameters
  ----------
  Z : ndarray (nZ,nTime)
    z coordinates for each vertical profile. If nTime>1 each profile will be
    interpolated separately.
  V : ndarray (nZ,nTime)
    values corresponding to each z point
  z : float, optional
    z coordinate to interpolate to
  k : int
    optionally, take k-th nodal value from bottom.
    k=1 stands for bottom, k=-1 stands for surface
  zRelToSurf : bool
    If True z coordinate is taken depth below free surface instead of
    static z coordinate (increasing upwards)
  """
  if k==None and z==None :
    raise Exception('either k or z must be given')
  # check that Z[0,:] is the bottom
  ZZ = Z.copy()
  VV = V.copy()
  if ZZ[-1,0] - ZZ[0,0] < 0 :
    # nope, flip first indices
    ZZ = ZZ[::-1,:]
    VV = VV[::-1,:]
  # ensure that nan masks match
  ixnan = ~np.isfinite(ZZ)
  ixnan = np.logical_or(ixnan,~np.isfinite(VV))
  nanPadded = ixnan.any()
  if nanPadded :
    ZZ[ixnan] = np.nan
    VV[ixnan] = np.nan
  nZ, nTime = ZZ.shape
  v = np.zeros((nTime,))
  if z != None :
    # interpolate in vertical
    z0 = z
    z_bot = np.nanmin(ZZ,axis=0)
    z_sur = np.nanmax(ZZ,axis=0)
    if zRelToSurf :
      z_target = z_sur - z
    else :
      z_target = z*np.ones((nTime))
    # check bounds
    z_target = np.minimum( np.maximum( z_target, z_bot ), z_sur )
    if nanPadded :
      for iTime in xrange(nTime) :
        ix = np.isfinite(ZZ[:,iTime])
        v[iTime] = interp1d( ZZ[ix,iTime], VV[ix,iTime], kind='linear', copy=False ) ( z_target[iTime] )
    else :
      for iTime in xrange(nTime) :
        v[iTime] = interp1d( ZZ[:,iTime], VV[:,iTime], kind='linear', copy=False ) ( z_target[iTime] )
  if k != None :
    # bottom: k=1 kk=0, surface: k=-1 kk=len(z)-1
    kk = k-1 if k>0 else nZ+k
    v = VV[kk,:]
    z_target = ZZ[kk,:]
  return v,z_target

def getTSFromProfile(runTag, location, var, z=None, k=None, 
                     zRelToSurf=False, rule=None) :
  """Extracts a time series from vertical profile at given depth
  
  Parameters
  ----------
  z : float
      z coordinate (<0 below datum) or
      depth (>0 below surface) if zRelToSurf=True
  k : int
      model vertical level k=1 is bottom, k=-1 is surface
  zRelToSurf : bool
    If True z coordinate is taken depth below free surface instead of
    static z coordinate (increasing upwards)

  """
  if var == 'elev' :
    return dtm.getDataContainer(location=location, tag=runTag,
                                dataType='timeseries', variable=var, rule=rule)

  prof = dtm.getDataContainer(location=location, tag=runTag,
                              dataType='profile', variable=var, rule=rule)
  if k==None and z==None :
    raise Exception('either k or z must be given')
  v = []
  nComp = prof.data.shape[1]
  for iComp in range(nComp) :
    v0,z_target = interpolateInVertical(prof.z,prof.data[:,iComp,:],z,k)
    #z = z_target.mean()
    v.append(v0)
  v = np.vstack(tuple(v))
  # make dataContainer
  meta=prof.getMetaData()
  meta['dataType'] = 'timeseries'
  meta['bracket'] = 'F' if zRelToSurf else 'A'
  meta['msldepth'] = str(int(round(abs(z)*100)))
  data = v[None,:,:]
  x = prof.x[0]
  y = prof.y[0]
  dc = dataContainer.dataContainer('', prof.time, x,y,z, data, prof.fieldNames, coordSys=prof.coordSys,metaData=meta)
  return dc

def read_mat_plume_file( path,var,saltThreshold ) :
  """Reads plume metrics from matlab file and converts data
  to python format"""
  varToChar = {'plume_area':'parea','plume_center':'pcenter',
               'plume_thickness':'pthicknes','plume_volume':'pvolume'}
  f = os.path.join(path,varToChar[var]+'_ms_'+str(int(saltThreshold))+'.mat')
  if not os.path.isfile(f) :
    raise IOError('file not found: '+f)
  print 'Reading',f
  d = loadmat(f)
  t = d['t'].flatten() # (1,nTime)
  data = d['data'] # (nVert,nTime)
  # convert time from Matlab datenum (in PST) to epoch (UTC)
  time = datenumPSTToEpoch(t)
  # round to nearest minute
  time = np.round(time/60.)*60.
  print '  Loaded data range: ',str(timeArray.epochToDatetime(time[0])),' -> ',  str(timeArray.epochToDatetime(time[-1]))
  return time,data

def getPlumeDC(runTag,dataDir,saltThreshold):
  """Reads plume metrics for given salinity threshold and returns a
  dataContainer.
  """
  varList = ['plume_area','plume_center','plume_thickness','plume_volume']
  values = []
  for var in varList :
    t,d = read_mat_plume_file(dataDir,var,saltThreshold)
    values.append(d)
  values = np.hstack(tuple(values))
  goodIx = np.isfinite(np.sum(values,axis=1))
  time = t[goodIx]
  values = values[goodIx,:]

  varNames = ['plume_area','plume_center',
              'plume_volume','plume_thickness']
  fieldIndices = [[0],[1,2], [3],[4]]
  fieldNames = {'plume_center':['plume_center_x','plume_center_y']}
  dcs=[]
  for i,var in enumerate(varNames):
    sthSuffix = '_{0:d}'.format(int(saltThreshold))
    data = np.swapaxes(values[:,fieldIndices[i]],0,1)[None,:,:] # (1,nStats,nTime)
    ta = timeArray.timeArray( time, 'epoch' )
    meta = {}
    meta['tag'] = runTag
    meta['location'] = 'plume'
    meta['instrument'] = 'surrogate'
    meta['variable'] = var+sthSuffix
    meta['dataType'] = 'plumemetrics'
    meta['saltThreshold'] = str(saltThreshold)
    x = y = z = 0
    fNames = [ fn+sthSuffix for fn in fieldNames.get(var,[var]) ]
    dc = dataContainer.dataContainer('', ta, x,y,z, data, fNames,
                        coordSys='spcs',metaData=meta)
    dcs.append(dc)
  return dcs

def gatherPlumeMetrics(runTag,baseDir,startTime,endTime):
  """Read plume metrics for all salt thresholds found."""
  dataDir = os.path.join(baseDir,'metrics')
  dcs = []
  pattern = os.path.join(dataDir,'parea'+'_ms_'+'*'+'.mat')
  saltThresholds = set()
  for f in glob.glob(pattern) :
    threshold = float( f.split('_')[-1].replace('.mat','') )
    saltThresholds.add( threshold )
  for threshold in saltThresholds :
    try :
      dcs2 = getPlumeDC(runTag,dataDir,threshold)
      for dc in dcs2 :
        dc = dc.timeWindow(startTime,endTime)
      dcs.extend(dcs2)
    except Exception as e :
      print 'Reading plume metrics failed:',threshold
      print e
  if not dcs :
    print 'no plume data found'
  return dcs

def gatherProfiles(runTag,baseDir,varList,stationsToExtract,startTime,endTime,
                   dataSetName='test',dataSetType='ms') :
  """Reads all profile time series data from matlab files. Returns list of
  dataContainers."""
  dataDir = os.path.join(baseDir,'stations')
  if not os.path.isdir(dataDir) :
    raise IOError('directory not found: '+dataDir)
  dcs = []
  for var in varList :
    varStations = uniqueList( v[0] for v in
                             stationsToExtract.getTuples(variable=var) )
    for loc in varStations :
      try :
        x = stationsToExtract.getX(loc)
        y = stationsToExtract.getY(loc)
        dc = genProfileDC(dataDir,runTag,loc,var,x,y,dataSetName,dataSetType)
        dc = dc.timeWindow(startTime,endTime)
        dcs.append(dc)
      except Exception as e :
        print 'Reading matlab data failed:',loc,var
        traceback.print_exc(file=sys.stdout)
        print '  ',e
  return dcs

def gatherSIL(runTag,baseDir,startTime,endTime) :
  """Reads all salt intrusion length time series for matlab files. Returns a
  list of dataContainers."""
  dataDir = os.path.join(baseDir,'metrics')
  if not os.path.isdir(dataDir) :
    raise IOError('directory not found: '+dataDir)
  dcs = []
  for i in range(1,6) :
    dc = genSILDC( dataDir,'sil_'+str(i)+'.mat',runTag )
    dc = dc.timeWindow(startTime,endTime)
    dcs.append(dc)
  return dcs

def extractTimeSeriesFromProfiles(runTag,stationsToExtract,
                                  rule=dtm.defaultTreeRule()):
  """Extracts time series at a certain depth from vertical profiles"""
  dcs = []
  for key in stationsToExtract.getTuples() :
    loc,x,y,z,zType,var = key
    if var == 'elev' : continue # nothing to be done for elevations
    print '* extracting',loc,z,zType,var
    try :
      dc = getTSFromProfile(runTag, loc, var, z=z,
                            zRelToSurf=zType=='depth',rule=rule)
      dcs.append(dc)
    except Exception as e :
      #traceback.print_exc(file=sys.stdout)
      print ' ',e
  return dcs

#-----------------------------------------------------------------------------
# Main routine
#-----------------------------------------------------------------------------
from files.csvStationFile import csvStationFile,csvStationFileWithDepth
from data.collection import uniqueList

def gatherAndExtract(runTag,dataDir,startTime,endTime,extractFile,
                     dataSetName='test',dataSetType='ms') :
  stationsToExtract = csvStationFileWithDepth()
  stationsToExtract.readFromFile(extractFile)

  rule = dtm.defaultTreeRule()

  # read profile data and store in dataContainer format
  try :
    dcs = gatherProfiles(runTag,dataDir, ['elev','salt','temp','hvel'],
                         stationsToExtract, startTime,endTime,
                         dataSetName,dataSetType)
    dtm.saveDataContainerInTree(dcs, path='', rule=rule, dtype=np.float32,
                                overwrite=True, compress=True)
  except Exception as e:
    print 'Could not read station profiles'
    print e

  # SIL
  try :
    dcs = gatherSIL(runTag,dataDir,startTime,endTime)
    dtm.saveDataContainerInTree(dcs, path='', rule=rule, dtype=np.float32,
                                overwrite=True)
  except Exception as e:
    print 'Could not read SIL'
    print e

  # plume metrics
  try :
    dcs = gatherPlumeMetrics(runTag,dataDir,startTime,endTime)
    dtm.saveDataContainerInTree(dcs, path='', rule=rule, dtype=np.float32,
                                overwrite=True)
  except Exception as e:
    print 'Could not read plume metrics'
    print e

  # Extract time series at known (x,y,z) locations and store
  dcs = extractTimeSeriesFromProfiles(runTag,stationsToExtract)
  dtm.saveDataContainerInTree(dcs, path='', rule=rule, dtype=np.float32,
                              overwrite=True)

#-------------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------------
def parseCommandLine() :
  
  from optparse import OptionParser

  parser = OptionParser()
  parser.add_option('-r', '--runTag', action='store', type='string',
                      dest='runTag', help='Run tag, used in directory tree and\
                      as a label in post-processing.')
  parser.add_option('-d', '--dataDirectory', action='store', type='string',
                      dest='dataDir',
                      help='directory where surrogate subdirectories\
                      \'metrics\' and \'stations\' are located')
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
  parser.add_option('-t', '--csvStationFile', action='store', type='string',
                      dest='csvStationFile',
                      help='file that defines station coordinates and\
                      variables for time series extraction')
  parser.add_option('', '--dataSetType', action='store', type='string',
                      dest='dataSetType',default='ms',
                      help='data set to load: \'ms\' for surrogate\
                      , \'svd\' for EOFs (default %default)')
  parser.add_option('', '--dataSetName', action='store', type='string',
                      dest='dataSetName',default='test',
                      help='name of the data set to load: \'test\', \'valid\',... \
                      (default %default)')

  (options, args) = parser.parse_args()

  runTag        = options.runTag
  dataDir       = options.dataDir
  startStr      = options.startStr
  endStr        = options.endStr
  csvStationFile = options.csvStationFile
  dataSetType = options.dataSetType
  dataSetName = options.dataSetName

  if not runTag :
    parser.print_help()
    parser.error('runTag undefined')
  if not dataDir :
    parser.print_help()
    parser.error('dataDir undefined')
  if not startStr :
    parser.print_help()
    parser.error('startStr undefined')
  if not endStr :
    parser.print_help()
    parser.error('endStr undefined')
  if not csvStationFile :
    parser.print_help()
    parser.error('csvStationFile undefined')

  startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
  endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')

  print 'Parsed options:'
  print ' - runTag',runTag
  print ' - dataDir',dataDir
  print ' - surrogate data set:',dataSetName,dataSetType
  print ' - time range:',str(startTime),'->', str(endTime)
  print ' - station data read from:',csvStationFile

  gatherAndExtract(runTag,dataDir,startTime,endTime,csvStationFile,
                   dataSetName,dataSetType)

if __name__=='__main__' :
  parseCommandLine()
