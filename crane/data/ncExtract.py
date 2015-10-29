"""
New implementation of netcdf output extractor routines.

Tuomas Karna 2013-11-06
"""
import os
import sys
import traceback
import datetime
import time as timeMod
from glob import glob

import numpy as np
from netCDF4 import Dataset as NetCDFFile

from crane.data import timeArray
from crane.data import dataContainer
from crane.data import gridUtils
from crane.data import meshContainer
from crane.data import selfeGridUtils
from crane.files import buildPoints
from crane.files import csvStationFile
from crane import fieldNameList, fieldNameToFilename

#-------------------------------------------------------------------------------
# High-level functions
#     all normal operations should go through these
#-------------------------------------------------------------------------------

def extractTrackForDataContainer(dataDir, trackDC, var, name, bpFile=None,
                                 verbose=False, stacks=None):
  """
  Extracts a track based on x,y,z,time in the given dataContaner.
  """
  # sanity check that trackDC is suitable
  if not trackDC.isTrack() :
    raise Exception( 'given dataContainer does not contain track information' )
  # track metadata
  bracket = trackDC.getMetaData('bracket')
  if not name :
    name = trackDC.getMetaData('location')
  if not var :
    var = trackDC.getMetaData('variable')
  varStr,fileTypeStr = splitVarToFileType(var)
  x = trackDC.x.flatten()
  y = trackDC.y.flatten()
  z = trackDC.z.flatten()

  # Overwrite x,y coords if using alternative coordinate system (open channels)
  if bpFile is not None:
    print 'Overwriting xy from buildpoint'
    bp = buildPoints.BuildPoint()
    bp.readFileFromDisk( bpFile )
    x = bp.getX()
    y = bp.getY()

  # expand scalar coordinates to time-dep track
  nx = max( max( len(x), len(y) ), len(z) )
  if len(x) == 1 : x = np.tile( x, (nx,))
  if len(y) == 1 : y = np.tile( y, (nx,))
  if len(z) == 1 : z = np.tile( z, (nx,))
  zRelToSurf = True if bracket == 'F' else False

  print ' * extracting track', trackDC.getMetaData('location'), varStr
  ee = selfeExtract(dataDir, var=varStr, fileTypeStr=fileTypeStr,
                    verbose=verbose)
  # divide track to smaller chunks to same memory in extraction
  dcs = []
  maxChunkLen = 3000
  chunkIndices = range(0,len(trackDC.time),1000)
  chunkIndices.extend( [len(trackDC.time)] )
  for i in range(len(chunkIndices)-1) :
    tt = trackDC.time.array[chunkIndices[i]:chunkIndices[i+1]]
    xx = x[chunkIndices[i]:chunkIndices[i+1]]
    yy = y[chunkIndices[i]:chunkIndices[i+1]]
    zz = z[chunkIndices[i]:chunkIndices[i+1]]
    try :
      dc = ee.extractTrack(varStr, xx, yy, zz, tt, name, zRelToSurf, stacks=stacks)
      if dc :
        dcs.append(dc)
    except Exception as e :
      print 'Extraction failed:'
      traceback.print_exc(file=sys.stdout)
  # merge chunks
  dc = dcs[0]
  for i in range(1,len(dcs)) :
    dc.mergeTemporal(dcs[i], acceptDuplicates=True)

  return dc

def extractTransectForCoords(x, y, dataDir, varList, startTime, endTime, name,
                             modelCoordSys='spcs',
                             wholeDays=True, stacks=None,
                             verbose=False) :
  """
  Extracts a transect defined by x,y coordinates and returns a dataContainer.
  """
  dcs = []
  for i,var in enumerate(varList) :
    try :
      varStr,fileTypeStr = splitVarToFileType(var)
      ee = selfeExtract(dataDir, var=varStr, fileTypeStr=fileTypeStr,
                        verbose=verbose)
      dc = ee.extractTransect(startTime, endTime, varStr, x,y,name,
                              wholeDays=wholeDays, stacks=stacks)
      dcs.append( dc )
    except Exception as e :
      print 'Extraction failed'
      traceback.print_exc(file=sys.stdout)
  return dcs

def extractTransectForBPFile(bpFile, dataDir, varList, startTime, endTime,
                             name, modelCoordSys=None,
                             wholeDays=True, stacks=None,
                             verbose=False) :
  """
  Extracts a transect defined in the bpFile and returns a dataContainer.
  """
  bpObj = buildPoints.BuildPoint()
  bpObj.readFileFromDisk(bpFile)
  x = bpObj.points[:,1]
  y = bpObj.points[:,2]

  return extractTransectForCoords(x, y, dataDir, varList, startTime, endTime,
                                  name, wholeDays=wholeDays, stacks=stacks,
                                  verbose=verbose)

def extractSlabForLevel(dataDir, varList, startTime, endTime, name,
                        z=None, k=None, zRelToSurf=False,
                        wholeDays=True, stacks=None,
                        verbose=False):
  """
  Extracts a slab at given z coordinate or k level for the given variables.
  """

  mcs = []
  for var in varList :
    varStr,fileTypeStr = splitVarToFileType(var)
    ee = selfeExtract(dataDir, var=varStr, fileTypeStr=fileTypeStr,
                      verbose=verbose)
    mc = ee.extractSlab(startTime, endTime, name, varStr,
                        z, k, zRelToSurf, wholeDays=wholeDays, stacks=stacks)
    mcs.append(mc)
  return mcs

def extractForXYZ( dataDir, var, startTime, endTime, x, y, z=None,
                   stationNames=None, profile=False, zRelToSurf=False,
                   wholeDays=True, stacks=None,
                   verbose=False) :
  """
  Extracts time series for given variable from stations defined by x,y,z.
  If profile=True, will extract profiles instead (z is ignored).
  """
  varStr,fileTypeStr = splitVarToFileType(var)
  ee = selfeExtract(dataDir,var=varStr,fileTypeStr=fileTypeStr,verbose=verbose)
  if not profile and z is None :
    raise Exception('z coordinates must be provided')
  try :
    if profile :
      dcs = ee.extractVerticalProfile(startTime, endTime, varStr, x, y,
                                      stationNames, stacks=stacks)
    else :
      dcs = ee.extractTimeSeries(startTime, endTime, varStr, x, y,
                                 stationNames, z, zRelToSurf=zRelToSurf,
                                 wholeDays=wholeDays, stacks=stacks)
    print ' * extracted'
  except Exception as e :
    print ' * extraction failed'
    traceback.print_exc(file=sys.stdout)
    dcs = []
  for dc in dcs :
    if profile :
      print ' '.join([ dc.getMetaData('location'), dc.getMetaData('variable') ] )
    else :
      print ' '.join([ dc.getMetaData('location'), dc.getMetaData('bracket'),
              dc.getMetaData('msldepth'), dc.getMetaData('variable') ] )
  return dcs

def extractForStations( dataDir, var, stationFile, startTime, endTime,
                        profile=False, wholeDays=True, stacks=None, verbose=False) :
  """
  Extracts time series for given variable from stations defined in stationFile.
  """
  # read station file, allow duplicate stations (with different depth)
  if not os.path.isfile(stationFile):
    raise Exception('File does not exist: '+stationFile)

  if profile:
    csvReader = csvStationFile.csvStationFile()
    csvReader.readFromFile(stationFile)
    tuples = csvReader.getTuples() # all entries (loc,x,y)
    stationNames = [ t[0] for t in tuples ]
    x = np.array([ t[1] for t in tuples ])
    y = np.array([ t[2] for t in tuples ])
    z = None
    print ' *** extracting profiles for stations *** '
    for i,s in enumerate(stationNames) :
      print s,x[i],y[i]
    return extractForXYZ(dataDir, var, startTime, endTime, x, y, z,
                         stationNames, profile=profile, wholeDays=wholeDays,
                         stacks=stacks, verbose=verbose)

  # not profile, depths defined in stationFile
  csvReader = csvStationFile.csvStationFileWithDepth()
  csvReader.readFromFile(stationFile)
  tuples = csvReader.getTuples() # all entries (loc,x,y,z,zType,var)
  stationNames = [ t[0] for t in tuples ]
  x = np.array([ t[1] for t in tuples ])
  y = np.array([ t[2] for t in tuples ])
  z = np.array([ t[3] for t in tuples ])
  zRelToSurf = np.array([ t[4]=='depth' for t in tuples ],dtype=bool)

  dcList = []
  for zIsDepth in [True,False] :
    # filter for z coordinate/depth cases
    ix = np.nonzero(zRelToSurf == zIsDepth)[0]
    if len(ix) == 0 : continue
    x_filt = x[ix]
    y_filt = y[ix]
    z_filt = z[ix]
    stationNames_filt = [stationNames[i] for i in ix ]
    print ' *** extracting for stations *** '
    for i,s in enumerate(stationNames_filt) :
      print s,x_filt[i],y_filt[i],z_filt[i], zRelToSurf[i]
    dcs = extractForXYZ(dataDir, var, startTime, endTime,
                        x_filt, y_filt, z_filt, stationNames_filt,
                        profile, zIsDepth, wholeDays=wholeDays,
                        stacks=stacks, verbose=verbose)
    dcList.extend(dcs)

  return dcList

def extractForOfferings( dataDir, var, offerings, startTime, endTime,
                        profile=False, stationFile=None, fileTypeStr=None,
                        wholeDays=True, verbose=False ) :
  """
  Extracts all stations and depths as in the list of offerings.
  If offerings list is not provided, it will be fetched from the database.

  Args:
    dataDir -- (str) path to model outputs directory
    var     -- (str) variable to extract, e.g. 'elev'. 'temp'
    offerings -- (list of dict) list of offerins from the database
                 each entry has keys 'location','msldepth','bracket',
                 'instrument','variable'
    startTime -- (datetime) first time stamp of extraction
    endTime   -- (datetime)  last time stamp of extraction
    profile -- (bool) if true, extracts vertical profile instead of a value at z
    stationFile   -- (str) a stations.csv file for reading coordinates

  Returns:
    dcList    -- (list of dataContainer) all dataContainers in a list
  """

  # read station file
  staReader = csvStationFile.csvStationFile()
  staReader.readFromFile(stationFile)

  # screen possible duplicates in the offerings (e.g. instrument can be ignored)
  uniqueOffrs = {}
  for o in offerings :
    if var.split('.')[0] in [ 'elev' ] or profile :
      # 2D variable, depth makes no difference
      key = ( o['location'] )
    else :
      key = ( o['location'],o['msldepth'],o['bracket'] )
    uniqueOffrs[ key ] = o
  offerings = uniqueOffrs.values()

  # extract
  dcList = []
  # split offerings to free surface ones ( bracket = 'F' ) and others
  bracketF = [ o for o in offerings if o['bracket'] == 'F' ]
  bracketA = [ o for o in offerings if o['bracket'] != 'F' ]
  for offerings,zRelToSurf in [ (bracketF,True), (bracketA,False)] :
    if not offerings :
      continue
    bracketStr = 'F' if zRelToSurf else 'A'
    print ' *** extracting for offerings: {0:s} bracket *** '.format(bracketStr)

    for o in offerings :
      print tuple( o[k] for k in ['location','msldepth','bracket','variable'] )
    stationNames = [ o['location'] for o in offerings ]
    # omit stations that are missing from sta file
    stationNames = list( set(staReader.getStations()).intersection( set(stationNames) ) )
    offerings = [ o for o in offerings if o['location'] in stationNames ]
    # station of each offering (may contain duplicates)
    stationNames = [ o['location'] for o in offerings ]
    x = np.array( [ staReader.getX(o['location']) for o in offerings ] )
    y = np.array( [ staReader.getY(o['location']) for o in offerings ] )
    zSign = 1 if zRelToSurf else -1 # zRelToSurf => depth below surface
    z = np.array( [ zSign*float(o['msldepth'])/100.0 for o in offerings ] )

    # execute
    dcs = extractForXYZ( dataDir, var, startTime, endTime, x, y, z,
                         stationNames, profile, zRelToSurf,
                         wholeDays=wholeDays,
                         verbose=verbose )
    dcList.extend( dcs )

  return dcList

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

vectorVars = ['hvel','dihv','dahv','wind','wist']

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def splitVarToFileType(var) :
  """Splits 'varname.ext' to 'var','ext'.
  If returns None as extension if cannot split."""
  if len(var.split('.'))==2 :
    return var.split('.')
  else :
    return var,None

def getNCVariableName( var ) :
  """Returns the array name used in netcdf files for given variable."""
  ncVarNames = {'dens':'conc'}
  if var in ncVarNames :
    # custom names
    return ncVarNames[var]
  elif var in fieldNameToFilename and fieldNameToFilename[var][:5] == 'trcr_' :
    # derived from tracer model,trcr_X
    return fieldNameToFilename[var].split('.')[0]
  else :
    # default: same name
    return var

#-------------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------

class selfeNCFile(object) :
  """
  A file that opens selfe netcdf files and provides the metadata (mesh etc.).
  """
  def __init__(self, filename, verbose=False) :
    self.filenamefull = filename
    self.path,self.filename = os.path.split(filename)
    self.fileTypeStr = self.filename.strip('.nc').split('.')[1] # e.g. '70'
    try:
        self.ncfile = NetCDFFile(self.filenamefull,'r')
    except Exception as e:
        print 'Could not read file', self.filename
        raise e
    # expose variables to users
    self.variables = self.ncfile.variables
    self.dimensions = self.ncfile.dimensions
    self.verbose = verbose
    self.headerIsRead = False

  def __delete__(self) :
    if self.verbose :
      print 'closing',self.filenamefull
    if self.ncfile != None :
      self.ncfile.close()
      self.variables = None

  def readHeader(self, meshSearchObj=None) :
    """
    Reads header of the netcdf file and prepares data structures.

    If ncfile is given, will read its header. Otherwise will search for first
    matching netcdf file in the path.
    """

    if self.verbose: print 'Opening file',self.filenamefull

    faceNOffset = self.ncfile.variables['face_nodes'].start_index
    self.faceNodes = self.ncfile.variables['face_nodes'][:].astype(int) - faceNOffset
    self.nodeX = self.ncfile.variables['node_x'][:]
    self.nodeY = self.ncfile.variables['node_y'][:]
    self.bath = self.ncfile.variables['depth'][:]

    self.edgeNodes = self.edgeX = self.edgeY = self.edgeLat = self.edgeLon = None
    if 'edge_nodes' in self.ncfile.variables :
      edgeNOffset = self.ncfile.variables['edge_nodes'].start_index
      self.edgeNodes = self.ncfile.variables['edge_nodes'][:] - edgeNOffset
    if 'edge_x' in self.ncfile.variables :
      self.edgeX = self.ncfile.variables['edge_x'][:]
      self.edgeY = self.ncfile.variables['edge_y'][:]
    if 'edge_lon' in self.ncfile.variables :
      self.edgeLon = self.ncfile.variables['edge_lon'][:]
      self.edgeLat = self.ncfile.variables['edge_lat'][:]

    self.faceX = self.faceY = self.faceLat = self.faceLon = None
    if 'face_x' in self.ncfile.variables :
      self.faceX = self.ncfile.variables['face_x'][:]
      self.faceY = self.ncfile.variables['face_y'][:]
    if 'face_lon' in self.ncfile.variables :
      self.faceLon = self.ncfile.variables['face_lon'][:]
      self.faceLat = self.ncfile.variables['face_lat'][:]

    self.kBottom = np.maximum( self.ncfile.variables['k_bottom'][:]-1, 0 ) # NOTE starts at 1
    #self.elev = self.ncfile.variables['elev'][:]
    self.nNodes = len(self.ncfile.dimensions['node'])
    self.nFaces = len(self.ncfile.dimensions['face'])
    if 'edge' in self.ncfile.dimensions :
      self.nEdges = len(self.ncfile.dimensions['edge'])
    else :
      self.nEdges = None
    self.nVert = len(self.ncfile.dimensions['layers'])
    self.nElemNodes = len(self.ncfile.dimensions['max_face_nodes']) # ==3 always
    self.nTime = len(self.ncfile.dimensions['time'])
    self.exportDt = np.diff(self.ncfile.variables['time'][:])[0]

    self.nDataNodes = self.ncfile.variables['elev'].shape[1]
    if self.nDataNodes == self.nFaces :
      self.discrType = 'face'
    elif self.nDataNodes == self.nEdges :
      self.discrType = 'edge'
    else :
      self.discrType = 'node'
    # whether data is located at half or full levels
    self.vertDiscrType = 'half' if self.fileTypeStr in ['70'] else 'full'
    
    if self.verbose :
      print 'discrType',self.discrType
      print 'nodes',self.nNodes
      print 'elems',self.nFaces
      print 'edges',self.nEdges
      print 'verts',self.nVert
      print 'elem nodes',self.nElemNodes
      print 'time stamps',self.nTime
      print 'export time step',self.exportDt
    timeStr = ' '.join(self.ncfile.variables['time'].base_date.split()[2:4])
    self.simulationStartTime = datetime.datetime.strptime( timeStr, '%Y-%m-%d %H:%M:%S' )

    self.node_lon = self.ncfile.variables['node_lon'][:]
    self.node_lat = self.ncfile.variables['node_lat'][:]

    # vertical coords related
    sigma = self.ncfile.variables['sigma'][:]
    Cs = self.ncfile.variables['Cs'][:]
    if 'nz' in self.ncfile.dimensions :
      nz = len(self.ncfile.dimensions['nz'])
      zz = np.zeros((nz,))
      zz[:nz] = self.ncfile.variables['z'][:]
      h_s = abs(zz[nz-1])
    else : # only S coordinates
      nz = 0
      h_s = 10000.0
      zz = np.zeros((1,))
      zz[0] = -h_s
    h_c = self.ncfile.variables['h_c'][0]
    if self.verbose :
      print 'nz',nz, 'h_s',h_s,'h_c',h_c
    if 'h0' in self.ncfile.variables :
      h0 = self.ncfile.variables['h0'][0]
    else :
      h0 = 0.01 # bug in old combine version: h0 missing

    # init vCoords object
    self.vCoords = selfeGridUtils.verticalCoordinates(self.nVert, nz, h_s, h_c,
                                                      ztot=zz, sigma=sigma, cs=Cs,
                                                      h0=h0)

    # construct mesh search object
    if meshSearchObj != None :
      self.meshSearch2d = meshSearchObj
    else :
      self.meshSearch2d = gridUtils.meshSearch2d( self.nodeX, self.nodeY,
                                                  self.faceNodes, self.edgeNodes )
    self.headerIsRead = True

  def getTime( self ) :
    """Returns time stamps from given netCDF file in epoch format."""
    startTime = ' '.join(self.ncfile.variables['time'].base_date.split()[2:4])
    startTime = datetime.datetime.strptime( startTime, '%Y-%m-%d %H:%M:%S' )
    time = timeArray.simulationToEpochTime( self.ncfile.variables['time'][:], startTime )
    return time

  def getStacks( self, startTime, endTime, wholeDays=False ) :
    """
    Returns a list of file stack numbers that covers the given
    time period [startTime,endTime].

    If wholeDays == True, returns the stacks that contain the whole days
    between startTime and endTime, ignoring hours,minutes etc.

    Simulation start time is read from the netcdf header.
    """
    nSpool = self.nTime # number of exports in each file
    spoolDt = nSpool*self.exportDt

    stackNumber = int(self.filename.split('_')[0])
    # first stack number can differ from 1
    stackOffset = stackNumber - 1
    # first time stamp first file can differ from 0
    timeOffset = float(self.ncfile.variables['time'][0] - self.exportDt)
    offset = (datetime.timedelta(seconds=timeOffset)-
              datetime.timedelta(days=stackOffset))

    if not wholeDays :
      # stacks so that time range covers the requested interval
      startDelta = (startTime - self.simulationStartTime - offset).total_seconds()
      endDelta = (endTime - self.simulationStartTime - offset).total_seconds()
      startDelta -= self.exportDt # include previous day if < 00:15:00
      endDelta -= 1.0 #           # include next     day if > 00:00:01
    else :
      # stacks so that days are included from 00:15 -> 00:00
      sDay = datetime.datetime(startTime.year,startTime.month,startTime.day)
      eDay = datetime.datetime(endTime.year,endTime.month,endTime.day)
      startDelta = (sDay - self.simulationStartTime - offset).total_seconds()
      endDelta = (eDay - self.simulationStartTime - offset).total_seconds()
    startStack = max(int(np.floor(startDelta/spoolDt)) + 1,1)
    endStack = int(np.ceil(endDelta/spoolDt))
    endStack = max( startStack, endStack )

    if endStack <= 0 :
      print endStack, endTime,self.simulationStartTime
      raise Exception('End day less than 1: requested extraction date earlier than simulation start date' )

    return range( startStack, endStack+1 )

  def variableIs3D(self, varStr) :
    """
    Tests whether given variable is defined on 3D mesh.
    """
    return len( self.ncfile.variables[getNCVariableName(varStr)].shape ) == 3

class selfeExtractBase(object) :
  """
  Base class for all SELFE netcdf extract objects
  """
  def __init__(self, path, var=None, verbose=False, fileTypeStr=None,
               meshSearchObj=None) :
    """Intializes reader object."""
    self.verbose = verbose
    self.path = path
    self.varStr = var
    self.externalMeshSeachObj = meshSearchObj
    if fileTypeStr is None :
      try :
        self.fileName = fieldNameToFilename[var]
        self.fileTypeStr = self.fileName.split('.')[1]
      # Try and guess the file name
      except :
          if 'DI' in var or 'DAVG' in var:
            self.fileName = var+'.61'
            self.fileTypeStr = '61'
          else:
            self.fileName = var+'.63'
            self.fileTypeStr = '63'
    else :
      self.fileName = self.varStr+'.'+fileTypeStr
      self.fileTypeStr = fileTypeStr
    self.headerIsRead = False
    self.dataFile = None # file that contains the header for target variable
    self.elevFile = None # file that contains the header for eta and nodal mesh
    self.initialize()

  def initialize(self) :
    """Loads data and elev files and identifies available stacks."""
    self.dataFile = self.getNCFile()
    self.dataFile.readHeader(meshSearchObj=self.externalMeshSeachObj)
    if self.dataFile.discrType == 'node' :
      self.elevFile = self.dataFile
    else :
      # read elev file explicitly to get access no nodal elevations
      self.elevFile = self.getNCFile(fileName='elev.61')
      self.elevFile.readHeader(meshSearchObj=self.dataFile.meshSearch2d)

  def generateFileName(self, iStack=None, fileName=None) :
    """Returns full path to the netcdf file for iStack.
    If iStack==None, returns a pattern with '*' as a wildcard."""
    if fileName is None : fileName = self.fileName
    stackStr = '*' if iStack is None else '{0:d}'.format(iStack) 
    fname = '{stack:s}_{typeStr:s}.nc'.format(
                        typeStr=fileName,stack=stackStr )
    return os.path.join(self.path,fname)
    
  def getNCFile( self, iStack=None, fileName=None ) :
    """Opens netcdf file corresponding to the given stack number.
    If no stack number is given opens first matching file."""
    if fileName is None : fileName = self.fileName
    f = self.generateFileName(iStack, fileName)
    if iStack==None :
      # try to find a file that matches file name pattern
      pattern=f
      files = sorted(glob(pattern))
      if len(files) == 0 :
        raise Exception('no files found in '+pattern)
      stacks = np.array([int(fn.split('/')[-1].split('_')[0]) for fn in files])
      sorted_ix = np.argsort(stacks)
      f = files[sorted_ix[0]]
    if not os.path.isfile(f) :
      raise IOError('File not found: '+f)
    else :
      return selfeNCFile(f,verbose=self.verbose)

  def getVerticalProfile(self, iStack, varStr, x, y, stationNames=None, horzInterp=None) :
    """
    Extracts vertical profiles for the given locations from the given ncfile.
    
    Parameters
    ----------
    iStack : int
           Stack number of the netCDF file to process
    varStr : string
           Variable to extract
    x,y    : array_like (nPoints,)
           Coordinates of the points where to extract
    stationNames : list of strings, optional
           Names of the points for debugging

    Returns
    -------
    time  : array_like (nTime,)
          Time stamps of the extracted data in epoch format
    vals  : array_like (nPoints, nVert, nTime)
          Values of the extracted profiles. For 2d variables nVert=1.
          Missing data is filled with NaNs.
    zcoords : array_like (nPoints, nVert, nTime)
          Z-coordinates for the vertical profiles
    is3d  : bool
          True if extracted variable is defined on 3D mesh
    """
    ncfile = self.getNCFile(iStack)
    if self.dataFile.discrType == 'node' :
      elevfile = ncfile
    else :
      elevfile = self.getNCFile(iStack,fileName='elev.61')

    if horzInterp==None :
      horzInterp = gridUtils.horizontalInterpolator(self.dataFile.meshSearch2d,
                                                    x,y,stationNames)
    nodesToRead = horzInterp.getParentNodes()
    parentElems = horzInterp.parentElems
    nodesToReadData = horzInterp.getParentNodes(self.dataFile.discrType)

    time = ncfile.getTime()
    nTime = len(time)

    if nTime == 0:
      raise Exception('File is corrupted, number of time steps is zero: '+ncfile.filenamefull)
    if nTime != len(elevfile.dimensions['time']):
      print nTime, len(elevfile.dimensions['time'])
      raise Exception('File is corrupted, wrong number of time steps: '+elevfile.filenamefull)

    V = ncfile.variables[getNCVariableName(varStr)]
    is3d = ncfile.variableIs3D(varStr)
    if not is3d : # 2D variable
      if self.verbose: print '2D var', V.shape
      nodalValues = V[:,nodesToReadData][:,None,:] # expand to (nTime,nVert,nUniqNodes)
    else :
      if self.verbose: print '3D var', V.shape
      nodalValues = V[:,:,nodesToReadData]

    nodalValues[nodalValues == -99] = np.nan
    nodalValues = np.swapaxes( nodalValues, 0, 2 ) # reshape to (nUniqNodes,nVert,nTime)
    vals = horzInterp.evaluateArrayFromParentNodes(nodalValues,
                                                   self.dataFile.discrType)
    if is3d :
      # read vertical coords from nodal discretization
      eta = elevfile.variables[getNCVariableName('elev')][:,nodesToRead].flatten()
      dep = np.tile( self.elevFile.bath[nodesToRead], (nTime,) ).flatten()
      Z, kbp2, iwet = self.elevFile.vCoords.computeVerticalCoordinates(eta,dep)
      # reshape from (nVert,nTime*nUniqNodes) to (nUniqNodes,nVert,nTime)
      Z = np.reshape( Z, (self.elevFile.nVert,nTime,-1) )
      Z = np.swapaxes(Z,0,1)
      Z = np.swapaxes(Z,0,2)
      zcoords = horzInterp.evaluateArrayFromParentNodes( Z )
      
      # if data is defined at half levels, convert to full levels
      if self.dataFile.vertDiscrType == 'half' :
        convType = 'naive' #'native'
        vals,zcoords = convertHalfLevelProfileToFullLevel(vals,zcoords,convType)
    else :
      zcoords = np.zeros_like(vals)

    vals = np.ma.masked_invalid(vals)
    zcoords = np.ma.masked_invalid(zcoords)
    return time, vals, zcoords, is3d

  def getVerticalProfileForStacks(self, stacks, varStr, x, y, stationNames=None) :
    """
    Extracts vertical profile for the given netcdf file stacks

    Parameters
    ----------
    stacks : list of int
           Stack numbers to process
    varStr : string
           Variable to extract
    x,y    : array_like (nPoints,)
           Coordinates of the points where to extract
    stationNames : list of strings, optional
           Names of the points for debugging

    Returns
    -------
    time  : array_like (nTime,)
          Time stamps of the extracted data in epoch format
    vals  : array_like (nPoints, nVert, nTime)
          Values of the extracted profiles. For 2d variables nVert=1.
          Missing data is filled with NaNs.
    zcoords : array_like (nPoints, nVert, nTime)
          Z-coordinates for the vertical profiles
    is3d  : bool
          True if extracted variable is defined on 3D mesh
    """

    # create interpolator object for recycling

    horzInterp = gridUtils.horizontalInterpolator(self.dataFile.meshSearch2d,
                                                  x,y,stationNames)

    time = []
    vals = []
    zcoords = []
    for stack in stacks :
      # extract for individual stacks
      try :
        ti,vi,zi,is3d = self.getVerticalProfile(stack,varStr,x,y,
                                                horzInterp=horzInterp)
        time.append(ti)
        vals.append(vi)
        zcoords.append(zi)
      except Exception as e :
        if isinstance(e,IOError) :
          # skip trace for common file not found errors
          print e
        else :
          print 'Extraction failed'
          traceback.print_exc(file=sys.stdout)
    # concatenate time axis
    time = np.concatenate(tuple(time),axis=0) # (nTime,)
    vals = np.concatenate(tuple(vals),axis=2) # (nProfiles,nVert,nTime)
    zcoords = np.concatenate(tuple(zcoords),axis=2) # (nProfiles,nVert,nTime)
    time = np.ma.masked_invalid(time)
    vals = np.ma.masked_invalid(vals)
    zcoords = np.ma.masked_invalid(zcoords)
    return time,vals,zcoords,is3d

  def getTimeSeriesFromProfiles(self, vals, zcoords, z=None, k=None,
                                zRelToSurf=None) :
    """
    Interpolates vertical profiles in vertical at given depth.

    Parameters
    ----------
    vals  : array_like (nPoints, nVert, nTime)
          Values of the extracted profiles. For 2d variables nVert=1.
          Missing data is filled with NaNs.
    zcoords : array_like (nPoints, nVert, nTime)
          Z-coordinates for the vertical profiles
    z : float, array_like (nProfiles,), optional
      z coordinate where each vertical profile is evaluated. z coordinates
      increase upwards.
    k : int, array_like (nProfiles,), optional
      Instead of interpolating, take k-th nodal value from bottom.
      k=1 stands for bottom, k=-1 stands for surface
    zRelToSurf : bool, array_like (nProfiles,), optional
      If True z coordinate is taken depth below free surface instead of
      static z coordinate

    Returns
    -------
    vals : array_like (nProfiles,nTime,)
        Interpolated values
    z_actual : array_like (nProfiles,nTime,)
        The z coordinate at which the interpolation actually took place
    """
    vertInterp = gridUtils.verticalInterpolator(z,k,zRelToSurf)
    v, z_actual = vertInterp.evaluateArray(zcoords,vals)
    v = np.ma.masked_invalid(v)
    return v,z_actual
  
  def getTimeSeries(self, iStack, varStr, x, y, stationNames=None, z=None, k=None, zRelToSurf=None) :
    """Extracts time series from the iStack netcdf file"""
    time, vprof, zcoords, is3d = self.getVerticalProfile(iStack, varStr, x, y, stationNames)
    if not is3d :
      vals = vprof[:,0,:]
      vals = np.ma.masked_invalid(vals)
      z_actual = np.zeros_like(vals)
    else :
      vals, z_actual = self.getTimeSeriesFromProfiles( vprof, zcoords, z, k, zRelToSurf)
    return time,vals,z_actual
  
  def getTimeSeriesForStacks(self, stacks, varStr, x, y, stationNames=None, z=None, k=None, zRelToSurf=None) :
    """
    Extracts time series for the given stacks.

    Parameters
    ----------
    stacks : list of int
           Stack numbers to process
    varStr : string
           Variable to extract
    x,y    : array_like (nPoints,)
           Coordinates of the points where to extract
    stationNames : list of strings, optional
           Names of the points for debugging
    z : float, array_like (nProfiles,), optional
      z coordinate where each vertical profile is evaluated. z coordinates
      increase upwards.
    k : int, array_like (nProfiles,), optional
      Instead of interpolating, take k-th nodal value from bottom.
      k=1 stands for bottom, k=-1 stands for surface
    zRelToSurf : bool, array_like (nProfiles,), optional
      If True z coordinate is taken depth below free surface instead of
      static z coordinate

    Returns
    -------
    time  : array_like (nTime,)
          Time stamps of the extracted data in epoch format
    vals : array_like (nProfiles,nTime,)
        Interpolated values
    z_actual : array_like (nProfiles,nTime,)
        The z coordinate at which the interpolation actually took place
    """
    time, vprof, zcoords, is3d = self.getVerticalProfileForStacks(stacks, varStr, x, y, stationNames)
    if not is3d :
      vals = vprof[:,0,:]
      vals = np.ma.masked_invalid(vals)
      z_actual = np.zeros_like(vals)
    else :
      vals, z_actual = self.getTimeSeriesFromProfiles( vprof, zcoords, z, k, zRelToSurf)
    return time,vals,z_actual

  def getSlab(self, iStack, var, z=None, k=None, zRelToSurf=None) :
    """
    Extracts a horizontal slice from the given ncfile.
    
    Parameters
    ----------
    iStack : int
           Stack number of the netCDF file to process
    var : string
           Variable to extract
    z      : float, optional
           z coordinate where each vertical profile is evaluated. z coordinates
           increase upwards.
    k      : int, optional
           Instead of interpolating, take k-th nodal value from bottom.
           k=1 stands for bottom, k=-1 stands for surface
    zRelToSurf : bool, optional
           If True z coordinate is taken depth below free surface instead of
           static z coordinate

    Returns
    -------
    time  : array_like (nTime,)
          Time stamps of the extracted data in epoch format
    vals  : array_like (nPoints, nTime)
          Values of the extracted horizontal slice.
    zSlab : array_like (nPoints, nTime)
          z coordinates of the extracted horizontal slice.
    """
    ncfile = self.getNCFile(iStack)
    etafile = self.getNCFile(iStack,fileName='elev.61')

    time = ncfile.getTime()
    nTime = len(time)

    varStr = getNCVariableName(var)
    ncVar = ncfile.variables[varStr]
    elev = etafile.variables['elev'][:]
    is3d = ncfile.variableIs3D(varStr)

    nTime = self.elevFile.nTime
    nVert = self.elevFile.nVert
    nNodes = self.elevFile.nNodes
    nElems = self.elevFile.nFaces

    if nTime == 0:
      raise Exception('File is corrupted, number of time steps is zero: '+ncfile.filenamefull)
    if nTime != len(etafile.dimensions['time']):
      print nTime, len(etafile.dimensions['time'])
      raise Exception('File is corrupted, wrong number of time steps: '+etafile.filenamefull)

    if not is3d :
      vals = ncVar[:].T
      zSlab = np.zeros_like(vals)
      dryNodes = np.zeros((nNodes,nTime),dtype=np.bool)
      dryElems = np.zeros((nElems,nTime),dtype=np.bool)
      for iT in range(nTime) :
        eta = elev[iT,:]
        dep = self.elevFile.bath
        dE,dN = self.elevFile.vCoords.computeDryElemMask(eta,dep,self.elevFile.faceNodes)
        dryNodes[:,iT]=dN
        dryElems[:,iT]=dE
    else :
      nDataNodes = self.dataFile.nDataNodes # dataFile may have different nodes
      vals = np.ma.ones((nDataNodes,nTime))*np.nan
      zSlab = np.ma.ones((nDataNodes,nTime))*np.nan
      if k is None and z is None :
        raise Exception('either k or z must be specified for 3D fields')
      # compute z coords for all time steps
      zDim = nVert-1 if self.dataFile.vertDiscrType == 'half' else nVert
      Z = np.zeros((nTime,zDim,nDataNodes),dtype=np.float32)
      dryNodes = np.zeros((nNodes,nTime),dtype=np.bool)
      dryElems = np.zeros((nElems,nTime),dtype=np.bool)
      for iT in range(nTime) : # NOTE computing vcoords is slow
        eta = elev[iT,:]
        dep = self.elevFile.bath
        Zi = self.elevFile.vCoords.computeVerticalCoordinates(eta,dep)[0].T
        if self.dataFile.discrType == 'face' :
          Zi = self.dataFile.meshSearch2d.convertNodalValuesToCenters(Zi)
        elif self.dataFile.discrType == 'edge' :
          Zi = self.dataFile.meshSearch2d.convertNodalValuesToEdges(Zi)
        if self.dataFile.vertDiscrType == 'half' :
          Zi = convertFullLevelProfileToHalfLevel(Zi.T).T
        dE,dN = self.elevFile.vCoords.computeDryElemMask(eta,dep,self.elevFile.faceNodes)
        Z[iT,:,:]=Zi.T
        dryNodes[:,iT]=dN
        dryElems[:,iT]=dE
      if self.dataFile.discrType == 'node' :
        dryDataNodes = dryNodes
      elif self.dataFile.discrType == 'face' :
        dryDataNodes = dryElems
      elif self.dataFile.discrType == 'edge' :
        dryEdges = self.dataFile.meshSearch2d.convertNodalValuesToEdges(dryNodes,np.min)
        dryDataNodes = dryEdges
      alongSLevel = z is None
      # do not read bottom value if data at half levels
      kOffset = 1 if self.dataFile.vertDiscrType == 'half' else 0
      if alongSLevel :
        t0 = timeMod.clock()
        if k < 0 :
          # k=-1 maps to nvrt
          kk = self.elevFile.nVert + k
          vals[:,:] = ncVar[:,kk,:].T # (nTime,nVert,nNodes)
          zSlab[:,:] = Z[:,kk-kOffset,:].T
        else :
          # k=+1 maps to bottom
          kArray = self.dataFile.kBottom + k-1
          uniqueK = np.unique( kArray )
          for kk in uniqueK :
            ix = np.nonzero(kArray == kk)[0]
            v = ncVar[:,kk+kOffset,:].T
            vals[ix,:] = v[ix,:]
            zi = Z[:,kk,:].T
            zSlab[ix,:] = zi[ix,:]
        print timeMod.clock() - t0,'s'
      else :
        # z specified, interpolate in vertical
        t0 = timeMod.clock()
        for iT in range(nTime) :
          Zi = Z[iT,:,:]
          zArr = z*np.ones((nDataNodes,))
          if zRelToSurf :
            # z treated as depth below surface
            zArr = Zi[-1,:] - z
          # find nodes where depth exists
          Zmin = np.nanmin( Zi, axis=0 )
          Zmax = np.nanmax( Zi, axis=0 )
          goodIx = np.logical_and( zArr >= Zmin, zArr <= Zmax )
          goodIx = np.nonzero(np.logical_and( goodIx, ~dryDataNodes[:,iT] ))[0]
          # Vi (nTime,nVert,nNodes)
          Vi = ncVar[iT,:,:]

          # find levels above and below z
          Zdiff = np.abs( Zi[:,goodIx] - zArr[goodIx] )
          absDiffIx = np.argsort( Zdiff, axis=0 )
          kup = absDiffIx[0,:] # indices for up (or down)
          kdw = absDiffIx[1,:] # indices for down (or up)

          Zdw = Zi[kdw,goodIx]
          Zup = Zi[kup,goodIx]
          Vdw = Vi[kdw+kOffset,goodIx]
          Vup = Vi[kup+kOffset,goodIx]
          # find coefficients
          #z = Zdw + x*(Zup-Zdw)
          coef = (zArr[goodIx]-Zdw)/(Zup-Zdw)
          vals[goodIx,iT] = (Vdw + coef*(Vup-Vdw)).T
          zSlab[goodIx,iT] = (Zdw + coef*(Zup-Zdw)).T

        print timeMod.clock() - t0,'s'

    # mask dry areas
    if self.dataFile.discrType == 'node' :
      vals[ dryNodes ] = np.nan
      vals = np.ma.masked_invalid( vals )
    elif self.dataFile.discrType == 'face' :
      vals[ dryElems ] = np.nan
      vals = np.ma.masked_invalid( vals )
      zSlab = self.dataFile.meshSearch2d.averageElemValuesToNodes(zSlab)
    elif self.dataFile.discrType == 'edge' :
      # convert nodal values to dg nodes
      vals = self.dataFile.meshSearch2d.convertEdgeValuesToDGNodes(vals)
      zSlab = self.dataFile.meshSearch2d.convertEdgeValuesToDGNodes(zSlab)
      mask = np.repeat(dryElems,3,axis=0)
      vals[ mask ] = np.nan
      vals = np.ma.masked_invalid( vals )
    return time, vals, zSlab

  def getSlabForStacks(self, stacks, varStr, z=None, k=None, zRelToSurf=None) :
    """
    Returns slab for the given stacks
    """
    time = []
    vals = []
    zcoords = []
    for stack in stacks :
      # extract for individual stacks
      try :
        ti,vi,zi = self.getSlab(stack,varStr,z, k, zRelToSurf)
        time.append(ti)
        vals.append(vi)
        zcoords.append(zi)
      except Exception as e :
        if isinstance(e,IOError) :
          # skip trace for common file not found errors
          print e
        else :
          print 'Extraction failed'
          traceback.print_exc(file=sys.stdout)
    if len(time) == 0 :
      raise Exception('Extraction Failed: no time steps were retrieved')
    # concatenate time axis
    time = np.concatenate(tuple(time),axis=0) # (nTime,)
    vals = np.concatenate(tuple(vals),axis=1) # (nProfiles,nTime)
    zcoords = np.concatenate(tuple(zcoords),axis=1) # (nProfiles,nTime)
    time = np.ma.masked_invalid(time)
    vals = np.ma.masked_invalid(vals)
    zcoords = np.ma.masked_invalid(zcoords)
    return time,vals,zcoords

  def getXYZT(self, varStr, X0, Y0, Z0, T0, zRelToSurf=False, stacks=None):
    """
    Extracts a track defined by the X,Y,Z,T arrays for the given variable.
    
    Parameters
    ----------
    varStr   : string
           Variable to extract
    X,Y,Z,T  : array_like (nPoints,)
           (x,y,z,t) coordinates of the track. Time is assumed to be in epoch format.
    zRelToSurf : bool, optional
             If True z coordinate is taken depth below free surface instead of
           static z coordinate

    Returns
    -------
    X,Y,Z,T  : array_like (nPoints,)
          Coordinates of the extracted track with missing values removed
    actualZ  : array_like (nPoints,)
          Z coordinates where the data was extracted. This is always in z space
          even if input Z would be depth relative to surface.
    data     : array_like (nPoints,)
          Values of the extracted track slice.
    """

    # figure out correct stacks
    if stacks is None:
        stacks = self.dataFile.getStacks(timeArray.epochToDatetime(T0[0]),
                                         timeArray.epochToDatetime(T0[-1]),
                                         wholeDays=False)
    # get profiles for all x,y locations (nPoints,nVert,nTime)
    time,vals,zcor,is3d = self.getVerticalProfileForStacks(stacks,varStr,X0,Y0)
    goodStaIx = np.any( np.any( np.isfinite(vals), axis=1 ), axis=1 )
    # exclude points outside the grid
    X = X0[goodStaIx]
    Y = Y0[goodStaIx]
    Z = Z0[goodStaIx]
    T = T0[goodStaIx]
    vals = vals[goodStaIx,:,:]
    zcor = zcor[goodStaIx,:,:]

    if time[0] > T[0] or time[-1] < T[-1] :
      print 'Extracted time range does not cover query range: cropping'
      if time[0] > T[0] :
        print 'start',timeArray.epochToDatetime(time[0]), '>', timeArray.epochToDatetime(T[0])
      if time[-1] < T[-1] :
        print 'end',timeArray.epochToDatetime(time[-1]), '<', timeArray.epochToDatetime(T[-1])
      goodTimeIx = np.logical_and( T >= time[0], T <= time[-1] )
      print goodTimeIx
      if not np.any(goodTimeIx):
        print 'extracted time range:', timeArray.epochToDatetime(time[0]), '->', timeArray.epochToDatetime(time[0])
        print 'query time range:', timeArray.epochToDatetime(T[0]), '->', timeArray.epochToDatetime(T[0])
        raise Exception('Time arrays do not overlap - cannot extract')
      X = X[goodTimeIx]
      Y = Y[goodTimeIx]
      Z = Z[goodTimeIx]
      T = T[goodTimeIx]

    # array for exported data
    data = np.zeros_like(X)

    actualZ = np.ones_like(Z)*np.nan
    nVert = vals.shape[1]
    nP = len(X)
    # interpolate in vertical
    if nVert==1 :
      # 2D variable, do only temporal interpolation
      for i in range(nP) :
        VV = vals[i,0,:]
        vtmp = interp1d( time,VV )( T[i] )
        if not np.isfinite(vtmp):
          # bad data
          continue
        actualZ[i] = 0.0
        data[i] = vtmp
    else :
      for i in range(nP) :
        # interpolate in time first (z depends on time)
        # z coordinates and values at X[i],Y[i],T[i]
        ZZ = zcor[i,:,:]
        VV = vals[i,:,:]
        goodIx = ~np.all(ZZ.mask,axis=1) # remove mask
        goodIx = np.logical_and( goodIx, ~np.all(np.isnan(VV),axis=1) )
        if not goodIx.any() :
          # all bad data
          continue
        ztmp = interp1d( time,ZZ[goodIx,:] )( T[i] )
        vtmp = interp1d( time,VV[goodIx,:] )( T[i] )
        # interpolate in vertical
        z = Z[i]
        if zRelToSurf :
          # staZ treated as depth below surface
          z = ztmp[-1] - Z[i]
        if z < ztmp[0] or z > ztmp[-1] :
          # z out of bounds
          if z < ztmp[0] :
            z = ztmp[0]
          if z > ztmp[-1] :
            z = ztmp[-1]
          print 'z coordinate changed: ',X[i],Y[i],Z[i],z
        actualZ[i] = z
        data[i] = interp1d( ztmp, vtmp ) ( z )
        if np.isnan(data[i]) :
          print zcor[:,:,i]
          print data[i]
          print z
          print ztmp
          print vtmp
          raise Exception('vertical interpolation failed: nan in result array')

    return X,Y,Z,T,actualZ,data

#-------------------------------------------------------------------------------
# Specialized extraction classes
#-------------------------------------------------------------------------------

class selfeExtract(selfeExtractBase) :
  """
  This class contains only high-level extraction routines and returns the
  data in dataContainer/meshContainer with metadata.
  """
  def __init__(self,path, var=None, verbose=False, fileTypeStr=None):
    selfeExtractBase.__init__(self,path,var,verbose=verbose,
                              fileTypeStr=fileTypeStr)

  def extractTimeSeries(self, startTime, endTime, varStr, staX, staY, stationNames,
                        staZ=None, k=None, zRelToSurf=False, wholeDays=True,
                        stacks=None):
    """
    Extracts time series for the given time range.
    """

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dcs = self.extractTimeSeries(startTime, endTime, var, staX, staY,
                                     stationNames, staZ=staZ, k=k,
                                     zRelToSurf=zRelToSurf,
                                     wholeDays=wholeDays, stacks=stacks)
        compDCs.append( dcs )
      # merge components
      for i in range(len(compDCs[0])) :
        compDCs[0][i].mergeFields( compDCs[1][i] )
        compDCs[0][i].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
        compDCs[0][i].setMetaData( 'variable',varStr )
      return compDCs[0]

    if stacks is None:
        stacks = self.dataFile.getStacks(startTime, endTime, wholeDays=wholeDays)
    zRelToSurfArr = np.ones_like(staX)*zRelToSurf
    time, vals, actualZ = self.getTimeSeriesForStacks(stacks, varStr, staX, staY,
                                                      stationNames, staZ, k,
                                                      zRelToSurfArr)
    # build dataContainer for each station
    dcs = []
    for iSta in range(len(staX)) :
      data = vals[iSta,:]
      # remove nans
      goodIx = np.logical_not( data.mask )
      if not goodIx.any() :
        # all bad data
        print 'all bad data:',stationNames[iSta]
        continue
      # NOTE data,time must be ndarray not masked array
      data = np.reshape( np.array(data[goodIx]), (1,1,-1) )
      t = np.array(time[goodIx])

      ta = timeArray.timeArray( t, 'epoch' )
      meta = {}
      meta['location'] = stationNames[iSta]
      meta['instrument'] = 'model'
      meta['variable'] = varStr
      alongSLevel = k != None
      if alongSLevel :
        meta['slevel'] = kLevel
      else :
        zSign = 1 if zRelToSurf else -1 # zRelToSurf => depth below surface
        zTarget = staZ[iSta]
        msldepth = str(int(round(zSign*zTarget*100)))
        meta['bracket'] = 'F' if zRelToSurf else 'A'
        meta['msldepth'] = msldepth
      meta['dataType'] = 'timeseries'
      z = np.mean(actualZ[iSta,:])
      x = staX[iSta]
      y = staY[iSta]
      dc = dataContainer.dataContainer('', ta, x,y,z, data, fieldNameList.get(varStr,[varStr]),
                          coordSys='spcs',metaData=meta)
      dcs.append(dc)
    return dcs

  def extractVerticalProfile(self,startTime, endTime, varStr, staX, staY,
                             stationNames=None, wholeDays=True, stacks=None):
    """
    Extracts vertical profiles for the given time period.
    """

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dcs = self.extractVerticalProfile(startTime, endTime, var, staX, staY,
                                          stationNames, wholeDays=wholeDays,
                                          stacks=stacks)
        compDCs.append( dcs )
      # merge components
      for i in range(len(compDCs[0])) :
        compDCs[0][i].mergeFields( compDCs[1][i] )
        compDCs[0][i].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
        compDCs[0][i].setMetaData( 'variable',varStr )
      return compDCs[0]

    if stacks is None:
        stacks = self.dataFile.getStacks(startTime,endTime,wholeDays=wholeDays)
    # time:(nTime,) vals zcoords:(nPoints, nVert, nTime)
    time,vals,zcoords,is3d = self.getVerticalProfileForStacks(stacks,varStr,
                                                              staX ,staY,
                                                              stationNames)
    
    if not is3d :
      raise Exception('This is a 2D variable, cannot extract vertical profile: {0:s}'.format(varStr))
    
    goodStaIx = np.any( np.any( np.isfinite(vals), axis=1 ), axis=1 )


    # build dataContainer for each station
    dcs = []
    for iSta in range(len(staX)) :
      staName = '' if stationNames is None else stationNames[iSta]
      staStr = '{x:f} {y:f} {name:s}'.format(x=staX[iSta],y=staY[iSta],name=staName)
      # remove time steps with all bad values
      goodIxTime = np.logical_and( ~np.all( vals[iSta,:,:].mask, axis=0 ),
                                   ~np.all( zcoords[iSta,:,:].mask, axis=0 ) )
      goodIxTime = np.nonzero( goodIxTime )[0]
      v = vals[iSta,:,:][:,goodIxTime]
      z = zcoords[iSta,:,:][:,goodIxTime]
      t = time[goodIxTime]
      # remove masked (below bottom) part (nVert,nTime) -> (nGoodVert,nTime)
      goodIxVert = np.logical_and( ~np.any( v.mask, axis=1 ),
                                   ~np.any( z.mask, axis=1 ) )

      if not goodIxVert.any() or not goodIxTime.any() :
        # all bad data
        print not goodIxVert.any(), not goodIxTime.any()
        print 'all bad data: ',staStr
        continue
      v = v[goodIxVert,:]
      z = z[goodIxVert,:]
      if v.mask.any() or z.mask.any() or t.mask.any() :
        print v.mask.any(), z.mask.any(), t.mask.any()
        print v.mask
        print z.mask
        print t.mask
        raise Exception('bad values remain: '+staStr)
      # to (nGoodVert,1,nTime)
      data = v[:,None,:]
      ta = timeArray.timeArray( np.array(t), 'epoch' )
      # to (nGoodVert,nTime)
      nZ = z.shape[0]
      x = staX[iSta]*np.ones((nZ,))
      y = staY[iSta]*np.ones((nZ,))
      meta = {}
      meta['location'] = stationNames[iSta]
      meta['instrument'] = 'model'
      meta['bracket'] = 'A'
      meta['variable'] = varStr
      meta['dataType'] = 'profile'
      dc = dataContainer.dataContainer('', ta, x,y,z, data, fieldNameList.get(varStr,[varStr]),
                          coordSys='spcs',metaData=meta)
      dcs.append(dc)
    return dcs
  
  def extractTransect(self,startTime, endTime,varStr,staX,staY,transName,
                      wholeDays=True, stacks=None):
    """
    Extracts a transect for the given (x,y) points and time range.
    """

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dc = self.extractTransect(startTime, endTime, var, staX, staY, transName,
                                  wholeDays=wholeDays, stacks=stacks)
        compDCs.append( dc )
      # merge components
      compDCs[0].mergeFields( compDCs[1] )
      compDCs[0].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
      compDCs[0].setMetaData('variable',varStr)
      return compDCs[0]

    if stacks is None:
        stacks = self.dataFile.getStacks(startTime,endTime,wholeDays=wholeDays)
    staX = np.array(staX)
    staY = np.array(staY)
    # time:(nTime,) vals zcoords:(nPoints, nVert, nTime)
    time,vals,zcoords,is3d = self.getVerticalProfileForStacks(stacks, varStr,
                                                              staX, staY)

    # reorganize transect in trasect format
    data = []
    X = []
    Y = []
    Z = []
    for iSta in range(len(staX)) :
      staStr = '{ix:d} {x:f} {y:f}'.format(ix=iSta, x=staX[iSta], y=staY[iSta])

      z = zcoords[iSta,:,:]
      v = vals[iSta,:,:]
      # remove masked (below bottom) part (nVert,nTime) -> (nGoodVert,nTime)
      goodIxVert = np.logical_and( ~np.all( v.mask, axis=1 ),
                                   ~np.all( z.mask, axis=1 ) )
      if not goodIxVert.any() :
        # all bad data
        print 'all bad data: ',staStr
        continue
      z = z[goodIxVert,:].filled(np.nan)
      v = v[goodIxVert,:].filled(np.nan)
      
      x = staX[iSta]*np.ones((z.shape[0],))
      y = staY[iSta]*np.ones((z.shape[0],))
      data.append( v ) # (nVert,nTime)
      X.append( x ) # (nVert,)
      Y.append( y ) # (nVert,)
      Z.append( z ) # (nVert,nTime)
    # concatenate
    X = np.concatenate( tuple(X), axis=0 )
    Y = np.concatenate( tuple(Y), axis=0 )
    Z = np.concatenate( tuple(Z), axis=0 )
    data = np.concatenate( tuple(data), axis=0 )
    # reshape for dataContainer
    # from (nVert*nSta,nTime) to (nVert*nSta,1,nTime)
    data = data[:,None,:]
    
    # build dataContainer
    ta = timeArray.timeArray( time, 'epoch' )
    meta = {}
    meta['location'] = transName
    meta['instrument'] = 'model'
    meta['bracket'] = 'A'
    meta['variable'] = varStr
    meta['dataType'] = 'transect'
    dc = dataContainer.dataContainer('', ta, X,Y,Z, data, fieldNameList.get(varStr,[varStr]),
                       coordSys='spcs', metaData=meta, acceptNaNs=True)
    return dc

  def extractTrack(self, varStr, X, Y, Z, T, trackName,
                   zRelToSurf=False, stacks=None):
    """Extracts track defined by X,Y,Z,T arrays for given stacks
       and variables."""

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dc = self.extractTrack(var,X,Y,Z,T,trackName,zRelToSurf)
        compDCs.append( dc )
      # merge components
      compDCs[0].mergeFields( compDCs[1] )
      compDCs[0].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
      compDCs[0].setMetaData('variable',varStr)
      return compDCs[0]

    XX,YY,ZZ,TT,actualZ,data = self.getXYZT(varStr, X, Y, Z, T,zRelToSurf, stacks)

    # create dataContainer
    ta = timeArray.timeArray(TT, 'epoch', acceptDuplicates=True)
    data = data[None,None,:]
    if not zRelToSurf :
      # export the actual z coordinate where data was extracted
      z = actualZ[None,:]
    else :
      z = ZZ[None,:]

    meta = {}
    meta['dataType'] = 'track'
    meta['location'] = trackName
    meta['instrument'] = 'model'
    meta['bracket'] = 'F' if zRelToSurf else 'A'
    meta['variable'] = varStr
    dc = dataContainer.dataContainer('', ta, XX[None,:],YY[None,:],z, data,
                       fieldNameList.get(varStr,[varStr]),
                       coordSys='spcs',metaData=meta,acceptNaNs=True)

    return dc

  def extractSlab(self,startTime, endTime, name, varStr, z=None, k=None,
                  zRelToSurf=None, wholeDays=True, stacks=None):
    """Extracts a horiontal slice for the given time range."""

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dc = self.extractSlab(startTime, endTime, name, var, z, k,
                              zRelToSurf, wholeDays=wholeDays, stacks=stacks)
        compDCs.append( dc )
      # merge components
      compDCs[0].mergeFields( compDCs[1] )
      compDCs[0].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
      compDCs[0].setMetaData('variable',varStr)
      return compDCs[0]

    if stacks is None:
        stacks = self.dataFile.getStacks(startTime,endTime,wholeDays=wholeDays)
    time,vals,zcoords = self.getSlabForStacks(stacks, varStr, z, k, zRelToSurf)
    ta = timeArray.timeArray( time, 'epoch' )
    vals = vals.filled(np.nan)
    data = vals[:,None,:]
    zcoords = zcoords.filled(np.nan)
    connectivity = self.dataFile.faceNodes
    x = self.dataFile.nodeX
    y = self.dataFile.nodeY
    if data.shape[0] == 3*self.dataFile.nFaces :
      # 3 nodes per element -> construct discontinuous mesh
      x = x[connectivity].ravel()
      y = y[connectivity].ravel()
      connectivity = np.reshape(np.arange(len(x)),(-1,3))
    msldepth = ''
    if k != None :
      msldepth = 'slev'+str(k)
      zArray = np.mean(zcoords,axis=1) # FIXME include time in z coords ?
    else :
      zSign = 1 if zRelToSurf else -1 # zRelToSurf => depth below surface
      msldepth = str(int(round(zSign*z*100)))
      zArray = z*np.ones_like(x)

    # need to reconstruct z coords?? hacky
    # make a special vcoords evaluator that reads elev.61 if necessary?
    # implement selfeNCFile class that reads header, comps vcoords, getsNodal data
    # selfeExtractBase has two files, dataFile and elevFile, latter always 63 file

    # make meshContainer
    meta = {}
    meta['dataType'] = 'slab'
    meta['location'] = name
    meta['instrument'] = 'model'
    meta['variable'] = varStr
    if k != None :
      meta['slevel'] = k
    else :
      meta['bracket'] = 'F' if zRelToSurf else 'A'
      meta['msldepth'] = msldepth
    mc = meshContainer.meshContainer('', ta, x,y,zArray, data, connectivity,
                       fieldNameList.get(varStr,[varStr]),coordSys='spcs',
                       metaData=meta,)
    return mc
