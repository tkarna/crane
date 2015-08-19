
import os
import sys
import traceback
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from scipy.spatial import KDTree, cKDTree
from scipy.interpolate import interp1d
from glob import glob
import time as timeMod

from crane.data import timeArray
from crane.data import dataContainer
from data.meshContainer import meshContainer
from data.extractStation import fieldNameList, fieldNameToFilename
from data.selfeGridUtils import *
from files.csvStationFile import csvStationFile, csvStationFileWithDepth
from files.buildPoints import BuildPoint

vectorVars = ['hvel','dihv','dahv','wind','wist']

ncVarNames = {'dens':'conc'}

def getNCVariableName( var ) :
  """Returns the array name used in netcdf files for given variable."""
  if var in ncVarNames :
    # custom names
    return ncVarNames[var]
  elif var in fieldNameToFilename and fieldNameToFilename[var][:5] == 'trcr_' :
    # derived from tracer model,trcr_X
    return fieldNameToFilename[var].split('.')[0]
  else :
    # default: same name
    return var

def extractForXYZ( dataDir, var, startTime, endTime, x, y, z=None,
                   stationNames=None, profile=False, zRelativeToSurf=False ) :
  """Extracts time series for given variable from stations defined by x,y,z.
  If profile=True, will extract profiles instead (z is ignored)."""

  if profile :
    ee = ncExtractProfile(dataDir,var=var)
  else :
    if z==None :
      raise Exception('z coordinates must be provided')
    ee = ncExtractStation(dataDir,var=var)
  try :
    if profile :
      dcs = ee.extractDates(startTime,endTime,var,x,y,stationNames)
    else :
      dcs = ee.extractDates(startTime,endTime, var, x, y, z, stationNames,
                            zRelativeToSurf=zRelativeToSurf)
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
                        profile=False, zRelativeToSurf=False ) :
  """Extracts time series for given variable from stations defined in stationFile."""
    # read station file, allow duplicate stations (with different depth)
  if not os.path.isfile(stationFile):
    raise Exception('File does not exist: '+stationFile)

  if profile:
    csvReader = csvStationFile()
    csvReader.readFromFile(stationFile)
    tuples = csvReader.getTuples() # all entries (loc,x,y)
    stationNames = [ t[0] for t in tuples ]
    x = np.array([ t[1] for t in tuples ])
    y = np.array([ t[2] for t in tuples ])
    z = None
    print ' *** extracting profiles for stations *** '
    for i,s in enumerate(stationNames) :
      print s,x[i],y[i]
    return extractForXYZ( dataDir, var, startTime, endTime, x, y, z,
                    stationNames, profile, zRelativeToSurf )

  # not profile, depths defined in stationFile
  csvReader = csvStationFileWithDepth()
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
                        profile, zIsDepth )
    dcList.extend(dcs)

  return dcList

def extractForOfferings( dataDir, var, offerings, startTime, endTime,
                        profile=False, stationFile=None ) :
  """
  Extracts all stations and depths as in the list of offerings.
  If offerings list is not provided, it will be fetched from the database.

  Args:
    dataDir -- (str) path to model outputs directory
    var     -- (str) variable to extract, e.g. 'elev'. 'temp'
    offerings -- (list of str) list of offerins from the database
                 each string is station.msldepth.bracket.instrument[.var]
    startTime -- (datetime) first time stamp of extraction
    endTime   -- (datetime)  last time stamp of extraction
    profile -- (bool) if true, extracts vertical profile instead of a value at z
    stationFile   -- (str) a stations.csv file for reading coordinates

  Returns:
    dcList    -- (list of dataContainer) all dataContainers in a list
  """

  # read station file
  staReader = csvStationFile()
  staReader.readFromFile(stationFile)

  # screen possible duplicates in the offerings (e.g. instrument can be ignored)
  uniqueOffrs = {}
  for o in offerings :
    if var in [ 'elev' ] or profile :
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
                         stationNames, profile, zRelToSurf )
    dcList.extend( dcs )

  return dcList


class ncExtractBase(object) :
  """Base class for all netcdf extract objects."""
  def __init__(self, path, fileTypeStr=None, var=None ) :
    """Intializes reader object."""
    self.path = path
    if not fileTypeStr :
      fileTypeStr = fieldNameToFilename[var]+'.nc'
    self.fileTypeStr = fileTypeStr
    self.headerIsRead = False

  def readHeader(self, ncfile=None, verbose=False) :
    """
    Reads header of the netcdf file and prepares data structures.

    If ncfile is given, will read its header. Otherwise will search for first
    matching netcdf file in the path.
    """
    if ncfile == None :
      # try to find a file that matches specs
      pattern = os.path.join(self.path,'*_'+self.fileTypeStr)
      files = sorted(glob(pattern))
      if len(files) == 0 :
        raise Exception('no files found in '+pattern)
      if verbose :
        print 'Reading header from',files[0]
      ncfile = NetCDFFile( files[0] )
    # read
    faceNOffset = ncfile.variables['face_nodes'].start_index
    self.faceNodes = ncfile.variables['face_nodes'][:].astype(int) - faceNOffset
    self.nodeX = ncfile.variables['node_x'][:]
    self.nodeY = ncfile.variables['node_y'][:]
    self.bath = ncfile.variables['depth'][:]

    self.kBottom = np.maximum( ncfile.variables['k_bottom'][:]-1, 0 ) # NOTE starts at 1
    #self.elev = ncfile.variables['elev'][:]
    self.nNodes = len(ncfile.dimensions['node'])
    self.nFaces = len(ncfile.dimensions['face'])
    self.nVert = len(ncfile.dimensions['layers'])
    self.nElemNodes = len(ncfile.dimensions['max_face_nodes']) # ==3 always
    self.nTime = len(ncfile.dimensions['time'])
    if verbose :
      print 'nodes',self.nNodes
      print 'elems',self.nFaces
      print 'verts',self.nVert
      print 'elem nodes',self.nElemNodes
      print 'time stamps',self.nTime
    timeStr = ' '.join(ncfile.variables['time'].base_date.split()[2:4])
    self.simulationStartTime = datetime.datetime.strptime( timeStr, '%Y-%m-%d %H:%M:%S' )

    self.node_lon = ncfile.variables['node_lon'][:]
    self.node_lat = ncfile.variables['node_lat'][:]

    if 'edge_nodes' in ncfile.variables :
      self.edgeNodes = ncfile.variables['edge_nodes'][:]
    if 'edge_x' in ncfile.variables :
      self.edge_x = ncfile.variables['edge_x'][:]
      self.edge_y = ncfile.variables['edge_y'][:]
    if 'edge_lon' in ncfile.variables :
      self.edge_lon = ncfile.variables['edge_lon'][:]
      self.edge_lat = ncfile.variables['edge_lat'][:]

    #from data.coordSys import spcs2lonlat, lonlat2spcs
    #for i in range(20) :
      #lo,la = spcs2lonlat( self.nodeX[i], self.nodeY[i] )
      #print self.node_lon[i],self.node_lat[i], lo, la
      #print '  ',self.node_lon[i]-lo,self.node_lat[i]-la
      #xx,yy = lonlat2spcs( self.node_lon[i],self.node_lat[i] )
      #print self.nodeX[i], self.nodeY[i]
      #print '  ',self.nodeX[i]-xx, self.nodeY[i]-yy
    #exit(0)

    # vertical coords related
    sigma = ncfile.variables['sigma'][:]
    Cs = ncfile.variables['Cs'][:]
    if 'nz' in ncfile.dimensions :
      nz = len(ncfile.dimensions['nz'])
      zz = np.zeros((nz,))
      zz[:nz] = ncfile.variables['z'][:]
      h_s = abs(zz[nz-1])
    else : # only S coordinates
      nz = 0
      h_s = 10000.0
      zz = np.zeros((1,))
      zz[0] = -h_s
    h_c = ncfile.variables['h_c'][0]
    if verbose :
      print nz, h_s,h_c
    if 'h0' in ncfile.variables :
      h0 = ncfile.variables['h0'][0]
    else :
      h0 = 0.01 # bug in old combine version: h0 missing

    # init vCoords object
    self.vCoords = verticalCoordinates(self.nVert, nz, h_s, h_c, ztot=zz, sigma=sigma, cs=Cs, h0=h0 )
    # construct KDTree
    self.tree = constructKDTree( self.nodeX, self.nodeY )
    # construct node to elem table
    self.node2elem = constructNodeToElemTable( self.faceNodes )
    self.headerIsRead = True

  def getTime( self, ncfile ) :
    """Returns time stamps from given netCDF file in epoch format."""
    nTime = len(ncfile.dimensions['time'])
    startTime = ' '.join(ncfile.variables['time'].base_date.split()[2:4])
    startTime = datetime.datetime.strptime( startTime, '%Y-%m-%d %H:%M:%S' )
    time = simulationToEpochTime( ncfile.variables['time'][:], startTime )
    return time

  def getVProf(self, ncfile, varStr, staX, staY, iTri, u, nix) :
    """Extracts vertical profiles for the given locations from the given ncfile.

    The return arrays are masked for missing values (e.g. below bottom part and
    dry periods).

    Args:
      ncfile -- (NetCDFFile) output file where data is extracted from
      varStr     -- (str) variable to extract, e.g. 'elev'. 'temp'
      staX,staY -- (ndarray (nSta,)) x,y coordinates for each station
      iTri -- (ndarray (nSta,)) Indices of parent triangles for each station
      u    -- (ndarray (nSta,3)) Barycentric coordnates for each station
      nix  -- (ndarray (nSta,3)) Triangle vertex indices for each station
    Returns:
      times -- (ndarray (nTime,)) Concatenated time stamps
      values -- (ndarray (nTime,nVert,nSta)) extracted values
      zcoords -- (ndarray (nTime,nVert,nSta)) z coordinates
    """
    if not self.headerIsRead :
      raise Exception('nc file header is not read')
    # brief sanity check
    nNodes = len(ncfile.dimensions['node'])
    nFaces = len(ncfile.dimensions['face'])
    nVert = len(ncfile.dimensions['layers'])
    elev = ncfile.variables['elev'][:]
    if nNodes != self.nNodes :
      raise Exception('number of nodes do not match',nNodes,self.nNodes)
    if nFaces != self.nFaces :
      raise Exception('number of faces do not match',nFaces,self.nFaces)
    if nVert != self.nVert :
      raise Exception('number of vertical layers do not match',nVert,self.nVert)

    # read time information for each file separately
    time = self.getTime( ncfile )
    nTime = len(time)
    nSta = len(staX)

    # read data only for unique nodes (deflate nix)
    uniNix, uniNixInvIx = np.unique( nix.flatten(), return_inverse=True )

    # read variable
    V = ncfile.variables[getNCVariableName(varStr)]
    is3d = len( V.shape ) == 3
    if not is3d : # 2D variable
      nVert = 1
      var = V[:,uniNix]
      if 0 in var.shape :
        raise Exception('Zero dimension in variable array')
      var[var == -99] = np.nan
      var = np.ma.masked_invalid( var )
      # inflate from (nTime, nUniqueNodes) to (nTime, nSta*3)
      var = var[:,uniNixInvIx]
      # reshape var from (nTime, nSta*3) to (nTime, nVert, nSta, 3)
      var = np.reshape( var, (var.shape[0],1,-1,3) )
    else : # 3D variable
      var = V[:,:,uniNix]
      if 0 in var.shape :
        raise Exception('Zero dimension in variable array')
      var[var == -99] = np.nan
      var = np.ma.masked_invalid( var )
      # inflate from (nTime, nVert, nUniqueNodes) to (nTime,nVert, nSta*3)
      var = var[:,:,uniNixInvIx]
      # reshape var from (nTime, nVert, nSta*3) to (nTime, nVert, nSta, 3)
      var = np.reshape( var, (var.shape[0],var.shape[1],-1,3) )

    # interpolate in horizontal
    vals = np.ones( (nTime,nVert,len(staX)) )*np.nan
    for i in range(len(staX)) :
      v = var[:,:,i,:].flatten()
      vali = interpolateTri(u[i,:], np.reshape( v, (-1,3)) )
      # to shape (nTime, nVert)
      vali = np.reshape( vali, (nTime,-1) )
      vals[:,:,i] = vali

    zcors = np.ma.ones( (nTime,nVert,len(staX)) )*np.nan
    if is3d :
      eta = elev[:,nix].flatten()
      dp = np.tile( self.bath[nix.flatten()], (nTime,) )
      # eta,dp in (nTime*nSta*3) shape
      Z, kbp2, iwet = self.vCoords.computeVerticalCoordinates(eta,dp)
      # reshape from (nVert, nTime*nSta*3) to (nTime, nVert, nSta, 3)
      Z = np.reshape( Z, (self.nVert,nTime, nSta,3) )
      Z = np.swapaxes(Z,0,1)
      for i in range(len(staX)) :
        z = Z[:,:,i,:].flatten()
        zi = interpolateTri(u[i,:], np.reshape( z, (-1,3)) )
        # to shape (nTime, nVert)
        zi = np.reshape( zi, (nTime,-1) )
        zcors[:,:,i] = zi
    else :
      zcors = np.zeros( (nTime,1,len(staX)) )

    vals = np.ma.masked_invalid(vals)
    zcors = np.ma.masked_invalid(zcors)

    return time, vals, zcors

  def getNCFile( self, stack ) :
    f = os.path.join( self.path, str(stack)+'_'+self.fileTypeStr )
    if not os.path.isfile(f) :
      print 'file not found',f
    else :
      return NetCDFFile(f,'r')

  def getProfForStacks(self, stacks, varStr, staX, staY, stationNames = None ) :
    """
    Extracts vertical profiles for given station locations and stacks.

    The return arrays are masked for missing values (e.g. below bottom part and
    dry periods).

    Args:
      varStr     -- (str) variable to extract, e.g. 'elev'. 'temp'
      staX,staY -- (ndarray (nSta,)) x,y coordinates for each station
      stationNames -- (list of string) station names (optional)
    Returns:
      times -- (masked array (nTime,)) Concatenated time stamps
      values -- (masked array (nTime,nVert,nOKSta)) extracted values
      zcoords -- (masked array (nTime,nVert,nOKSta)) z coordinates
      goodStaIx -- (masked array (nOKSta,)) indices of good stations (inside the grid)
    """
    staX = np.array(staX)
    staY = np.array(staY)
    iTri = u = nix = staList = None
    nStacks = len(stacks)
    times=values=zcoords=None # allocate later
    for i in range(len(stacks)) :
      stack = stacks[i]
      f = os.path.join( self.path, str(stack)+'_'+self.fileTypeStr )
      if not os.path.isfile(f) :
        print 'file not found, skipping',f
        continue
      ncfile = NetCDFFile(f,'r')
      if not self.headerIsRead :
        self.readHeader(ncfile)
      if iTri == None :
        # find parent triangles
        iTri, u, nix, goodStaIx = findStations( self.tree, self.node2elem,
                                          self.nodeX, self.nodeY, self.faceNodes,
                                          staX, staY, stationNames )
        staX = staX[goodStaIx]
        staY = staY[goodStaIx]
      if times==None :
        nSta = len(staX)
        times = np.ma.ones( (self.nTime*nStacks,) )*np.nan
        values = np.ma.ones( (self.nTime*nStacks,self.nVert,nSta) )*np.nan
        zcoords = np.ma.ones( (self.nTime*nStacks,self.nVert,nSta) )*np.nan
      try :
        ti,vi,zi = self.getVProf( ncfile, varStr, staX, staY, iTri, u, nix )
        times[i*self.nTime:(i+1)*self.nTime] = ti
        values[i*self.nTime:(i+1)*self.nTime,:,:] = vi
        zcoords[i*self.nTime:(i+1)*self.nTime,:,:] = zi
      except Exception as e :
        print 'Extraction failed',f
        traceback.print_exc(file=sys.stdout)
    if times==None :
      raise Exception('no data could be extracted.')
    # remove skipped stacks
    goodIx = np.logical_not( np.isnan( times ) )

    return times[goodIx],values[goodIx,:,:],zcoords[goodIx,:,:],goodStaIx

  def getSlab(self, ncfile, varStr, zCoord=None, kLevel=None,
              zRelativeToSurf=False, iTri=None) :
    """Extracts horizontal slabs from the given ncfile.
    """
    if not self.headerIsRead :
      raise Exception('nc file header is not read')
    # brief sanity check
    nNodes = len(ncfile.dimensions['node'])
    nFaces = len(ncfile.dimensions['face'])
    nVert = len(ncfile.dimensions['layers'])
    elev = ncfile.variables['elev'][:]
    if nNodes != self.nNodes :
      raise Exception('number of nodes do not match',nNodes,self.nNodes)
    if nFaces != self.nFaces :
      raise Exception('number of faces do not match',nFaces,self.nFaces)
    if nVert != self.nVert :
      raise Exception('number of vertical layers do not match',nVert,self.nVert)

    time = self.getTime( ncfile )
    nTime = len(time)

    # quick'n dirty implementation: reads the whole file to memory
    is3d = len( ncfile.variables[varStr].shape ) == 3
    if is3d :
      if kLevel == None and zCoord == None :
        raise Exception('either kLevel or zCoord must be specified for 3D fields')
      alongSLevel = zCoord == None
      if alongSLevel :
        if kLevel > 0 :
          t0 = timeMod.clock()
          # kLevel=+1 maps to bottom
          kArray = self.kBottom + kLevel-1
          var = np.ones((nTime,nNodes))
          ks = np.unique( kArray )
          for k in ks :
            ix = np.nonzero(kArray == k)[0]
            # TODO test which is faster...
            #var[:,ix] = ncfile.variables[varStr][:,k,ix] # (nTime,nVert,nNodes)
            v = ncfile.variables[varStr][:,k,:] # (nTime,nVert,nNodes)
            var[:,ix] = v[:,ix]
          print timeMod.clock() - t0,'s'
        else :
          # kLevel=-1 maps to nvrt
          k = self.nVert + kLevel
          var = ncfile.variables[varStr][:,k,:] # (nTime,nVert,nNodes)
      else :
        # zCoord specified, interpolate in vertical
        # TODO optimize, store kup,kdw from previous step
        var = np.ma.ones((nTime,nNodes))*np.nan
        t0 = timeMod.clock()
        for iT in range(nTime) :
          eta = elev[iT,:]
          dp = self.bath[:]
          # Z (nVert, nNodes)
          Z, kbot, iwet = self.vCoords.computeVerticalCoordinates(eta,dp)
          z = zCoord*np.ones((nNodes,))
          if zRelativeToSurf :
            # staZ treated as depth below surface
            z = Z[-1,:] - zCoord
          # find nodes where depth exists
          Zmin = np.min( Z, axis=0 )
          Zmax = np.max( Z, axis=0 )
          goodIx = np.logical_and( z >= Zmin, z <= Zmax )
          goodIx = np.nonzero(np.logical_and( goodIx, iwet ))[0]
          # V (nTime,nVert,nNodes)
          V = ncfile.variables[varStr][iT,:,:]

          # find levels above and below z
          Zdiff = np.abs( Z[:,goodIx] - z[goodIx] )
          six = np.argsort( Zdiff, axis=0 )
          kup = six[0,:] # indices for up (or down)
          kdw = six[1,:] # indices for down (or up)

          Zdw = Z[kdw,goodIx]
          Zup = Z[kup,goodIx]
          Vdw = V[kdw,goodIx]
          Vup = V[kup,goodIx]
          # find coefficients
          #z = Zdw + x*(Zup-Zdw)
          coef = (z[goodIx]-Zdw)/(Zup-Zdw)
          var[iT,goodIx] = Vdw + coef*(Vup-Vdw)

          ## check
          #iii = 3
          #print z[goodIx[iii]],Zdw[iii] + coef[iii]*(Zup[iii]-Zdw[iii])
          #print var[iT,goodIx[iii]]
          ##print Z[:,goodIx[iii]]
          #print kup[iii],Z[kup[iii],goodIx[iii]],V[kup[iii],goodIx[iii]]
          #print kdw[iii],Z[kdw[iii],goodIx[iii]],V[kdw[iii],goodIx[iii]]
        var = np.ma.masked_invalid( var )
        print timeMod.clock() - t0,'s'
    else :
      var = ncfile.variables[varStr][:]

    # mask dry areas
    dryMask = self.vCoords.computeDryMask(elev,self.bath)
    var = np.ma.MaskedArray( var, mask=dryMask )
    return time, var

  def getSlabForStacks(self, stacks, varStr, zCoord=None, kLevel=None,
                        zRelativeToSurf=False, iTri=None) :
    """Extract """
    times=[]
    values=[]
    for i in range(len(stacks)) :
      stack = stacks[i]
      f = os.path.join( self.path, str(stack)+'_'+self.fileTypeStr )
      if not os.path.isfile(f) :
        print 'file not found, skipping',f
        continue
      ncfile = NetCDFFile(f,'r')
      if not self.headerIsRead :
        self.readHeader(ncfile)
      # (nTime,) (nTime,nNodes)
      ti, vi = self.getSlab( ncfile, varStr, zCoord, kLevel, zRelativeToSurf, iTri)
      times.append(ti)
      values.append(vi)

    times = np.ma.concatenate( tuple(times), axis=0 )
    values = np.ma.concatenate( tuple(values), axis=0 )

    return times, values

  def getStacks( self, startTime, endTime, ncfile=None, exactBounds=False ) :
    """Returns a list of stacks for the given time period bounded by
    startTime,endTime.

    Simulation start time is read from the netcdf header.
    NOTE: startTime,endTime is converted to midnight of the given date, and
    the days in between are considered.
    """
    # deduce correct stacks
    if not self.headerIsRead :
      self.readHeader( ncfile )
    startStack = (startTime - self.simulationStartTime).days + 1
    endStack = (endTime - self.simulationStartTime).days
    endDay = datetime.datetime(endTime.year,endTime.month,endTime.day)
    pastMidnight = endTime > endDay
    if exactBounds and pastMidnight :
      # last time step is during the next day
      endStack += 1
    if exactBounds and (startTime.hour==0 and startTime.minute<15) :
      # first time step is before 00:15 thus in previous stack
      startStack -= 1
    endStack = max( startStack, endStack )
    if startStack < 0 :
      print startStack, self.simulationStartTime
      raise Exception('Negative start day: requested extraction date earlier than simulation start date' )
    return range( startStack, endStack+1 )

  def extract( self, *args, **kwargs ) :
    raise NotImplementedError('This method must be implemented in the derived class' )

  def extractDates( self, startTime,endTime, *args, **kwargs ) :
    """Computes the correct stacks for the given time range and calls the extract routine. args and kwargs are passed to the extract routine."""
    stacks = self.getStacks( startTime,endTime )
    return self.extract( stacks, *args, **kwargs )

class ncExtractStation(ncExtractBase) :
  """A class for extracting station time series from SELFE netCDF4
  output files."""

  def extract(self, stacks, varStr, staX, staY, staZ, stationNames, zRelativeToSurf=False ) :

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dcs = self.extract( stacks, var, staX, staY, staZ, stationNames,
                            zRelativeToSurf )
        compDCs.append( dcs )
      # merge components
      for i in range(len(compDCs[0])) :
        compDCs[0][i].mergeFields( compDCs[1][i] )
        compDCs[0][i].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
        compDCs[0][i].setMetaData( 'variable',varStr )
      return compDCs[0]

    staX = np.array(staX)
    staY = np.array(staY)
    staZ = np.array(staZ)
    actualZ = np.ones_like(staZ)*np.nan
    iTri = u = nix = staList = None
    nStacks = len(stacks)
    times=values=None
    for i in range(len(stacks)) :
      stack = stacks[i]
      f = os.path.join( self.path, str(stack)+'_'+self.fileTypeStr )
      if not os.path.isfile(f) :
        print 'file not found, skipping',f
        continue
      ncfile = NetCDFFile(f,'r')
      if not self.headerIsRead :
        self.readHeader(ncfile)
      if iTri == None :
        # find parent triangles
        iTri, u, nix, goodStaIx = findStations( self.tree, self.node2elem,
                                          self.nodeX, self.nodeY, self.faceNodes,
                                          staX, staY, stationNames )
        # discard information on bad stations
        stationNames = list(np.array(stationNames)[goodStaIx])
        staX = staX[goodStaIx]
        staY = staY[goodStaIx]
        staZ = staZ[goodStaIx]

      nSta = len(staX)

      if i == 0 :
        times = np.ma.ones( (self.nTime*nStacks,) )*np.nan
        values = np.ma.ones( (self.nTime*nStacks,nSta) )*np.nan

      try :
        ti,vi,zi = self.getVProf( ncfile, varStr, staX, staY, iTri, u, nix )
      except Exception as e :
        print 'Extraction failed',e
        traceback.print_exc(file=sys.stdout)
        continue
      nTime = vi.shape[0]
      nVrt = vi.shape[1]
      is3d = nVrt > 1
      vals = np.ma.ones( (self.nTime,nSta) )*np.nan
      if not is3d :
        try :
          vals[:,:] = vi[:,0,:]
        except Exception as e :
          print 'Extraction failed',e
          traceback.print_exc(file=sys.stdout)
          continue
      else :
        for iSta in range(nSta) :
          # interpolate in vertical
          # vi,zi is (nTime,nVert)
          for iTime in range(nTime) :
            ZZ = zi[iTime,:,iSta]
            VV = vi[iTime,:,iSta]
            goodIx = np.logical_and( ~ZZ.mask, ~VV.mask )
            ZZ = ZZ[goodIx]
            VV = VV[goodIx]
            if len(ZZ)==0 or len(VV)==0 :
              continue
            # TODO could be optimized with interp2d ??
            z = staZ[iSta]
            if zRelativeToSurf :
              # staZ treated as depth below surface
              z =ZZ[-1] - staZ[iSta]
            if z < ZZ[0] or z > ZZ[-1] :
              # z out of bounds
              if z < ZZ[0] :
                z = ZZ[0]
              if z > ZZ[-1] :
                z = ZZ[-1]
              staStr = stationNames[iSta] if stationNames else ''
              actualZ[iSta] = z
            try :
              vals[iTime,iSta] = interp1d( ZZ, VV ) ( z )
            except Exception as e :
              print 'Vertical interpolation failed',e
              print ZZ.shape, VV.shape
              traceback.print_exc(file=sys.stdout)
              continue
          if not np.isnan( actualZ[iSta] ) :
            print 'z coordinate changed: ',staStr,staZ[iSta],actualZ[iSta]
      try :
        # ti (nTime,)
        times[i*self.nTime:(i+1)*self.nTime] = ti
        # vals (nTime,nSta)
        values[i*self.nTime:(i+1)*self.nTime,:] = vals
      except Exception as e :
        print 'Extraction failed',e
        traceback.print_exc(file=sys.stdout)
        continue

    # remove skipped stacks
    goodIx = np.logical_not( np.isnan( times ) )

    times = np.ma.masked_invalid( times[goodIx] )
    values = np.ma.masked_invalid( values[goodIx,:])

    # build dataContainer for each station
    dcs = []
    for iSta in range(len(staX)) :
      data = values[:,iSta]
      # remove nans
      goodIx = np.logical_not( data.mask )
      if not goodIx.any() :
        # all bad data
        print 'all bad data',stationNames[iSta]
        continue
      # NOTE data,time must be ndarray not masked array
      data = np.reshape( np.array(data[goodIx]), (1,1,-1) )
      t = np.array(times[goodIx])

      ta = timeArray( t, 'epoch' )
      zSign = 1 if zRelativeToSurf else -1 # zRelToSurf => depth below surface
      msldepth = str(int(round(zSign*staZ[iSta]*100)))
      meta = {}
      meta['location'] = stationNames[iSta]
      meta['instrument'] = 'model'
      meta['variable'] = varStr
      meta['bracket'] = 'F' if zRelativeToSurf else 'A'
      meta['msldepth'] = msldepth
      meta['dataType'] = 'timeseries'
      z = actualZ[iSta] if not np.isnan(actualZ[iSta]) else staZ[iSta]
      x = staX[iSta]
      y = staY[iSta]
      dc = dataContainer('', ta, x,y,z, data, fieldNameList.get(varStr,[varStr]),
                          coordSys='spcs',metaData=meta)
      dcs.append(dc)
    return dcs

class ncExtractProfile(ncExtractBase) :
  """A class for extracting vertical profile time series from SELFE netCDF4
  output files."""

  def extract(self,stacks,varStr,staX,staY,stationNames) :

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dcs = self.extract( stacks, var, staX, staY, stationNames)
        compDCs.append( dcs )
      # merge components
      for i in range(len(compDCs[0])) :
        compDCs[0][i].mergeFields( compDCs[1][i] )
        compDCs[0][i].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
        compDCs[0][i].setMetaData( 'variable',varStr )
      return compDCs[0]

    staX = np.array(staX)
    staY = np.array(staY)
    time,vals,zcor,goodStaIx = self.getProfForStacks(stacks,varStr,
                                                      staX,staY,stationNames)
    # discard information on bad stations
    stationNames = list(np.array(stationNames)[goodStaIx])
    staX = staX[goodStaIx]
    staY = staY[goodStaIx]
    # build dataContainer for each station
    dcs = []
    for iSta in range(len(staX)) :
      # remove time steps with all bad values
      goodIxTime = np.logical_and(
                    np.logical_not( np.all( vals[:,:,iSta].mask, axis=1 ) ),
                    np.logical_not( np.all( zcor[:,:,iSta].mask, axis=1 ) ) )
      v = vals[goodIxTime,:,iSta]
      z = zcor[goodIxTime,:,iSta]
      t = time[goodIxTime]
      # remove masked (below bottom) part (nTime,nVert) -> (nTime,nGoodVert)
      goodIxVert = np.logical_and( np.logical_not( np.any( v.mask, axis=0 ) ),
                                   np.logical_not( np.any( z.mask, axis=0 ) ) )
      if not goodIxVert.any() or not goodIxTime.any() :
        # all bad data
        print not goodIxVert.any(), not goodIxTime.any()
        print 'all bad data',stationNames[iSta]
        continue
      v = v[:,goodIxVert]
      z = z[:,goodIxVert]
      if v.mask.any() or z.mask.any() or t.mask.any() :
        print v.mask.any(), z.mask.any(), t.mask.any()
        print v.mask
        print z.mask
        print t.mask
        raise Exception('bad values remain '+stationNames[iSta])
      # to (nGoodVert,1,nTime)
      data = np.swapaxes(np.array(v),0,1)[:,None,:]
      ta = timeArray( np.array(t), 'epoch' )
      # to (nGoodVert,nTime)
      z = np.swapaxes(np.array(z),0,1)
      nZ = z.shape[0]
      x = staX[iSta]*np.ones((nZ,))
      y = staY[iSta]*np.ones((nZ,))
      meta = {}
      meta['location'] = stationNames[iSta]
      meta['instrument'] = 'model'
      meta['bracket'] = 'A'
      meta['variable'] = varStr
      meta['dataType'] = 'profile'
      dc = dataContainer('', ta, x,y,z, data, fieldNameList.get(varStr,[varStr]),
                          coordSys='spcs',metaData=meta)
      dcs.append(dc)
    return dcs

def extractTransectForCoords( x, y, dataDir, varList, startTime, endTime, name, modelCoordSys='spcs' ) :
  """Extracts a transect defined by x,y coordinates and stores to disk in netCDF format."""
  dcs = []
  for var in varList :
    try :
      ee = ncExtractTransect(dataDir, var=var)
      dc = ee.extractDates( startTime, endTime, var, x,y,name )
      dcs.append( dc )
    except Exception as e :
      print 'Extraction failed'
      traceback.print_exc(file=sys.stdout)
  return dcs

def extractTransectForBPFile( bpFile, dataDir, varList, startTime, endTime,
                             name, modelCoordSys=None ) :
  """Extracts a transect defined in the bpFile and stores to disk in netCDF format."""
  bpObj = BuildPoint()
  bpObj.readFileFromDisk(bpFile)
  x = bpObj.points[:,1]
  y = bpObj.points[:,2]

  return extractTransectForCoords( x, y, dataDir, varList, startTime, endTime,
                                   name )

class ncExtractTransect(ncExtractBase) :
  """A class for extracting transect time series from SELFE netCDF4
  output files."""

  def extract(self,stacks,varStr,staX,staY,transName) :

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dc = self.extract( stacks, var, staX, staY, transName)
        compDCs.append( dc )
      # merge components
      compDCs[0].mergeFields( compDCs[1] )
      compDCs[0].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
      compDCs[0].setMetaData('variable',varStr)
      return compDCs[0]

    staX = np.array(staX)
    staY = np.array(staY)
    time,vals,zcor,goodStaIx = self.getProfForStacks(stacks,varStr,staX,staY)
    # discard information on bad stations
    staX = staX[goodStaIx]
    staY = staY[goodStaIx]

      #times -- (masked array (nTime,)) Concatenated time stamps
      #values -- (masked array (nTime,nVert,nSta)) extracted values
      #zcoords -- (masked array (nTime,nVert,nSta)) z coordinates
      #goodStaIx -- (masked array (nOKSta,)) indices of good stations (inside the grid)

    # reorganize transect in track format
    data = []
    X = []
    Y = []
    Z = []
    for iSta in range(len(staX)) :
      z = zcor[:,:,iSta]
      v = vals[:,:,iSta]
      # remove masked (below bottom) part (nTime,nVert) -> (nTime,nGoodVert)
      goodIxVert = np.logical_and( ~np.all( v.mask, axis=0 ),
                                   ~np.all( z.mask, axis=0 ) )
      z = z[:,goodIxVert].filled(np.nan)
      v = v[:,goodIxVert].filled(np.nan)
      
      x = staX[iSta]*np.ones((z.shape[1],))
      y = staY[iSta]*np.ones((z.shape[1],))
      data.append( v ) # (nTime,nVert)
      X.append( x ) # (nVert,)
      Y.append( y ) # (nVert,)
      Z.append( z ) # (nTime,nVert)
    # concatenate
    X = np.concatenate( tuple(X), axis=0 )
    Y = np.concatenate( tuple(Y), axis=0 )
    Z = np.concatenate( tuple(Z), axis=1 )
    data = np.concatenate( tuple(data), axis=1 )
    # reshape for dataContainer
    # from (nTime,nVert*nSta) to (nVert*nSta,1,nTime)
    data = np.swapaxes( data,0,1 )[:,None,:]
    Z = np.swapaxes( Z,0,1 )[:,:]
    
    # build dataContainer
    ta = timeArray( time, 'epoch' )
    meta = {}
    meta['location'] = transName
    meta['instrument'] = 'model'
    meta['bracket'] = 'A'
    meta['variable'] = varStr
    meta['dataType'] = 'transect'
    dc = dataContainer('', ta, X,Y,Z, data, fieldNameList.get(varStr,[varStr]),
                       coordSys='spcs', metaData=meta, acceptNaNs=True)
    return dc

def extractTrackForDataContainer( dataDir, trackDC, var, name ) :
  """Extracts a track based on x,y,z,time in the given dataContaner."""
  # sanity check that trackDC is suitable
  if not trackDC.isTrack() :
    raise Exception( 'given dataContainer does not contain track information' )
  # track metadata
  bracket = trackDC.getMetaData('bracket')
  if not name :
    name = trackDC.getMetaData('location')
  if not var :
    var = trackDC.getMetaData('variable')
  x = trackDC.x.flatten()
  y = trackDC.y.flatten()
  z = trackDC.z.flatten()
  # expand scalar coordinates to time-dep track
  nx = max( max( len(x), len(y) ), len(z) )
  if len(x) == 1 : x = np.tile( x, (nx,))
  if len(y) == 1 : y = np.tile( y, (nx,))
  if len(z) == 1 : z = np.tile( z, (nx,))
  zRelativeToSurf = True if bracket == 'F' else False

  print ' * extracting track', trackDC.getMetaData('location'), var
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
      ee = ncExtractTrack(dataDir, var=var)
      dc = ee.extract(var,xx,yy,zz,tt,name,zRelativeToSurf)
      if dc :
        dcs.append(dc)
    except Exception as e :
      print 'Extraction failed:'
      traceback.print_exc(file=sys.stdout)
  # merge chunks
  dc = dcs[0]
  for i in range(1,len(dcs)) :
    dc.mergeTemporal( dcs[i] )

  return dc

class ncExtractTrack(ncExtractBase) :
  """A class for extracting x,y,z,t track from SELFE netCDF4
  output files."""

  def extract(self,varStr,X,Y,Z,T,trackName,zRelativeToSurf=False) :
    """Extracts track defined by X,Y,Z,T arrays for given stacks
       and variables."""
    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dc = self.extract(var,X,Y,Z,T,trackName,zRelativeToSurf)
        compDCs.append( dc )
      # merge components
      compDCs[0].mergeFields( compDCs[1] )
      compDCs[0].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
      compDCs[0].setMetaData('variable',varStr)
      return compDCs[0]

    # figure out correct stacks
    stacks = self.getStacks( epochToDatetime( T[0] ),  epochToDatetime( T[-1] ), exactBounds=True )
    # get profiles for all x,y locations
    time,vals,zcor,goodStaIx = self.getProfForStacks(stacks,varStr,X,Y)
    # exclude points outside the grid
    X = X[goodStaIx]
    Y = Y[goodStaIx]
    Z = Z[goodStaIx]
    T = T[goodStaIx]

    if time[0] > T[0] or time[-1] < T[-1] :
      print 'Extracted time range does not cover query range: cropping'
      if time[0] > T[0] :
        print 'start',epochToDatetime(time[0]), '>', epochToDatetime(T[0])
      if time[-1] < T[-1] :
        print 'end',epochToDatetime(time[-1]), '<', epochToDatetime(T[-1])
      goodTimeIx = np.logical_and( T >= time[0], T <= time[-1] )
      X = X[goodTimeIx]
      Y = Y[goodTimeIx]
      Z = Z[goodTimeIx]
      T = T[goodTimeIx]

    # array for exported data
    data = np.zeros_like(X)

    actualZ = np.ones_like(Z)*np.nan
    nTime = vals.shape[0]
    nP = len(X)
    # interpolate in vertical
    for i in range(nP) :
      # interpolate in time first (z depends on time)
      # z coordinates and values at X[i],Y[i],T[i]
      #print 'bad z',len(np.nonzero(~np.isfinite(zcor))[0])
      ZZ = np.swapaxes(zcor[:,:,i],0,1)
      VV = np.swapaxes(vals[:,:,i],0,1)
      goodIx = ~np.all(ZZ.mask,axis=1) # remove mask
      goodIx = np.logical_and( goodIx, ~np.all(np.isnan(VV),axis=1) )
      if not goodIx.any() :
        # all bad data
        continue
      ztmp = interp1d( time,ZZ[goodIx,:] )( T[i] )
      vtmp = interp1d( time,VV[goodIx,:] )( T[i] )
      # interpolate in vertical
      z = Z[i]
      if zRelativeToSurf :
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

    # create dataContainer
    ta = timeArray(T,'epoch')
    data = data[None,None,:]

    meta = {}
    meta['dataType'] = 'track'
    meta['location'] = trackName
    meta['instrument'] = 'model'
    meta['bracket'] = 'F' if zRelativeToSurf else 'A'
    meta['variable'] = varStr
    dc = dataContainer('', ta, X[None,:],Y[None,:],Z[None,:], data,
                       fieldNameList.get(varStr,[varStr]),
                       coordSys='spcs',metaData=meta,acceptNaNs=True)

    return dc

def extractSlabForLevel( dataDir, varList, startTime, endTime, name, zCoord=None, kLevel=None, zRelativeToSurf=False ) :

  mcs = []
  for var in varList :
    ee = ncExtractSlab(dataDir, var=var)
    mc = ee.extractDates( startTime, endTime, var, name,
                          zCoord, kLevel, zRelativeToSurf )
    mcs.append( mc )
  return mcs

class ncExtractSlab(ncExtractBase) :
  """A class for extracting slab time series from SELFE netCDF4
  output files."""

  def extract(self,stacks,varStr,slabName,zCoord=None, kLevel=None,
              zRelativeToSurf=False,bBox=None) :

    if varStr in vectorVars :
      # recursion: if a vector field requested, extract components separately
      varList = [varStr+'_x', varStr+'_y'] # hvel_x, hvel_y
      compDCs = []
      for var in varList :
        dc = self.extract( stacks, var, slabName, zCoord, kLevel,
                            zRelativeToSurf, bBox )
        compDCs.append( dc )
      # merge components
      compDCs[0].mergeFields( compDCs[1] )
      compDCs[0].fieldNames = fieldNameList.get(varStr,[varStr]) # hvel: ['u','v']
      compDCs[0].setMetaData('variable',varStr)
      return compDCs[0]

    # (nTime, ) (nTime,nNodes)
    t,v = self.getSlabForStacks(stacks, varStr, zCoord, kLevel,
                                zRelativeToSurf, iTri=None)
    alongSLevel = zCoord==None
    # build meshContainer
    ta = timeArray(t, 'epoch')
    x = self.nodeX
    y = self.nodeY
    connectivity = self.faceNodes
    data = np.swapaxes(v,0,1)[:,None,:]
    # replace masked values with nans
    data = data.filled(np.nan)
    msldepth = ''
    if alongSLevel :
      msldepth = 'slev'+str(kLevel)
      z = 0.0*np.ones_like(x)
    else :
      zSign = 1 if zRelativeToSurf else -1 # zRelToSurf => depth below surface
      msldepth = str(int(round(zSign*zCoord*100)))
      z = zCoord*np.ones_like(x)
    meta = {}
    meta['dataType'] = 'slab'
    meta['location'] = slabName
    meta['instrument'] = 'model'
    meta['variable'] = varStr
    if alongSLevel :
      meta['slevel'] = kLevel
    else :
      meta['bracket'] = 'F' if zRelativeToSurf else 'A'
      meta['msldepth'] = msldepth
    mc = meshContainer('', ta, x,y,z, data, connectivity, fieldNameList.get(varStr,[varStr]), coordSys='spcs',metaData=meta)
    return mc

def test() :
  # test for nc extration
  #path = '/home/tuomas/workspace/cmop/projects/netcdf_outputs/run144/outputs/ncfiles/'
  path = '/home/tuomas/temp/db29skill/db29-2012/combined_nc'
  #fname = 'elev.nc'
  fname = 'salt.63.nc'
  #fname = 'hvel.64.nc'
  #sT = datetime.datetime(2012,5,8)
  #eT = datetime.datetime(2012,5,10)
  sT = datetime.datetime(2012,11,9)
  eT = datetime.datetime(2012,11,12)

  staX = [288458.4, 344208, 800000, 380769.3, 358500, 337022]
  staY = [416759.7, 287180, 800000, 289948.7, 285529, 294843]
  staZ = [-1.5, -1.5, -1.5, -1.5, -1.5, -1.5]
  stationNames = ['nbd41','saturn03','fake','tnslh', 'foreverdry', 'jetta']

  #varStr = 'elev'
  varStr = 'salt'
  #varStr = 'hvel'

  # test time series
  se = ncExtractStation(path,fname)
  dcs = se.extractDates(sT,eT, varStr, staX, staY, staZ, stationNames, zRelativeToSurf=False)
  for dc in dcs :
    print ' *** '
    print dc
  from plotting.timeSeriesPlot import timeSeriesPlotDC, plt
  kwargs = {}
  kwargs['title'] = 'timeSeriesPlotDC example'
  kwargs['ylim'] = [-0.5,34]
  dia = timeSeriesPlotDC( 'Elevation', unit='m')
  for dc in dcs :
    dia.addSample(dc, label=dc.getMetaData('location'), linestyle='-')
    #dia.addSample(dc.extractField('v'), label=dc.getMetaData('location'), linestyle='-')
  dia.makePlot(**kwargs)
  plt.show()

  # test profile
  pe = ncExtractProfile(path,fname)
  dcs = pe.extractDates(sT,eT, varStr, staX, staY, stationNames)
  for dc in dcs :
    print ' *** '
    print dc
  from plotting.profilePlot import stackProfileTimeSeriesDC, plt
  dia = stackProfileTimeSeriesDC(clabel=varStr,unit='psu',clim=[0,34])
  for dc in dcs :
    dia.addSample( dc.getMetaData('location'), dc, title='test' )
    dia.addTitle( dc.getMetaData('location'), tag=dc.getMetaData('location') )
  dia.showColorBar()
  plt.show()

  # test transect
  dcs = extractTransectForBPFile('nchannel_fine.bp', path, [varStr], sT, eT, 'nchannel_fine')
  dc = dcs[0]
  print dc
  from plotting.transectPlot import stackTransectPlotDC, plt
  dia = stackTransectPlotDC()
  dia.addSample( '1', dc,0 , N=35, clabel=varStr,unit='unit')#, clim=[0,34] )
  dia.addSample( '2', dc,10, N=35, clabel=varStr,unit='unit')#, clim=[0,34] )
  dia.addSample( '3', dc,20, N=35, clabel=varStr,unit='unit')#, clim=[0,34] )
  dia.addTitle('stackTransectPlotDC example')
  plt.show()

  # test slab
  dcs = extractSlabForLevel(path,[varStr],sT,eT,'slab',kLevel=-1)
  dc = dcs[0]
  print dc
  dc = dc.cropGrid( [326000,352000,284000,302000] )

  from plotting.slabPlot import slabSnapshotDC, plt
  fig = plt.figure()
  dia = slabSnapshotDC(clabel=varStr,unit='unit')
  dia.setupAxes( fig )
  dia.addSample( dc, 20 )
  dia.showColorBar()
  plt.show()

if __name__ == '__main__' :
  test()
