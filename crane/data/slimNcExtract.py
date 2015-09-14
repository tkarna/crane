
import os
import sys
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from scipy.spatial import KDTree, cKDTree
from scipy.interpolate import interp1d
from glob import glob
import time as timeMod

from data.timeArray import *
from data.dataContainer import *
from data.meshContainer import meshContainer
from data.extractStation import fieldNameList
from data.selfeGridUtils import *
from files.stationFile import StationFile
from files.buildPoints import BuildPoint

def constructCGConnectivity(x,y,connectivity,kdtree) :
  """Constructs connectivity table by collapsing overlapping nodes
  If the given connectivity table is discontinuous, will return equivalent
  continuous grid."""
  MAX_NB_SIBLINGS = 15 # max nb element touching the same node
  TOL = 1e-3 # distance tolerance for collapsing nodes
  # find n nearest neighbors for each node
  d,i = cElem = kdtree.query( np.concatenate( (x[:,None],y[:,None]), axis=1 ), k=MAX_NB_SIBLINGS )
  # discard nodes outside distance tolerance
  i[d > TOL] = 1e10
  # for each list of overlapping nodes, take the minimum as temp index
  nodeToCGNode = np.min( i, axis=1 )
  # assing new indices ranging from 0 to nb_CG_nodes
  uni, inv = np.unique(nodeToCGNode, return_inverse=True)
  newIx = np.arange(len(uni))
  nodeToCGNode = newIx[inv]
  connectivityCG = nodeToCGNode[connectivity]
  cgNodeToNode = np.unique(nodeToCGNode,return_index=True)[1]
  return connectivityCG, nodeToCGNode, cgNodeToNode

def constructNodeToElemTable(conn) :
  """Constructs node to element table for the given connectivity array

  Parameters
  ----------
    conn : (nElems,3) array
           mesh connectivity array
  Returns:
    node2elem : (dict of lists)
                a dictionary that maps node index to list of elements that share that node
"""
  node2elem = {}
  nElems = conn.shape[0]
  for iE in range(nElems) :
    for iN in conn[iE,:] :
      node2elem.setdefault(iN,[]).append( iE )
  return node2elem

def constructElemToNeighborTable(conn,node2elem) :
  """Constructs element to neighbor element table for the given connectivity array

  Parameters
  ----------
    conn : (nElems,3) array
           mesh connectivity array
  Returns:
    elem2neigh : (nElems, 3) array
               table that maps each element to its 3 neighbors
"""
  elem2neigh = np.ones_like(conn,dtype=int)*-1
  nElems, nNodes = conn.shape
  for iE in range(nElems) :
    for iEdge in range(nNodes) :
      n1 = conn[iE,iEdge]
      n2 = conn[iE,(iEdge+1)%nNodes]
      # all possible neighbors
      ix = list(node2elem[n1])
      ix += node2elem[n2]
      # find elements that appear twice
      multiplicity = {}
      for i in ix :
        multiplicity.setdefault(i,[]).append(i)
      iNeigh = [ i for i in multiplicity if len(multiplicity[i]) == 2 and i!=iE ]
      if len(iNeigh) > 0 :
        elem2neigh[iE,iEdge]=iNeigh[0]
    if elem2neigh[iE,:].max() == -1 :
      print 'Warning: element with no neighbors: ',i
  return elem2neigh

def isInsideTri(x,y,u=None,nodeX=None,nodeY=None,) :
  """ Tests whether (x,y) is in a triangle  whose vertices are nodeX,nodeY.
  Parameters
  ----------
    NodeX,nodeY -- array (nTri,3)
    x,y         -- scalar
  Returns:
    result      -- bool array (nTri,)
  """
  if u==None :
    u = barycentricCoords(nodeX,nodeY,x,y)
  return  np.logical_and( u>=0 , u<=1 ).all( axis=1 )

def interpolateTri(u,nodalValues) :
  """Interpolates nodal values in a location given by barycentric coordinates u
  for all the elements.
  
  Parameters
  ----------
    nodalValues : array_like (nTri,3)
        Nodal values of the field at the triangle vertices
    u           : array_like (1,3) or (nTri,3)
        Element barycentric coordinates for the interpolation point
  Returns
  -------
    values      : array_like (nTri,1)
        Interpolated values
  """
  nTri = nodalValues.shape[0]
  values = np.ma.ones( (nTri,) )*np.nan
  nv = np.ma.masked_invalid( nodalValues )
  goodIx =  np.all( np.logical_not(nv.mask), axis=1 )
  if not goodIx.any() :
    return np.ma.masked_invalid( values[:,None] )
  nv = nv[goodIx,:]
  if u.shape[0] == 1 and nodalValues.shape[0] > 1 :
    U = np.tile(u,(nv.shape[0],1))
  else :
    U = u
  values[goodIx,:] = np.sum(u*nv,axis=1)
  return np.ma.masked_invalid( values[:,None] )

def barycentricCoords( nodeX, nodeY, x, y ) :
  """Returns barycentric coordinates for (x,y) in triangle whose vertices are
  nodeX,nodeY.
  Parameters
    NodeX,nodeY -- array (nTri,3)
    x,y         -- scalar
  Returns:
    u           -- array (nTri,3)
  """
  detT = ( (nodeY[:,1]-nodeY[:,2])*(nodeX[:,0]-nodeX[:,2]) -
           (nodeX[:,1]-nodeX[:,2])*(nodeY[:,0]-nodeY[:,2]) )
  u1 = ( (nodeY[:,1]-nodeY[:,2])*((x-nodeX[:,2])) - (nodeX[:,1]-nodeX[:,2])*(y-nodeY[:,2]) ) / detT
  u2 = ( (nodeY[:,2]-nodeY[:,0])*((x-nodeX[:,2])) - (nodeX[:,2]-nodeX[:,0])*(y-nodeY[:,2]) ) / detT
  u3 = 1.0 - u1 - u2
  u = np.array([u1,u2,u3]).T
  return u

class meshSearch2d(object) :
  """A class for finding mesh elements that contain query points"""
  def __init__(self, node_x, node_y, connectivity) :
    """Initialize search tree"""
    self.node_x = node_x
    self.node_y = node_y
    self.connectivity = connectivity
    # build 2d kd-tree
    self.tree = cKDTree( np.concatenate( (node_x[:,None],node_y[:,None]), axis=1 ) )
    a,b,c = constructCGConnectivity(self.node_x,self.node_y,self.connectivity,self.tree)
    self.connectivityCG = a
    self.nodeToCGNode = b
    self.cgNodeToNode = c
    # build node to element table
    self.node2elem = constructNodeToElemTable(self.connectivityCG)
    # build elem to neighbor table
    self.elem2neigh = constructElemToNeighborTable(self.connectivityCG,self.node2elem)
    #print self.node2elem[1]
    #print self.elem2neigh[:10,:]

    #import matplotlib.pyplot as plt
    #import matplotlib.tri as tri
    #print self.cgNodeToNode.shape
    #xx = self.node_x[self.cgNodeToNode]
    #yy = self.node_y[self.cgNodeToNode]
    #tri = tri.Triangulation(xx, yy, self.connectivityCG)
    #plt.triplot(tri, 'b-')
    #plt.show()
    #exit(0)

    # compute element centers
    self.center_x = np.mean(self.node_x[connectivity],axis=1)
    self.center_y = np.mean(self.node_y[connectivity],axis=1)
    self.tree_center = cKDTree( np.concatenate( (self.center_x[:,None],self.center_y[:,None]), axis=1 ) )
    
  def findParentElement(self, x, y) :
    """
    Returns element that contains the given point (x,y).
    
    Parameters:
    ----------
    x,y  : float
        coordinates of the query point

    Returns:
    --------
    parentId : int
        Index of the parent element that cointains the query point.
        Index is -1 if no parent element is found.
    """
    # find nearest element center
    cElem = self.tree_center.query( np.array([[x,y]]) )[1][0]
    # list of all possible elements: center elem + its neighbors
    elems = self.elem2neigh[cElem,:]
    elems = elems[elems>=0] # discard -1s
    elems = np.hstack(([cElem],elems))
    # test if inside an element
    e_x = self.node_x[self.connectivity[elems],:]
    e_y = self.node_y[self.connectivity[elems],:]
    goodElems = isInsideTri(x,y,nodeX=e_x,nodeY=e_y)
    if np.any(goodElems) :
      parent = elems[goodElems][0]
    else :
      parent = -1
    return parent
    
  def findParentElementMultiple(self, x, y) :
    """
    Returns elements that contain the given points (x,y).
    
    Parameters:
    ----------
    x,y  : array_like (nPoints,)
        coordinates of the query points

    Returns:
    --------
    parentId : array_like (nPoints,)
        Indices of the parent elements that cointains the query points.
        Index is set to -1 if no parent element is found.
    """
    # find nearest element center
    query = np.vstack((x,y)).T
    cElem = self.tree_center.query( query )[1]
    parents = np.ones_like(x,dtype=int)*-1
    for i in range(len(x)) :
      # list of all possible elements: center elem + its neighbors
      elems = self.elem2neigh[cElem[i],:]
      elems = np.hstack(([cElem[i]],elems))
      # test if inside an element
      e_x = self.node_x[self.connectivity[elems],:]
      e_y = self.node_y[self.connectivity[elems],:]
      goodElems = isInsideTri(x[i],y[i],nodeX=e_x,nodeY=e_y)
      if np.any(goodElems) :
        parents[i] = elems[goodElems][0]
    return parents    

  def evaluateOnMesh(self, x, y, nodalValues, fillNaNsWithNearest = False):
    """
    Evaluates a field, defined by the nodal values, on the mesh at locations (x,y)
    
    Parameters:
    ----------
    x,y  : array_like (nPoints,)
        coordinates of the query points
    nodalValues : array_like (nNodes,)
        values of the field to evaluate at the mesh nodes
    
    Returns:
    --------
    vals : array_like (nPoints)
        Interpolated values at the points (x,y)

    """
    parentElems = self.findParentElementMultiple(x,y)
    goodElems = parentElems >= 0
    parentElems = parentElems[goodElems]
    xx = x[goodElems]
    yy = y[goodElems]
    
    vals = np.ones_like(x)*np.nan
    
    elemX = self.node_x[self.connectivity[parentElems,:]]
    elemY = self.node_y[self.connectivity[parentElems,:]]
    elemVals = nodalValues[self.connectivity[parentElems,:]]
    u = barycentricCoords( elemX, elemY, xx, yy )
    vals[goodElems] = interpolateTri( u, elemVals )

    if fillNaNsWithNearest :
      # fill gaps with nearest node value
      # assuming corrupted nodes are the same for all components
      badIx = np.nonzero(~np.isfinite( vals ))[0]
      # nearest nodes
      nix = self.tree.query( np.vstack((x[badIx],y[badIx])).T )[1]
      vals[badIx,:] = nodalValues[nix]

    return vals
  
class horizontalInterpolator(object) :
  """
  A class that stores horizontal interpolation information for quickly
  interpolating over multiple fields
  """
  def __init__(self, meshSearch2d, x, y, stationNames=None):
    """Initialized the interpolation routine for the given (x,y) points"""
    self.meshSearch2d = meshSearch2d
    parentElems = self.meshSearch2d.findParentElementMultiple(x,y)
    self.goodElems =  parentElems >= 0
    self.parentElems = parentElems[self.goodElems]
    self.x = x
    self.y = y
    self.x_good = x[self.goodElems]
    self.y_good = y[self.goodElems]
    self.elem_x = self.meshSearch2d.node_x[self.meshSearch2d.connectivity[self.parentElems,:]]
    self.elem_y = self.meshSearch2d.node_y[self.meshSearch2d.connectivity[self.parentElems,:]]
    self.u = barycentricCoords( self.elem_x, self.elem_y, self.x_good, self.y_good )
    if not np.all(self.goodElems) :
      badIx = np.nonzero(~self.goodElems)[0]
      print 'Warning: Point(s) out of grid'
      for ix in badIx :
        if stationNames :
          print '  ', stationNames[ix], self.x[ix],self.y[ix]
        else :
          print '  ', self.x[ix],self.y[ix]

  def evaluate(self, nodalValues) :
    """Interpolates the field defined by the nodal values
    
    Parameters
    ----------
    nodalValues : array_like (nNodes,)
                  Nodal values of a field defined on the mesh
                  
    Returns
    -------
    values      : array_like (nPoints,)
                  Interpolated values
    """
    elemVals = nodalValues[self.meshSearch2d.connectivity[self.parentElems,:]]
    vals = np.ones_like(self.x)*np.nan
    vals[self.goodElems] = interpolateTri( self.u, elemVals )
    return vals
  
  def evaluateArray(self, nodalValues) :
    """Interpolates the field defined by the nodal values
    
    Parameters
    ----------
    nodalValues : array_like (nNodes,nDim2,nDim3)
                  Nodal values of a field defined on the mesh for multiple
                  depths and time instances

    Returns
    -------
    values      : array_like (nPoints,nDim2,nDim3)
                  Interpolated values
    """
    nPoints = len(self.x)
    nNodes, nDim2, nDim3 = nodalValues.shape
    vals = np.ones((nPoints,nDim2,nDim3))
    for i in xrange(nDim2):
      for j in xrange(nDim3):
        vals[:  ,i,j] = self.evaluate(nodalValues[:,i,j])
    
    return vals

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
    
  Returns
  -------
  vals : array_like (nTime,)
       Interpolated values
  z_target : array_like (nTime,)
       The z coordinate at which the interpolation actually took place
  """
  if k==None and z==None :
    raise Exception('Either k or z must be defined')
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
  vals = np.zeros((nTime,))
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
        vals[iTime] = interp1d( ZZ[ix,iTime], VV[ix,iTime], kind='linear', copy=False ) ( z_target[iTime] )
    else :
      for iTime in xrange(nTime) :
        vals[iTime] = interp1d( ZZ[:,iTime], VV[:,iTime], kind='linear', copy=False ) ( z_target[iTime] )
  if k != None :
    # bottom: k=1 kk=0, surface: k=-1 kk=len(z)-1
    kk = k-1 if k>0 else nZ+k
    vals = VV[kk,:]
    z_target = ZZ[kk,:]
  return vals,z_target

class verticalInterpolator(object) :
  """A class that handles interpolates data in vertical."""
  def __init__(self,z=None,k=None,zRelToSurf=None):
    """
    Parameters
    ----------
    z : float, array_like (nProfiles,), optional
      z coordinate where each vertical profile is evaluated. z coordinates
      increase upwards.
    k : int, array_like (nProfiles,), optional
      Instead of interpolating, take k-th nodal value from bottom.
      k=1 stands for bottom, k=-1 stands for surface
    zRelToSurf : bool, array_like (nProfiles,), optional
      If True z coordinate is taken depth below free surface instead of
      static z coordinate
      
    """
    if k==None and z==None :
      raise Exception('Either k or z must be defined')
    self.doZInterpolation = k == None
    self.z = z
    self.k = k
    if self.doZInterpolation and zRelToSurf == None :
      zRelToSurf = np.zeros_like(self.z,dtype=bool)
    self.zRelToSurf = zRelToSurf
  
  def evaluateArray(self, Z,V ):
    """
    Parameters
    ----------
    Z : ndarray (nProfiles,nZ,nTime)
      z coordinates for each vertical profile. If nTime>1 each profile will be
      interpolated separately.
    V : ndarray (nProfiles,nZ,nTime)
      values corresponding to each z point

    Returns
    -------
    vals : array_like (nProfiles,nTime,)
        Interpolated values
    z_actual : array_like (nProfiles,nTime,)
        The z coordinate at which the interpolation actually took place
    """
    nProfiles, nZ, nTime = Z.shape
    vals = np.ones((nProfiles,nTime))*np.nan
    z_actual = np.ones((nProfiles,nTime))*np.nan
    
    goodProfiles = ~np.all(np.all(~np.isfinite(V),axis=2),axis=1)
    goodProfiles = np.nonzero(goodProfiles)[0]
    for i in goodProfiles:
      if self.doZInterpolation :
        vi,zti = interpolateInVertical(Z[i,:,:],V[i,:,:], z=self.z[i],
                                       zRelToSurf=self.zRelToSurf[i])
      else :
        vi,zti = interpolateInVertical(Z[i,:,:],V[i,:,:], k=self.k[i])
      vals[i,:] = vi
      z_actual[i,:] = zti
    return vals,z_actual

# use consistent field names throughout the skill assessment package

VARS2D = ['elev','dahv']

def getNCVariableName(varStr) :
  nc_names = { 'elev':'eta',
               'salt':'S',
              }
  return nc_names.get(varStr,varStr)

class slimExtractBase(object) :
  """Base class for all netcdf extract objects."""
  # TODO morf into generic base class, leave all model dependent stuff undefined
  def __init__(self, path, var=None, verbose=False ) :
    """Intializes reader object."""
    fieldNameToFilename = { 'temp':'T',
                            'elev':'eta',
                            'salt':'S',
                            'kine':'',
                            'vdff':'',
                            'tdff':'',
                            'mixl':'',
                            'hvel':'',
                            'vert':'',
                            'dens':'',
                            'trcr_1':'',
                            'trcr_2':'',
                            'trcr_3':'',
                            'trcr_4':'',
                            'turbidity':''} # TODO belongs to derived class
    self.path = path
    self.component = 0
    self.fileTypeStr = fieldNameToFilename[var]
    self.headerIsRead = False
    self.verbose = verbose

  def generateFileName(self, iStack=None, fileTypeStr=None) :
    """Returns full path to the netcdf file for iStack.
    If iStack==None, returns a pattern with '*' as a wildcard."""
    # TODO raise NotImplementedError('This method must be defined in the derived class')
    if fileTypeStr == None : fileTypeStr = self.fileTypeStr
    stackStr = '*' if iStack == None else '{0:05d}'.format(iStack) 
    fname = '{typeStr:s}_{stack:s}_COMP_{comp:d}.nc'.format(
                        typeStr=fileTypeStr,stack=stackStr,
                        comp=self.component )
    return os.path.join(self.path,fname)
    
  def getNCFile( self, iStack=None, fileTypeStr=None ) :
    """Opens netcdf file corresponding to the given stack number.
    If no stack number is given opens first matching file."""
    if fileTypeStr == None : fileTypeStr = self.fileTypeStr
    f = self.generateFileName(iStack, fileTypeStr)
    if iStack==None :
      # try to find a file that matches file name pattern
      pattern=f
      files = sorted(glob(pattern))
      if len(files) == 0 :
        raise Exception('no files found in '+pattern)
      f = files[0]
    if not os.path.isfile(f) :
      raise IOError('File not found: '+f)
    else :
      if self.verbose: print 'Opening file',f
      return NetCDFFile(f,'r')

  def readHeader(self, ncfile=None) :
    """
    Reads header of the netcdf file and prepares data structures.

    If ncfile is given, will read its header. Otherwise will search for first
    matching netcdf file in the path.
    """
    # TODO raise NotImplementedError('This method must be defined in the derived class')
    ncfileGiven = ncfile != None
    if self.verbose : 
      print 'Reading header'
    if not ncfileGiven :
      ncfile = self.getNCFile()

    # read
    faceNOffset = 0 #ncfile.variables['face_nodes'].start_index
    self.faceNodes = ncfile.variables['face_nodes'][:].astype(int) - faceNOffset
    self.nodeX = ncfile.variables['node_x'][:]
    self.nodeY = ncfile.variables['node_y'][:]

    self.nNodes = len(ncfile.dimensions['node'])
    self.nFaces = len(ncfile.dimensions['face'])
    self.nElemNodes = len(ncfile.dimensions['nFaceNodes']) # ==3 always
    self.nTime = len(ncfile.dimensions['time'])
    if 'layers' in ncfile.dimensions :
      self.nVert = len(ncfile.dimensions['layers'])
    else :
      self.nVert = 0
    if self.verbose :
      print 'nodes',self.nNodes
      print 'elems',self.nFaces
      print 'verts',self.nVert
      print 'elem nodes',self.nElemNodes
      print 'time stamps',self.nTime
    timeStr = ' '.join(ncfile.variables['time'].base_date.split()[2:4])
    self.simulationStartTime = datetime.datetime.strptime( timeStr, '%Y-%m-%d %H:%M:%S' )

    if 'node_lon' in ncfile.variables :
      self.node_lon = ncfile.variables['node_lon'][:]
    if 'node_lat' in ncfile.variables :
      self.node_lat = ncfile.variables['node_lat'][:]

    if 'edge_nodes' in ncfile.variables :
      self.edgeNodes = ncfile.variables['edge_nodes'][:]
    if 'edge_x' in ncfile.variables :
      self.edge_x = ncfile.variables['edge_x'][:]
      self.edge_y = ncfile.variables['edge_y'][:]
    if 'edge_lon' in ncfile.variables :
      self.edge_lon = ncfile.variables['edge_lon'][:]
      self.edge_lat = ncfile.variables['edge_lat'][:]

    # construct mesh search object
    self.meshSearch2d = meshSearch2d( self.nodeX, self.nodeY, self.faceNodes )
    self.headerIsRead = True
    if not ncfileGiven :
      ncfile.close()

  def getTime( self, ncfile ) :
    """Returns time stamps from given netCDF file in epoch format."""
    # TODO raise NotImplementedError('This method must be defined in the derived class')
    nTime = len(ncfile.dimensions['time'])
    startTime = ' '.join(ncfile.variables['time'].base_date.split()[2:4])
    startTime = datetime.datetime.strptime( startTime, '%Y-%m-%d %H:%M:%S' )
    time = simulationToEpochTime( ncfile.variables['time'][:], startTime )
    return time

  def getZCoordinates(self, iStack) :
    """Returns vertical coordinates for the given stack."""
    # TODO raise NotImplementedError('This method must be defined in the derived class')
    ncfile = self.getNCFile( iStack, 'z' )
    Z = ncfile.variables[getNCVariableName('z')][:]
    ncfile.close()
    return Z

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
    """

    if horzInterp==None :
      horzInterp = horizontalInterpolator( self.meshSearch2d,x,y,stationNames )

    ncfile = self.getNCFile(iStack)
    time = self.getTime(ncfile)

    V = ncfile.variables[getNCVariableName(varStr)]
    is3d = len( V.shape ) == 3
    if not is3d : # 2D variable
      if self.verbose: print '2D var', V.shape
      nodalValues = V[:][:,None,:] # expand to (nNodes,nVert,nTime)
    else :
      if self.verbose: print '3D var', V.shape
      nodalValues = V[:] # NOTE reading the whole array may be inefficient?
    ncfile.close()

    vals = horzInterp.evaluateArray( nodalValues )
    if is3d :
      Z = self.getZCoordinates(iStack)
      zcoords = horzInterp.evaluateArray( Z )
    else :
      zcoords = np.zeros_like(vals)

    vals = np.ma.masked_invalid(vals)
    zcoords = np.ma.masked_invalid(zcoords)
    return time, vals, zcoords

  def getVerticalProfileForStacks(self, stacks, varStr, x, y, stationNames=None) :
    """Extracts vertical profile for the given netcdf file stacks

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
    """
    # create interpolator object for recycling
    horzInterp = horizontalInterpolator( self.meshSearch2d,x,y,stationNames )

    time = []
    vals = []
    zcoords = []
    for stack in stacks :
      # extract for individual stacks
      try :
        ti,vi,zi = self.getVerticalProfile(stack,varStr,x,y, horzInterp=horzInterp)
        time.append(ti)
        vals.append(vi)
        zcoords.append(zi)
      except Exception as e :
        print e
    # concatenate time axis
    time = np.concatenate(tuple(time),axis=0) # (nTime,)
    vals = np.concatenate(tuple(vals),axis=2) # (nProfiles,nVert,nTime)
    zcoords = np.concatenate(tuple(zcoords),axis=2) # (nProfiles,nVert,nTime)
    time = np.ma.masked_invalid(time)
    vals = np.ma.masked_invalid(vals)
    zcoords = np.ma.masked_invalid(zcoords)
    return time,vals,zcoords

  def getTimeSeriesFromProfiles(self, vals, zcoords, z=None, k=None, zRelToSurf=None ) :
    """Interpolates vertical profiles in vertical at given depth.

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
    vertInterp = verticalInterpolator(z,k,zRelToSurf)
    v, z_actual = vertInterp.evaluateArray(zcoords,vals)
    v = np.ma.masked_invalid(v)
    return v,z_actual
  
  def getTimeSeries(self, iStack, varStr, x, y, stationNames=None, z=None, k=None, zRelToSurf=None) :
    """Extracts time series from the iStack netcdf file"""
    time, vprof, zcoords = self.getVerticalProfile(iStack, varStr, x, y, stationNames)
    if varStr in VARS2D :
      vals = vprof[:,0,:]
      vals = np.ma.masked_invalid(vals)
      z_actual = np.zeros_like(vals)
    else :
      vals, z_actual = self.getTimeSeriesFromProfiles( vprof, zcoords, z, k, zRelToSurf)
    return time,vals,z_actual
  
  def getTimeSeriesForStacks(self, stacks, varStr, x, y, stationNames=None, z=None, k=None, zRelToSurf=None) :
    """Extracts time series for the given stacks.

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
    time, vprof, zcoords = self.getVerticalProfileForStacks(stacks, varStr, x, y, stationNames)
    if varStr in VARS2D :
      vals = vprof[:,0,:]
      vals = np.ma.masked_invalid(vals)
      z_actual = np.zeros_like(vals)
    else :
      vals, z_actual = self.getTimeSeriesFromProfiles( vprof, zcoords, z, k, zRelToSurf)
    return time,vals,z_actual

  def getSlab(self, iStack, varStr, z=None, k=None, zRelToSurf=None) :
    """
    Extracts a horizontal slice from the given ncfile.
    
    Parameters
    ----------
    iStack : int
           Stack number of the netCDF file to process
    varStr : string
           Variable to extract
    z      : float, array_like (nProfiles,), optional
           z coordinate where each vertical profile is evaluated. z coordinates
           increase upwards.
    k      : int, array_like (nProfiles,), optional
           Instead of interpolating, take k-th nodal value from bottom.
           k=1 stands for bottom, k=-1 stands for surface
    zRelToSurf : bool, array_like (nProfiles,), optional
           If True z coordinate is taken depth below free surface instead of
           static z coordinate

    Returns
    -------
    time  : array_like (nTime,)
          Time stamps of the extracted data in epoch format
    vals  : array_like (nPoints, nTime)
          Values of the extracted horizontal slice.
    """

    ncfile = self.getNCFile(iStack)
    time = self.getTime(ncfile)

    V = ncfile.variables[getNCVariableName(varStr)]
    is3d = len( V.shape ) == 3
    if not is3d : # 2D variable
      if self.verbose: print '2D var', V.shape
      vals = V[:] # take
      z_actual = np.zeros_like(vals)
    else :
      if self.verbose: print '3D var', V.shape
      nodalValues = V[:] # NOTE reading the whole array may be inefficient?
      zcoords = self.getZCoordinates(iStack)
      nNodes,nZ,nTime = nodalValues.shape
      # nProfiles,nZ,nTime
      # TODO add vertical interpolation
      if k != None :
        # bottom: k=1 kk=0, surface: k=-1 kk=len(z)-1
        kk = k-1 if k>0 else nZ+k
        vals = nodalValues[:,kk,:]
        z_actual = zcoords[:,kk,:]
      else : # interpolate in vertical
        zArray = np.ones((nodalValues.shape[0],))*z
        zRelArray = np.ones((nodalValues.shape[0],),dtype=int)*zRelToSurf
        vertInterp = verticalInterpolator(zArray,None,zRelArray)
        vals, z_actual = vertInterp.evaluateArray(zcoords,nodalValues)
      vals = np.ma.masked_invalid(vals)
    ncfile.close()

    return time, vals, z_actual

  def getSlabForStacks(self, stacks, varStr, z=None, k=None, zRelToSurf=None) :
    """Returns slab for the given stacks"""
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
        print e
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

  def getStacks( self, startTime, endTime, ncfile=None, firstPointIncluded=True ) :
    """Returns a list of file stack numbers that covers the given
    time period [startTime,endTime].

    Simulation start time is read from the netcdf header.
    """
    # deduce correct stacks
    ncfileGiven = ncfile != None
    if self.verbose : 
      print 'Reading header'
    if not ncfileGiven :
      ncfile = self.getNCFile()

    if not self.headerIsRead :
      self.readHeader( ncfile )
    time = self.getTime( ncfile )
    nSpool = 1 # number of exports in each file
    exportDt = 15*60 # time interval between exports
    spoolDt = nSpool*exportDt

    startDelta = (startTime - self.simulationStartTime).total_seconds()
    endDelta = (endTime - self.simulationStartTime).total_seconds()
    if not firstPointIncluded and nSpool>1: startDelta -= exportDt
    if firstPointIncluded and nSpool>1 : endDelta += exportDt
    startStack = int(np.floor(startDelta/spoolDt)) + 1
    endStack = int(np.ceil(endDelta/spoolDt))
    if not ncfileGiven :
      ncfile.close()
    return range( startStack, endStack+1 )

class slimExtract(slimExtractBase) :
  """
  This class contains only high-level extraction routines and returns the
  data in dataContainer with metadata.
  """
  def __init__(self,path, var=None, verbose=False):
    slimExtractBase.__init__(self,path,var,verbose)
    
  def extractTimeSeries(self,startTime,endTime,var,staX,staY,stationNames, staZ=None, k=None, zRelToSurf=None):
    """Extracts time series for the given time range."""
    stacks = self.getStacks(startTime,endTime)
    time,vals,actualZ = self.getTimeSeriesForStacks(stacks,var,staX,staY,stationNames,
                                                  staZ,k,zRelToSurf)

    # build dataContainer for each station
    dcs = []
    for iSta in range(len(staX)) :
      data = vals[iSta,:]
      # remove nans
      goodIx = np.logical_not( data.mask )
      if not goodIx.any() :
        # all bad data
        print 'all bad data',stationNames[iSta]
        continue
      # NOTE data,time must be ndarray not masked array
      data = np.reshape( np.array(data[goodIx]), (1,1,-1) )
      t = np.array(time[goodIx])

      ta = timeArray( t, 'epoch' )
      meta = {}
      meta['location'] = stationNames[iSta]
      meta['instrument'] = 'model'
      meta['variable'] = var
      alongSLevel = False # FIXME
      if alongSLevel :
        meta['slevel'] = kLevel
      else :
        zSign = 1 if zRelToSurf else -1 # zRelToSurf => depth below surface
        zTarget = 0.0 if var in VARS2D else staZ[iSta]
        msldepth = str(int(round(zSign*zTarget*100)))
        meta['bracket'] = 'F' if zRelToSurf else 'A'
        meta['msldepth'] = msldepth
      meta['dataType'] = 'timeseries'
      z = np.mean(actualZ[iSta,:])
      x = staX[iSta]
      y = staY[iSta]
      dc = dataContainer('', ta, x,y,z, data, fieldNameList.get(var,[var]),
                          coordSys='spcs',metaData=meta)
      dcs.append(dc)
    return dcs

  def extractVerticalProfile(self,startTime, endTime, var, staX, staY, stationNames=None):
    """Extracts vertical profiles for the given time period."""
    #TODO TEST
    stacks = self.getStacks(startTime,endTime)
    time,vals,zcoords = self.getVerticalProfileForStacks(stacks, var, staX, staY, stationNames)

    #time  : array_like (nTime,)
    #vals  : array_like (nPoints, nVert, nTime)
    #zcoords : array_like (nPoints, nVert, nTime)
    # build dataContainer for each station
    dcs = []
    for iSta in range(len(staX)) :
      staName = '' if stationNames == None else stationNames[iSta]
      staStr = '{x:f} {y:f} {name:s}'.format(x=staX[iSta],y=staY[iSta],name=staName)
      # remove time steps with all bad values
      goodIxTime = np.logical_and( ~np.all( vals[iSta,:,:].mask, axis=0 ),
                                   ~np.all( zcoords[iSta,:,:].mask, axis=0 ) )
      goodIxTime = np.nonzero( goodIxTime )[0]
      v = vals[iSta,:,:][:,goodIxTime]
      z = zcoords[iSta,:,:][:,goodIxTime]
      t = time[goodIxTime]
      print '1',z.shape
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
      ta = timeArray( np.array(t), 'epoch' )
      # to (nGoodVert,nTime)
      nZ = z.shape[0]
      x = staX[iSta]*np.ones((nZ,))
      y = staY[iSta]*np.ones((nZ,))
      meta = {}
      meta['location'] = stationNames[iSta]
      meta['instrument'] = 'model'
      meta['bracket'] = 'A'
      meta['variable'] = var
      meta['dataType'] = 'profile'
      dc = dataContainer('', ta, x,y,z, data, fieldNameList.get(var,[var]),
                          coordSys='spcs',metaData=meta)
      dcs.append(dc)
    return dcs
  
  def extractTransect(self,startTime, endTime,var,staX,staY,transName):
    """Extracts a transect for the given (x,y) points and time range."""
    stacks = self.getStacks(startTime,endTime)
    staX = np.array(staX)
    staY = np.array(staY)
    time,vals,zcoords = self.getVerticalProfileForStacks(stacks, var, staX, staY)

    #time  : array_like (nTime,)
    #vals  : array_like (nPoints, nVert, nTime)
    #zcoords : array_like (nPoints, nVert, nTime)

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
    ta = timeArray( time, 'epoch' )
    meta = {}
    meta['location'] = transName
    meta['instrument'] = 'model'
    meta['bracket'] = 'A'
    meta['variable'] = var
    meta['dataType'] = 'transect'
    dc = dataContainer('', ta, X,Y,Z, data, fieldNameList.get(var,[var]),
                       coordSys='spcs', metaData=meta, acceptNaNs=True)
    return dc

  def extractTransectForBPFile(self,startTime,endTime,var,bpFile,transName):
    """Extracts a transect for the given build point ASCII file and time range."""
    from files.buildPoints import BuildPoint
    bpObj = BuildPoint()
    bpObj.readFileFromDisk(bpFile)
    return self.extractTransect(startTime,endTime,var,
                                bpObj.getX(),bpObj.getY(),transName)
    
  def extractTrack(self):
    """Extracts a (x,y,z,t) track for the given time range."""
    raise NotImplementedError('This feature has not been implemented yet.')
  
  def extractSlab(self,startTime, endTime, name, var, z=None, k=None, zRelToSurf=None):
    """Extracts a horiontal slice for the given time range."""
    stacks = self.getStacks(startTime,endTime)
    time,vals,zcoords = self.getSlabForStacks(stacks, var, z, k, zRelToSurf)
    ta = timeArray( time, 'epoch' )
    data = vals[:,None,:]
    data = data.filled(np.nan)
    connectivity = self.faceNodes
    x = self.nodeX
    y = self.nodeY
    msldepth = ''
    if k != None :
      msldepth = 'slev'+str(k)
      z = zcoords[:,0] # FIXME include time in z coords ?
    else :
      zSign = 1 if zRelativeToSurf else -1 # zRelToSurf => depth below surface
      msldepth = str(int(round(zSign*z*100)))
      zArray = z*np.ones_like(x)

    # make meshContainer
    meta = {}
    meta['dataType'] = 'slab'
    meta['location'] = name
    meta['instrument'] = 'model'
    meta['variable'] = var
    if k != None :
      meta['slevel'] = k
    else :
      meta['bracket'] = 'F' if zRelativeToSurf else 'A'
      meta['msldepth'] = msldepth
    mc = meshContainer('', ta, x,y,z, data, connectivity, fieldNameList[var], coordSys='spcs',metaData=meta)
    return mc
