#!/usr/bin/python
"""
Generic implementation for a mesh container.
Can represent 1D/2D/3D unstructured grid data.
Contains data, necessary information to interpret it, and some basic
data manipulation methods.

Tuomas Karna 2012-11-26
"""
import numpy as np
from crane.data import dataContainer
from crane.data import netcdfIO


class meshContainer(dataContainer.dataContainer) :
  """
  Generic mesh container.
  Contains an array of node coordinates,
  and an array of elements (connectivity array).
  Cannot represent mixed meshes (very well).
  """
  # TODO add a tag for element types, pol. order, discretization type (continuous, discontinuos, non-conf. etc) ?
  # TODO element type (etc) could be an array for mixed mesh
  # TODO test
  def __init__(self, description, time, x,y,z, data, connectivity, fieldNames, coordSys='',metaData={}, dtype=np.float64) :
    """Creates new meshData object

    description -- string to identify the data array
    time -- timeArray object
    x,y,z -- (nPoints,3) array of coordinates
    connectivity -- (Nelems, Nvertices) element connectivity
    data -- (nPoints,nFields,nTime) array of p fields
    fieldNames -- list of p strings
    coordSys   -- string representing the used coordinate system
    """
    nElem = connectivity.shape[0]
    self.boundaries = []
    checkDataDim = True
    self.dataByElement = False
    if data.shape[0] == nElem :
      self.dataByElement = True
      checkDataDim = False
    if isinstance(time,(int,float,type(None))) :
      time = timeArray.timeArray(np.array([0]),'epoch')
    super(meshContainer, self).__init__(description, time, x,y,z, data,
                           fieldNames, coordSys, metaData, acceptNaNs=True,
                           checkDataXDim=checkDataDim,dtype=dtype)
    self.connectivity = connectivity.astype(np.int32)

  @classmethod
  def fromDataContainer( cls, dc, connectivity ) :
    """Creates meshContainer from dataContainer, enriched by the given connectivity table."""
    return cls( dc.description, dc.time, dc.x, dc.y, dc.z, dc.data,
                connectivity, dc.fieldNames, dc.coordSys, dc.metaData)

  def __eq__(self, other) :
    """True if all data is equal to other."""
    if not super(meshContainer,  self).__eq__(other) :
      return False
    if not np.array_equal( self.connectivity, other.connectivity ) :
      return False
    return True

  def __ne__(self, other) :
    """True if some data is equal to other."""
    return not self.__eq__(other)

  def addBoundary( self, bnd ) :
    """Adds meshBoundary object in this container."""
    if not isinstance( bnd, meshBoundary ) :
      raise Exception( 'boundary must be meshBoundary object')
    n = bnd.nodes
    if n.max() > self.x.shape[0]-1 :
      raise Exception('Boundary node index exceed number of nodes in the mesh ({0:d} > {1:d})'.format(n.max(),self.x.shape[0]-1) )
    self.boundaries.append( bnd )

  def extractFields( self, *fields ) :
    """
    Returns a meshContainer containing only the requested field.

    Parameters
    ----------
    fields : str or int
        fieldName to extract. If int, the index of the field to extract.
    """
    indices = []
    names = []
    for f in fields :
      if isinstance( f, str ) or isinstance( f, unicode ) :
        # deduce index
        i = self.fieldNames.index( f )
      else :
        i = f
      indices.append( i )
      names.append( self.fieldNames[i] )
    data = self.data[:,indices,:]

    return meshContainer( self.description, self.time, self.x, self.y, self.z, data, self.connectivity, names, self.coordSys, self.metaData )

  def duplicateWithFields( self, fields, fieldNames ) :
    """Creates a new meshContainer by replacing the field data with given array(s)."""
    nTime = len(self.time)
    nPoints = self.x.shape[0]
    for i in range(len(fields)) :
      if len(fields[i].shape) < 2 :
        fields[i] = fields[i][:,None].copy()
      fields[i] = fields[i][:,None,:].copy()
    data = np.hstack( fields )
    return meshContainer( self.description, self.time, self.x, self.y, self.z, data,
                self.connectivity, fieldNames, self.coordSys, self.metaData)

  def convertToNodalMC( self ) :
    """If data is at elements, returns a new object with average value at each node."""
    if self.dataByElement :
      # compute mean aroung each element TODO weighted by elem size
      nodal_data = np.zeros((nNodes,nFields,nTime))
      nodal_data = np.zeros((nNodes,nFields,nTime))
      conn = self.connectivity
      for iElem in range(self.connectivity.shape[0]) :
        n = conn[iElem,:]
        nodal_data[ n,:,: ] += self.data[iElem,:,:]
        nodal_multiplicity[ n ] += 1
      nodal_data = nodal_data/nodal_multiplicity

      mc = meshContainer( self.description, self.time, self.x, self.y, self.z,
                        nodal_data, self.connectivity, self.fieldNames,
                        self.coordSys, self.metaData)
      mc.boundaries = list(self.boundaries)
      return mc
    else :
      return self.copy()

  @classmethod
  def loadFromNetCDF( cls, filename, startTime=None, endTime=None, includeEnd=False ) :
    """Creates a new dataContainer from netCDF file.
    """
    nc = netcdfIO.netcdfIO(filename)
    mc = nc.readToDataContainer( startTime, endTime, includeEnd=includeEnd )
    return mc

  def interpolateInTime(self, newTime, acceptNaNs=False) :
    """
    Interpolates data to newTime. newTime is a timeArray object.
    Returns a new dataContainer object.
    Note: x,y,z are referenced rather than copied!
    """
    dc = super(meshContainer, self).interpolateInTime(newTime, acceptNaNs)
    return meshContainer.fromDataContainer(dc, self.connectivity)

  def cropGrid( self, boundingBox ) :
    """Returns a new meshContainer where all elements outside the boundingBox
    have been discarded. boundingBox = [xmin,xmax,ymin,ymax]."""
    
    ## find all elements whose center is in the box
    #xCntr = np.mean(self.x[self.connectivity],axis=1)
    #yCntr = np.mean(self.y[self.connectivity],axis=1)
    #goodElems = np.logical_and( xCntr >= boundingBox[0], xCntr <= boundingBox[1] )
    #goodElems = np.logical_and( goodElems ,
                #np.logical_and( yCntr >= boundingBox[2], yCntr <= boundingBox[3] ) )

    # find all elements with at least one node in the box
    goodX = np.logical_and( self.x >= boundingBox[0], self.x <= boundingBox[1] )
    goodY = np.logical_and( self.y >= boundingBox[2], self.y <= boundingBox[3] )
    goodElems = np.logical_and( goodX[self.connectivity].sum(axis=1),
                                goodY[self.connectivity].sum(axis=1) )

    conn = self.connectivity[goodElems,:]
    # TODO generalize for moving mesh
    nodes = np.sort(np.unique(conn))
    x = self.x[nodes]
    y = self.y[nodes]
    z = self.z[nodes]
    if self.dataByElement :
      data = self.data[goodElems,:,:]
    else :
      data = self.data[nodes,:,:]
    # inverse mapping, old node to new node ix
    newNodeIndex = -1*np.ones( (len(self.x)), dtype=int )
    newNodeIndex[nodes] = range(len(nodes))
    # express connectivity with new node ix
    conn = newNodeIndex[conn]

    newMC = meshContainer(self.description, self.time, x,y,z, data, conn,
                           self.fieldNames, self.coordSys, self.metaData)
    return newMC

  def copy(self) :
    """Deep copy, all numpy arrays are copied instead of referenced."""
    return meshContainer(self.description, self.time.copy(),
                         self.x.copy(), self.y.copy(), self.z.copy(),
                         self.data.copy(), self.connectivity.copy(), list(self.fieldNames), self.coordSys, dict(self.metaData))

  def timeWindow(self, startDate, endDate,includeEnd=False) :
    """Returns a meshContainer, restricted to the given interval.
    startDate and endDate are datetime objects. Data in the returned array is refenreced, not copied.
    To obtain a copy use timeWindow(startDate, endDate).copy()
    """
    dc = super(meshContainer, self).timeWindow(startDate, endDate, includeEnd)
    return meshContainer.fromDataContainer(dc, self.connectivity )

  def subsample( self, timeStamps=None, skipFactor=None, targetDt=None, currentDt=None, gapFactor=5) :
    """Subsamples the data with the given skipFactor or targetDt.
    If timeStamps is not given, appropriate time indices will be estimated.
    Returns a new dataContainer.

    Args:
    timeStamps -- (ndarray) array of time indices to include
    skipFactor -- (int) subsampling factor, accept every skipFactor data point.
    targetDt -- (float) alternatively estimate skipFactor based on targetDt (sec).
    currentDt -- (float) data sampling period. If None, taken as a mean step between data points. Passed to detectGaps.
    gapFactor -- (float) A factor to determine a minimum gap: gapFactor*dt. Passed to detectGaps.

    Returns:
    newDC      -- (array) subsampled version of this dataContainer
    """
    dc = super(meshContainer, self).subsample(timeStamps, skipFactor, targetDt, currentDt, gapFactor)
    return meshContainer.fromDataContainer(dc, self.connectivity)

  def fixTriangleArea( self ) :
    """Permutate triangle vertices so that element area is positive (in SELFE)"""
    areas = self.computeAreas()
    for i in range(self.connectivity.shape[0]) :
      tri = self.connectivity[i,:]
      area = areas[i]
      if area < 0 :
        # swap first and last node
        tmp = self.connectivity[i,0]
        self.connectivity[i,0] = self.connectivity[i,2]
        self.connectivity[i,2] = tmp

  def computeAreas(self) :
    """Compute the signed area of each element"""
    return computeAreas(self.connectivity, self.x, self.y)

  def computeEquivalentRadius(self) :
    """Compute the equivalent radius of each element"""
    return computeEquivalentRadius(self.connectivity, self.x, self.y)

  def computeInnerRadius(self) :
    """Compute the radius of incircle of each element"""
    return computeInnerRadius( self.connectivity, self.x, self.y)

  def computeTriangleQuality(self) :
    """Computes triangle quality measure for each element
    
    Quality measure is defined as
    Q = 4*A*sqrt(3)/(a**2+b**2+c**2)
    where A is the area and a,b,c are the side lengths.
    Typically q>=0.6 is acceptable. q=1 for equilateral triangle.
    """
    return computeTriangleQuality(self.connectivity, self.x, self.y)

  def computeMaxAngle(self) :
    """Computes the maximum angle in each element (in degrees)"""
    return computeMaxAngle(self.connectivity, self.x, self.y)

  def computeMaxEdgeLen(self) :
    """Computes the maximum edge length in each element"""
    return computeMaxEdgeLen(self.connectivity, self.x, self.y)

def computeAreas(connectivity,x,y) :
  """Compute the signed area of each element"""
  x1 = x[connectivity[:,0]]
  x2 = x[connectivity[:,1]]
  x3 = x[connectivity[:,2]]
  y1 = y[connectivity[:,0]]
  y2 = y[connectivity[:,1]]
  y3 = y[connectivity[:,2]]
  return ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3)) * 0.5

def computeEquivalentRadius(connectivity,x,y) :
  """Compute the equivalent radius of each element"""
  return 2.0 * np.sqrt(computeAreas(connectivity,x,y)) / np.pi

def computeInnerRadius(connectivity,x,y) :
  """Compute the radius of incircle of each element"""
  x1 = x[connectivity[:,0]]
  x2 = x[connectivity[:,1]]
  x3 = x[connectivity[:,2]]
  y1 = y[connectivity[:,0]]
  y2 = y[connectivity[:,1]]
  y3 = y[connectivity[:,2]]
  # side lengths
  a = np.sqrt( (x1-x2)**2 + (y1-y2)**2 )
  b = np.sqrt( (x2-x3)**2 + (y2-y3)**2 )
  c = np.sqrt( (x1-x3)**2 + (y1-y3)**2 )
  s = 0.5*(a+b+c) # semiperimeter
  return np.sqrt( ( (s-a)*(s-b)*(s-c) ) / s )

def computeTriangleQuality(connectivity,x,y):
  """Computes triangle quality measure for each element
  
  Quality measure is defined as
  Q = 4*A*sqrt(3)/(a**2+b**2+c**2)
  where A is the area and a,b,c are the side lengths.
  Typically q>=0.6 is acceptable. q=1 for equilateral triangle.
  """
  x1 = x[connectivity[:,0]]
  x2 = x[connectivity[:,1]]
  x3 = x[connectivity[:,2]]
  y1 = y[connectivity[:,0]]
  y2 = y[connectivity[:,1]]
  y3 = y[connectivity[:,2]]
  A = ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3)) * 0.5
  # side lengths
  a = np.sqrt( (x1-x2)**2 + (y1-y2)**2 )
  b = np.sqrt( (x2-x3)**2 + (y2-y3)**2 )
  c = np.sqrt( (x1-x3)**2 + (y1-y3)**2 )
  return 4*A*np.sqrt(3.0)/ (a**2+b**2+c**2)

def computeMaxAngle(connectivity,x,y):
  """Computes the maximum angle in each element (in degrees)"""
  x1 = x[connectivity[:,0]]
  x2 = x[connectivity[:,1]]
  x3 = x[connectivity[:,2]]
  y1 = y[connectivity[:,0]]
  y2 = y[connectivity[:,1]]
  y3 = y[connectivity[:,2]]
  # side lengths
  a = np.sqrt( (x1-x2)**2 + (y1-y2)**2 )
  b = np.sqrt( (x2-x3)**2 + (y2-y3)**2 )
  c = np.sqrt( (x1-x3)**2 + (y1-y3)**2 )
  # max angle is opposite to longest side
  L = np.vstack((a,b,c)).T
  L = np.sort(L,axis=1)
  theta = np.arccos( (L[:,0]**2 + L[:,1]**2 - L[:,2]**2)/(2*L[:,0]*L[:,1]) )
  return theta/2.0/np.pi*360.0

def computeMaxEdgeLen(connectivity,x,y):
  """Computes the maximum angle in each element (in degrees)"""
  x1 = x[connectivity[:,0]]
  x2 = x[connectivity[:,1]]
  x3 = x[connectivity[:,2]]
  y1 = y[connectivity[:,0]]
  y2 = y[connectivity[:,1]]
  y3 = y[connectivity[:,2]]
  # side lengths
  a = np.sqrt( (x1-x2)**2 + (y1-y2)**2 )
  b = np.sqrt( (x2-x3)**2 + (y2-y3)**2 )
  c = np.sqrt( (x1-x3)**2 + (y1-y3)**2 )
  # max angle is opposite to longest side
  L = np.vstack((a,b,c)).T
  L = np.max(L,axis=1)
  return L


class meshBoundary(object) :
  """Object for storing mesh boundary information"""
  def __init__( self, type, tag, nodes ) :
    """Create a new boundary.
    Arguments:
    type -- ('open' or 'land') type of the boundary
    tag -- (str) String that identifies the boundary (e.g. 'pacific', 'beaver')
    nodes -- (singleton ndarray) ordered array of node indices that belong to the boundary.
    """
    if type not in ['open','land','island'] :
      raise Exception('boundary type must be either land or open')
    self.type = type
    self.tag = tag
    self.nodes = nodes

  @classmethod
  def openBoundary( cls, tag, nodes ) :
    """Shorthand for creating open boundary"""
    return cls( 'open', tag, nodes )
  
  @classmethod
  def landBoundary( cls, tag, nodes ) :
    """Shorthand for creating land boundary"""
    return cls( 'land', tag, nodes )
  
  def __str__(self) :
    outStr = '{type:s} {tag:s}\n'.format(type=self.type, tag=self.tag)
    outStr += str(self.nodes)
    return outStr
  
  def copy(self) :
    """Returns a deep copy of this object"""
    return meshBoundary( self.type, self.tag, self.nodes.copy() )
