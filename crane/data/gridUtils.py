"""
Methods for reading and interpolating data on unstructured triangular grids.
Only linear or constant basis functions are supported.

Tuomas Karna 2014-07-14
"""
import numpy as np
from scipy.spatial import KDTree, cKDTree
from scipy.interpolate import interp1d

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def constructCGConnectivity(x,y,connectivity,kdtree) :
  """
  Constructs connectivity table by collapsing overlapping nodes
  If the given connectivity table is discontinuous, will return equivalent
  continuous grid.
  """
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
  """
  Constructs element to neighbor element table for the given connectivity array

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
  """
  Tests whether (x,y) is in a triangle  whose vertices are nodeX,nodeY.
  Parameters
  ----------
    NodeX,nodeY -- array (nTri,3)
    x,y         -- scalar
  Returns:
    result      -- bool array (nTri,)
  """
  if u is None :
    u = barycentricCoords(nodeX,nodeY,x,y)
  return  np.logical_and( u>=0 , u<=1 ).all( axis=1 )

def interpolateTri(u,nodalValues) :
  """
  Interpolates nodal values in a location given by barycentric coordinates u
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
  values[goodIx] = np.sum(u*nv,axis=1)
  return np.ma.masked_invalid( values[:,None] )

def barycentricCoords( nodeX, nodeY, x, y ) :
  """
  Returns barycentric coordinates for (x,y) in triangle whose vertices are
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
  if k is None and z is None :
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
  vals = np.ones((nTime,))*np.nan
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
    for iTime in xrange(nTime) :
      ix = np.isfinite(ZZ[:,iTime])
      if hasattr(ZZ,'mask'):
        ix = np.logical_and(ix, ~ZZ.mask[:,iTime])
      if np.sum(ix) > 2:
        vals[iTime] = interp1d( ZZ[ix,iTime], VV[ix,iTime], kind='linear', copy=False ) ( z_target[iTime] )
  if k != None :
    # bottom: k=1 kk=0, surface: k=-1 kk=len(z)-1
    kk = k-1 if k>0 else nZ+k
    vals = VV[kk,:]
    z_target = ZZ[kk,:]
  return vals,z_target

def convertHalfLevelProfileToFullLevel(vals,zcoords,profileType='naive'):
  """
  Converts profile data extracted from a 70 file to full levels.
  
  Parameters
  ----------
  
  vals     : array_like (nProfile,nVert,nTime)
      Profile value array at half levels extracted from *.70 netcdf file
  zcoords  : array_like (nProfile,nVert,nTime)
      Profile z coordinates at whole levels computed from *.61
      (or other nodal discretization) file.
  profileType : string 'naive'|'native'
      If 'native', profile data is converted to array where each element has
      a constant value, with a discontinuous jump in between. The jump is
      actually a tiny element whose height is given by Z_DELTA.
      Otherwise no conversion is performed and the original cell center data
      and  z coordinates of cell centers is returned.
  
  Returns
  -------
  vals2    : array_like (nProfile,nVertNew,nTime)
      Profile data converted to whole levels
  zcoords2 : array_like (nProfile,nVertNew,nTime)
      Corresponding zcoordinate data
  """
  goodProfiles = ~np.all(np.all(~np.isfinite(vals),axis=1),axis=1)
  goodProfiles = np.nonzero(goodProfiles)[0]
  nProf,nVert,nTime = vals.shape
  
  if profileType == 'native':
    # native discretization, constant elements with jumps
    Z_DELTA = 1.0e-6 # small gap to mark jumps between elements
    vals2 = np.ones((nProf,2*(nVert-1),nTime))*np.nan
    zcoords2 = np.ones((nProf,2*(nVert-1),nTime))*np.nan
    for iP in goodProfiles :
      v = vals[iP,:,:] # (nVert,nTime)
      z = zcoords[iP,:,:]
      # find lowest node with real value -> bottom
      iBot = np.nonzero(~np.all(~np.isfinite(v)))[0].min()+1
      v = v[iBot:,:] # element values only
      v = np.repeat(v,2,axis=0)
      z_jump = np.repeat(z[iBot:-1,:],2,axis=0)
      z_jump[1::2] += Z_DELTA
      z = np.concatenate((z[[iBot],:],z_jump,z[[-1],:]),axis=0)
      vals2[iP,:,:] = v
      zcoords2[iP,:,:] = z
  else :  
    # naive interpolation between cell centers
    vals2 = np.ones((nProf,(nVert+1),nTime))*np.nan
    zcoords2 = np.ones((nProf,(nVert+1),nTime))*np.nan
    for iP in goodProfiles :
      v = vals[iP,:,:] # (nVert,nTime)
      z = zcoords[iP,:,:]
      # find lowest node with real value -> bottom
      iBot = np.nonzero(~np.all(~np.isfinite(v)))[0].min()+1
      vout = v[iBot:,:] # element values only
      zout = 0.5*(z[1:,:]+z[:-1,:])
      # copy bottom/surface values to ends
      vout = np.concatenate((vout[[0],:],vout,vout[[-1],:]),axis=0)
      zout = np.concatenate((z[[0],:],zout,z[[-1],:]),axis=0)
      vals2[iP,:,:] = vout
      zcoords2[iP,:,:] = zout
    
  
  vals = np.ma.masked_invalid(vals)
  zcoords = np.ma.masked_invalid(zcoords)

  return vals2,zcoords2

def convertFullLevelProfileToHalfLevel(vals) :
  """
  Converts data defined at full levels to half levels.
  
  Parameters
  ----------
  vals  : array_like (nVert,nDim1)
       Values at full levels of the vertical grid
       
  Returns
  -------
  vals2 : array_like (nVert-1,nDim1)
       Values at half levels
  """
  return 0.5*(vals[1:,...] + vals[:-1,...])

def constructEdgeNodeArray(faceNodes):
    """
    Given triangle connectivity matrix, returns element edge
    connectivity matrix. The edge/node ordering is the same as in SELFE.
    """
    import time as timeMod
    # all individual edges in the selfe order
    edge1 = faceNodes[:,[1,2]]
    edge2 = faceNodes[:,[2,0]]
    edge3 = faceNodes[:,[0,1]]
    allEdges = np.concatenate((edge1,edge2,edge3), axis=1)
    allEdges = np.reshape(allEdges, (-1,2))
    # discard duplicates ignoring order of the edge nodes
    maxNode = faceNodes.max()
    hashedEdges = np.sort(allEdges,axis=1)
    hashedEdges = hashedEdges[:,1]*2*maxNode + hashedEdges[:,0]
    # get first occurrence of any unique edge (numpy.unique can be any order)
    order = hashedEdges.argsort(kind='mergesort')
    hashedEdges = hashedEdges[order]
    diff = hashedEdges[1:] != hashedEdges[:-1]
    diff = np.concatenate([[True], diff])
    uniqueHashed = hashedEdges[diff]
    unique_ix = order[diff]
    # retain original ordering of the edges
    unique_ix = np.sort(unique_ix)
    edgeNodes = allEdges[unique_ix,:]
    return edgeNodes

#-------------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------

class meshSearch2d(object) :
  """A class for finding mesh elements that contain query points"""
  def __init__(self, node_x, node_y, faceNodes, edgeNodes=None) :
    """Initialize search tree"""
    self.node_x = node_x
    self.node_y = node_y
    self.faceNodes = faceNodes
    self.edgeNodes = edgeNodes
    self.nNodes = self.node_x.shape[0]
    self.nFaces = self.faceNodes.shape[0]
    # build 2d kd-tree
    self.tree = cKDTree( np.concatenate( (node_x[:,None],node_y[:,None]), axis=1 ) )
    a,b,c = constructCGConnectivity(self.node_x,self.node_y,self.faceNodes,self.tree)
    self.faceNodesCG = a
    self.nodeToCGNode = b
    self.cgNodeToNode = c
    # build node to element table
    self.node2elem = constructNodeToElemTable(self.faceNodesCG)
    # build elem to neighbor table
    self.elem2neigh = constructElemToNeighborTable(self.faceNodesCG,self.node2elem)
    # compute element centers
    self.center_x = np.mean(self.node_x[faceNodes],axis=1)
    self.center_y = np.mean(self.node_y[faceNodes],axis=1)
    self.tree_center = cKDTree( np.concatenate( (self.center_x[:,None],self.center_y[:,None]), axis=1 ) )

    self.edge2elem = self.elem2edge = None
    if edgeNodes != None :
      # build sturctures for edge search
      maxNodeIx = self.faceNodes.max()
      nEdge = edgeNodes.shape[0]
      nElem = self.faceNodes.shape[0]
      # edges (using SELFE convention: 1st edge opposite of 1st node etc)
      elemEdgeNodes1 = self.faceNodes[:,[1,2]]
      elemEdgeNodes2 = self.faceNodes[:,[2,0]]
      elemEdgeNodes3 = self.faceNodes[:,[0,1]]
      # hash each edge to unigue int
      w = np.array([[1,maxNodeIx]],dtype=np.uint64)
      hashEdge1 = np.sum(np.sort(elemEdgeNodes1,axis=1).astype(np.uint64)*w,axis=1)
      hashEdge2 = np.sum(np.sort(elemEdgeNodes2,axis=1).astype(np.uint64)*w,axis=1)
      hashEdge3 = np.sum(np.sort(elemEdgeNodes3,axis=1).astype(np.uint64)*w,axis=1)
      # hash edge faceNodes with same method
      hashEdgeNodes = np.sum(np.sort(self.edgeNodes,axis=1).astype(np.uint64)*w,axis=1)
      hashEdgeNodesInv = dict( zip(hashEdgeNodes,np.arange(len(hashEdgeNodes)) ) )
      # map each element to its 3 edges
      elemToHashEdge = np.vstack((hashEdge1,hashEdge2,hashEdge3)).T
      self.elem2edge = np.zeros((nElem,3),dtype=int)
      self.edge2elem = np.ones((nEdge,2),dtype=int)*-1
      for i in xrange(nElem) :
        self.elem2edge[i,0] = hashEdgeNodesInv[hashEdge1[i]]
        self.elem2edge[i,1] = hashEdgeNodesInv[hashEdge2[i]]
        self.elem2edge[i,2] = hashEdgeNodesInv[hashEdge3[i]]
        for j in xrange(3):
            ix = np.sum(self.edge2elem[self.elem2edge[i, j], :] != -1)
            self.edge2elem[self.elem2edge[i, j], ix] = i
      #elemEdge1 = np.array( [hashEdgeNodesInv[i] for i in hashEdge1] )
      #elemEdge2 = np.array( [hashEdgeNodesInv[i] for i in hashEdge2] )
      #elemEdge3 = np.array( [hashEdgeNodesInv[i] for i in hashEdge3] )
      #self.elem2edge = np.vstack((elemEdge1,elemEdge2,elemEdge3)).T
      #for iElem in xrange(nElem) :
        #self.edge2elem[self.elem2edge[iElem,:]] = iElem

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
    # construct list of all possible elements
    # center elem + its neighbors
    elems1 = self.elem2neigh[cElem,:]
    elems1 = elems1[elems1>=0] # discard -1s
    # nearest node + all elements touching it
    nearestNode = self.tree.query( np.array([[x,y]]) )[1][0]
    elems2 = self.node2elem[nearestNode]
    elems = np.hstack(([cElem],elems1,elems2))
    # test if inside an element
    e_x = self.node_x[self.faceNodes[elems,:]]
    e_y = self.node_y[self.faceNodes[elems,:]]
    goodElems = isInsideTri(x,y,nodeX=e_x,nodeY=e_y)
    if np.any(goodElems) :
      parent = elems[goodElems][0]
    else :
      parent = -1
    return parent
    
  def findParentElementArray(self, x, y) :
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
    cElems = self.tree_center.query( query )[1]
    nearestNodes = self.tree.query( query )[1]
    parents = np.ones_like(x,dtype=int)*-1
    for i in range(len(x)) :
      # list of all possible elements
      # center elem + its neighbors
      elems1 = self.elem2neigh[cElems[i],:]
      # all elements touching nearest node
      elems2 = self.node2elem[self.nodeToCGNode[nearestNodes[i]]]
      elems = np.hstack(([cElems[i]],elems1,elems2))
      # test if inside an element
      e_x = self.node_x[self.faceNodes[elems,:]]
      e_y = self.node_y[self.faceNodes[elems,:]]
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
    parentElems = self.findParentElementArray(x,y)
    goodElems = parentElems >= 0
    parentElems = parentElems[goodElems]
    xx = x[goodElems]
    yy = y[goodElems]
    
    vals = np.ones_like(x)*np.nan
    
    elemX = self.node_x[self.faceNodes[parentElems,:]]
    elemY = self.node_y[self.faceNodes[parentElems,:]]
    elemVals = nodalValues[self.faceNodes[parentElems,:]]
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

  def convertNodalValuesToEdges(self,vals,reductionOp=np.mean) :
    """Converst data defined at element vertices to edge centers.
    
    Parameters
    ----------
    vals  : array_like (nNodes,nDim1)
          Values at nodes of continuous mesh.
    
    Returns
    -------
    vals2 : array_like (nEdges,nDim1)
         Values on edge centers.
    """
    # convert to nodes of each edge (nEdge,2,nDim1)
    valsAtEdgeNodes = vals[self.edgeNodes,:]
    # take mean over each edge
    return reductionOp(valsAtEdgeNodes,axis=1)

  def convertNodalValuesToCenters(self,vals) :
    """Converst data defined at element vertices to element centers.
    
    Parameters
    ----------
    vals  : array_like (nNodes,nDim1)
          Values at nodes of continuous mesh.
    
    Returns
    -------
    vals2 : array_like (nElems,nDim1)
         Values at element centers.
    """
    # convert to nodes of each edge (nElem,3,nDim1)
    valsAtElemNodes = vals[self.faceNodes,:]
    # take mean over each edge
    return np.mean(valsAtElemNodes,axis=1)
    
  def averageElemValuesToNodes(self,vals) :
    """Converts data defined at element centers to element vertices by averaging.
    
    Parameters
    ----------
    vals  : array_like (nElems,nDim1)
        Values at element centers.
    
    Returns
    -------
    vals2 : array_like (nNodes,nDim1)
         Values at mesh nodes.
    """
    nNodes = len(self.node_x)
    nDim1 = vals.shape[1]
    vals2 = np.zeros((nNodes,nDim1),dtype=vals.dtype)
    multp = np.zeros((nNodes,nDim1),dtype=int)
    for iN in xrange(nNodes) :
      ix = self.node2elem[iN] # list
      vals2[iN,:] = np.mean( vals[ix,:], axis=0 )
    return vals2

  def convertEdgeValuesToDGNodes(self, vals) :
    """Converst data defined at element edge centers to each element vertex.
    Resulting mesh is discontinuous.
    
    Parameters
    ----------
    vals  : array_like (nEdges,nTime)
          Values at element edge centers.
    
    Returns
    -------
    vals2 : array_like (3*nElems,nTime)
         Values on discontinuous mesh. Values in axis 0 are laid out as
         [elem0_v0, elem0_v1, elem0_v2, elem1_v0, ... ]
    """
    # a) convert edges to edges for each element
    # from (nDataNodes,nTime) -> (nElem,3,nTime)
    nEdges,nTime = vals.shape
    vals = vals[self.elem2edge,:]
    # b) convert from edges to nodes in DG mesh (nTime,nElem,1) each
    n1 = (np.sum(vals[:,[1,2],:],axis=1) - vals[:,0,:]).T[...,None]
    n2 = (np.sum(vals[:,[2,0],:],axis=1) - vals[:,1,:]).T[...,None]
    n3 = (np.sum(vals[:,[0,1],:],axis=1) - vals[:,2,:]).T[...,None]

    return np.reshape(np.concatenate((n1,n2,n3),axis=2),(nTime,-1)).T

class horizontalInterpolator(object) :
  """
  A class that stores horizontal interpolation information for quickly
  interpolating over multiple fields
  """
  def __init__(self, meshSearch2d, x, y, stationNames=None):
    """Initialized the interpolation routine for the given (x,y) points"""
    self.meshSearch2d = meshSearch2d
    parentElems = self.meshSearch2d.findParentElementArray(x,y)
    self.goodElems =  parentElems >= 0
    self.parentElems = parentElems[self.goodElems]
    self.x = x
    self.y = y
    self.x_good = x[self.goodElems]
    self.y_good = y[self.goodElems]
    # nodes of all parent elements
    self.parentNodes = self.meshSearch2d.faceNodes[self.parentElems,:]
    # unique list of nodes in the parent elements, and inverse mapping
    uniqNodes, uniqNodesInv = np.unique( self.parentNodes.flatten(),
                                        return_inverse=True )
    uniqNodesInv = np.reshape(uniqNodesInv,(-1,3))
    self.uniqParentNodes = uniqNodes
    self.uniqParentNodesInv = uniqNodesInv
    # unique list of parent elements
    uniqElems, uniqElemsInv = np.unique( self.parentElems.flatten(),
                                        return_inverse=True )
    self.uniqParentElems = uniqElems
    self.uniqParentElemsInv = uniqElemsInv
    # unigue list of parent edges
    if self.meshSearch2d.elem2edge != None :
      self.parentEdges = self.meshSearch2d.elem2edge[self.parentElems,:]
      uniqEdges, uniqEdgesInv = np.unique( self.parentEdges.flatten(),
                                           return_inverse=True )
      uniqEdgesInv = np.reshape(uniqEdgesInv,(-1,3))
      self.uniqParentEdges = uniqEdges
      self.uniqParentEdgesInv = uniqEdgesInv
    else:
      self.parentEdges = self.uniqParentEdges = self.uniqParentEdgesInv = None

    self.elem_x = self.meshSearch2d.node_x[self.parentNodes]
    self.elem_y = self.meshSearch2d.node_y[self.parentNodes]
    self.u = barycentricCoords( self.elem_x, self.elem_y, self.x_good, self.y_good )
    if not np.all(self.goodElems) :
      badIx = np.nonzero(~self.goodElems)[0]
      print 'Warning: Point(s) out of grid'
      for ix in badIx :
        if stationNames :
          print '  ', stationNames[ix], self.x[ix],self.y[ix]
        else :
          print '  ', self.x[ix],self.y[ix]
    self.goodElems = np.nonzero(self.goodElems)[0]

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
    elemVals = nodalValues[self.parentNodes]
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
    
    for i in self.goodElems : # FIXME TEST
      elemVals = nodalValues[self.parentNodes[i,:],:,:]
      # from shape (3,nDim2,nDim3) to (nDim3,nDim2,3)
      elemVals = np.swapaxes(elemVals,0,2)
      # to (nDim3*nDim2,3)
      elemVals = np.reshape(elemVals,(-1,3))
      vi = interpolateTri(self.u[i,:], elemVals )
      # from  to shape (nDim2, nDim3)
      vi = np.reshape( vi, (nDim3,nDim2) ).T
      vals[i,:,:] = vi

    return vals
  
  def getParentNodes(self, discretizationType='node') :
    """
    Returns a list of unique node indices that form the parent elements.
    Can be used to limit reading to only minimum number of nodal values.

    Returns
    -------
    values      : array_like (nUniqParentNodes,)
                  Unique parent nodes

    """
    if discretizationType == 'node' :
      #print 'parent nodes (node)',self.uniqParentNodes
      return self.uniqParentNodes
    elif discretizationType == 'face' :
      #print 'parent nodes (face)',self.uniqParentElems
      return self.uniqParentElems
    elif discretizationType == 'edge' :
      #print 'parent nodes (edge)',self.uniqParentEdges
      if self.parentEdges is None :
        raise Exception('This object doesn\'t have edge information.')
      return self.uniqParentEdges

  def mapDataToTriVertices(self, data, discretizationType='node') :
    """
    Maps data in native discretization to triangle vertices.
    
    Parameters
    ----------
    data : array_like (nUniqParentNodes,)
           Nodal data defined at unique parent nodes
    
    Returns
    -------
    vals : array_like (nPoints,3)
           Nodal data converted to triangle vertices for each (x,y) point
    """
    if discretizationType == 'node' :
      # map unique nodal values to each element
      return data[self.uniqParentNodesInv]
    elif discretizationType == 'face' :
      # copy center value to all nodes
      return np.tile(data[self.uniqParentElemsInv],(3,1))
    elif discretizationType == 'edge' : 
      if self.parentEdges is None :
        raise Exception('This object doesn\'t have edge information.')
      # map unique edges to edges of each element (nPoints,3)
      dataEdges = data[self.uniqParentNodesInv]
      # map edge values to values at vertices (sum of nearest edges - opposite)
      vals1 = np.sum(dataEdges[:,[1,2]],axis=1) - dataEdges[:,0]
      vals2 = np.sum(dataEdges[:,[2,0]],axis=1) - dataEdges[:,1]
      vals3 = np.sum(dataEdges[:,[0,1]],axis=1) - dataEdges[:,2]
      return np.vstack((vals1,vals2,vals3)).T

  def mapDataArrToTriVertices(self, data, discretizationType='node') :
    """
    Maps data in native discretization to triangle vertices.
    
    Parameters
    ----------
    data : array_like (nUniqParentNodes,dim2,dim3)
           Nodal data defined at unique parent nodes (of native discretization)
    
    Returns
    -------
    
    vals : array_like (nPoints,3,dim2,dim3)
           Nodal data converted to triangle vertices for each (x,y) point
    
    """
    if discretizationType == 'node' :
      return data[self.uniqParentNodesInv,:,:]
    elif discretizationType == 'face' :
      # copy center value to all nodes
      return np.tile(data[self.uniqParentElemsInv,:,:][:,None,:,:],(1,3,1,1))
    elif discretizationType == 'edge' :
      if self.parentEdges is None :
        raise Exception('This object doesn\'t have edge information.')
      # map unique edges to edges of each element (nPoints,3,dim2,dim3)
      dataEdges = data[self.uniqParentEdgesInv,:,:]
      # map edge values to values at vertices (sum of nearest edges - opposite)
      vals1 = np.sum(dataEdges[:,[1,2],:,:],axis=1) - dataEdges[:,0,:,:]
      vals2 = np.sum(dataEdges[:,[2,0],:,:],axis=1) - dataEdges[:,1,:,:]
      vals3 = np.sum(dataEdges[:,[0,1],:,:],axis=1) - dataEdges[:,2,:,:]
      return np.concatenate((vals1[:,None,:,:],vals2[:,None,:,:],
                             vals3[:,None,:,:]),axis=1)

  def evaluateFromParentNodes(self, nodalValues, discretizationType='node') :
    """
    Interpolates the field defined by the nodal values of only the necessary
    nodes. The necessary nodes are stored in self.uniqParentNodes member.
    
    Parameters
    ----------
    nodalValues : array_like (nParentNodes,)
                  Nodal values of all necessary nodes
                  
    Returns
    -------
    values      : array_like (nPoints,)
                  Interpolated values
    """
    elemVals = self.mapDataToTriVertices( nodalValues, discretizationType )
    vals = np.ones_like(self.x)*np.nan
    vals[self.goodElems] = interpolateTri( self.u, elemVals )
    return vals

  def evaluateArrayFromParentNodes(self, nodalValues, discretizationType='node') :
    """
    Interpolates the field defined by the nodal values of only the necessary
    nodes. The necessary nodes are stored in self.uniqParentNodes member.

    
    Parameters
    ----------
    nodalValues : array_like (nParentNodes,nDim2,nDim3)
                  Nodal values of a field defined on the mesh for multiple
                  depths and time instances

    Returns
    -------
    values      : array_like (nPoints,nDim2,nDim3)
                  Interpolated values
    """
    nPoints = len(self.x)
    nNodes, nDim2, nDim3 = nodalValues.shape
    vals = np.ones((nPoints,nDim2,nDim3))*np.nan

    elemValsArr = self.mapDataArrToTriVertices(nodalValues, discretizationType)

    for i,iE in enumerate(self.goodElems) :
      elemVals = elemValsArr[i,:,:,:]
      # from shape (3,nDim2,nDim3) to (nDim3,nDim2,3)
      elemVals = np.swapaxes(elemVals,0,2)
      # to (nDim3*nDim2,3)
      elemVals = np.reshape(elemVals,(-1,3))
      vi = interpolateTri(self.u[i,:], elemVals )
      # from  to shape (nDim2, nDim3)
      vi = np.reshape( vi, (nDim3,nDim2) ).T
      vals[iE,:,:] = vi

    return vals

class verticalInterpolator(object) :
  """A class that interpolates data in vertical."""
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
    if k is None and z is None :
      raise Exception('Either k or z must be defined')
    self.doZInterpolation = k is None
    self.z = z
    self.k = k
    if self.doZInterpolation and zRelToSurf is None :
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
    if hasattr(V,'mask'):
      goodProfiles = np.logical_and(goodProfiles,
                                    ~np.all(np.all(V.mask,axis=2),axis=1))
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

