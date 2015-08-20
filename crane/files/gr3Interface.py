#!/usr/bin/env python
"""
Collection of methods for reading/writing gr3 unstructured grid files.

Tuomas Karna 2012-12-04
"""
import string
import os
import sys
import numpy as np

from crane.data import timeArray
from crane.data import meshContainer

def readNodeBlock(infile, nPoints) :
  # iN x y data
  N = 4
  nodes = np.fromfile(infile,dtype=float,count=N*nPoints,sep=' ').reshape((nPoints,N))
  return nodes

def readElemBlock(infile, nElems) :
  # iT,elemType,n1,n2,n3
  N = 5
  try :
    elems = np.fromfile(infile,dtype=int,count=N*nElems,sep=' ').reshape((nElems,N))
    return elems
  except :
    raise Exception( 'Could not read elements' )

def readNodalValues(inputFilename) :
  """Reads gr3 file, returns x,y,val arrays"""
  print 'Reading '+inputFilename+' ...'
  infile = open(inputFilename,'r')
  description = infile.readline().strip() # remove leading/trailing whitespace
  tmpStr = infile.readline()
  nTriangles, nNodes = ( int(s) for s in tmpStr.split() )
  print '  nTriangles={0:d} nNodes={1:d}'.format(nTriangles,nNodes)

  ### nodes
  nodeArray = readNodeBlock(infile,nNodes)
  nodenum   = np.array( nodeArray[:,0].flatten(), dtype=int )
  nodexyz = np.zeros((nNodes,3))
  nodexyz[:,:2] = nodeArray[:,1:3]
  nodalValues = nodeArray[:,3]

  return nodexyz[:,0],nodexyz[:,1],nodalValues

def readGR3FileToMC(inputFilename, fieldName='depth') :
  x,y,vals,conn,boundaries,description = readGR3File(inputFilename)
  data = np.reshape(vals, (-1,1,1) )
  z = np.zeros_like(x)
  ta = timeArray.timeArray( np.array([0]), 'epoch' )
  mc = meshContainer.meshContainer( description, ta,x,y,z,data,conn,[fieldName],coordSys='spcs')
  # append boundary information
  for bnd in boundaries :
    print bnd.type, bnd.tag, len(bnd.nodes)
    mc.addBoundary( bnd )
  
  return mc

def readGR3File(inputFilename) :
  """Reads gr3 file, returns data in separate arrays."""
  print 'Reading '+inputFilename+' ...'
  infile = open(inputFilename,'r')
  description = infile.readline().strip() # remove leading/trailing whitespace
  tmpStr = infile.readline()
  nTriangles, nNodes = ( int(s) for s in tmpStr.split() )
  print '  nTriangles={0:d} nNodes={1:d}'.format(nTriangles,nNodes)

  ### nodes
  nodeArray = readNodeBlock(infile,nNodes)
  nodenum   = np.array( nodeArray[:,0].flatten(), dtype=int )
  nodexyz = np.zeros((nNodes,3))
  nodexyz[:,:2] = nodeArray[:,1:3]
  nodalValues = nodeArray[:,3]

  print '  Nodal values min={0:g} max={1:g}'.format(min(nodalValues),max(nodalValues) )

  ### triangular elements
  triArray = readElemBlock(infile, nTriangles)

  trinum   = triArray[:,0].flatten()
  tritype  = triArray[0,1]
  trinodes = triArray[:,-3:]-1 # three last columns, 0-based indexing
  #triangles = meshElements(trinodes,trinum,tritype)

  x = nodexyz[:,0]
  y = nodexyz[:,1]

  tmpStr = infile.readline()
  boundaries = []
  if len(tmpStr) > 0 :
    ### boundary information, if not end of file
    nOpenBndSegments = int(tmpStr.split()[0])
    nOpenBndNodesTot = int(infile.readline().split()[0])
    print '  nOpenBndSegments={0:d} nOpenBndNodesTot={1:d}'.format(nOpenBndSegments, nOpenBndNodesTot)
    for iBnd in range(nOpenBndSegments) :
      bndHeader = infile.readline().split()
      nBndNodes = int(bndHeader[0])
      tag = bndHeader[-1]
      if tag.isdigit() :
        tag = 'open'+tag
      print '   open bnd {0:d} {1:s}: {2:d} nodes'.format(iBnd+1,tag,nBndNodes)
      tmpList = []
      for iN in range(nBndNodes) :
        tmpList.append( int(infile.readline()) )
      nodes = np.array(tmpList,dtype=int)-1
      boundaries.append( meshContainer.meshBoundary( 'open', tag, nodes ) )
    nLandBndSegments = int(infile.readline().split()[0])
    nLandBndNodesTot = int(infile.readline().split()[0])
    landBndTags = range(nOpenBndSegments+1,nOpenBndSegments+nLandBndSegments+1)
    print '  nLandBndSegments={0:d} nLandBndNodesTot={1:d}'.format(nLandBndSegments, nLandBndNodesTot)
    for iBnd in range(nLandBndSegments) :
      bndHeader = infile.readline().split()
      nBndNodes = int(bndHeader[0])
      try :
        landType = int(bndHeader[1])
      except :
        print """Land boundary type missing in gr3 file. Add 0/1 (land/island) after number of nodes in each land boudary, e.g.
        1002 = Total number of closed boundary nodes
        501 0 = Number of nodes in closed boundary 1"""
        raise Exception('Could not parse land boundary type (0/1 - land/island)\n')
      landType = 'island' if landType == 1 else 'land'
      tag = landType+bndHeader[-1]
      print '   land bnd {0:d} {1:s}: {2:d} nodes'.format(iBnd+1,tag,nBndNodes)
      tmpList = []
      for iN in range(nBndNodes) :
        tmpList.append( int(infile.readline()) )
      #tmpList = fromfile(infile,dtype=int,count=nBndNodes,sep=' ')
      nodes = np.array(tmpList,dtype=int)-1
      boundaries.append( meshContainer.meshBoundary( landType, tag, nodes ) )

  infile.close()

  # for better interpolation, round coordinates to 1e-4
  nDig = 4
  x = np.round(x,nDig)
  y = np.round(y,nDig)

  return x,y,nodalValues,trinodes,boundaries,description

def findElementsInRegion(nodalValues,triangleList,targetTag) :
  #goodElems = []
  #for ie in range(len(triangleList)) :
    #regionTag = max( [ nodalValues[i-1] for i in triangleList[ie] ] )
    #if regionTag == targetTag :
      #goodElems.append( ie )
  #return goodElems
  ix = np.nonzero( np.sum( nodalValues[triangleList] == targetTag, axis=1 ).flatten() == 3 )
  if len(ix) == 0 :
    print targetTag
    print nodalValues.min(),nodalValues.max()
    raise Exception('No region found')
  return ix[0]

def writeMCToGR3File( filename, mc ) :
  """Writes MeshContainer to a gr3 file"""
  nodes = np.vstack( (mc.x,mc.y) ).T
  nodalValues = mc.data[:,0,0].squeeze()[:,None]
  connectivity = mc.connectivity
  openBndNodes = []
  landBndNodes = []
  writeGR3File( filename, '' ,nodes, nodalValues, connectivity, mc.boundaries)

def writeGR3File( filename, description ,nodes, nodalValues, connectivity, boundaries=None, writeTriangles=True) :
  print 'Writing '+filename+' ...'
  triangles = connectivity+1
  outfile = open(filename,'w')
  if description == None or len(description) == 0 :
    description = os.path.split(filename)[1] # take filename without path
  outfile.write( description+'\n' )
  ne = triangles.shape[0]
  nPoints = nodes.shape[0]
  outfile.write( '{ne:d} {nPoints:d}\n'.format(ne=ne,nPoints=nPoints) )
  if nodalValues == None :
     nodalValues = nodes[:,2]
  triType = 3
  for i in range(nPoints) :
    line = repr(i+1) + ' ' + string.join( nodes[i,:2].astype('|S30'),' ') + ' ' + repr(nodalValues[i,0])
    outfile.write( line+'\n' )
  if writeTriangles :
    for i in range(ne) :
      line = repr(i+1) +' '+ repr(triType) +' '+ string.join( triangles[i,:].astype('|S30'),' ')
      outfile.write( line+'\n' )
  # n open bnds
  if boundaries :
    openBnds = []
    landBnds = []
    for bnd in boundaries :
      if bnd.type == 'open' :
        openBnds.append( bnd )
      else :
        landBnds.append( bnd )
    # sort land bnds before island bnds
    newLandBnds = []
    for bnd in landBnds :
      if bnd.type == 'land' :
        newLandBnds.append( bnd )
    for bnd in landBnds :
      if bnd.type != 'land' :
        newLandBnds.append( bnd )
    landBnds = newLandBnds
    #outfile.write( '{0:d} = Number of open boundaries\n'.format(len(openBnds)) )
    #totBndNodes = 0
    #for bnd in openBnds :
      #totBndNodes += len(bnd.nodes)
    ## tot n open bnd nodes
    #outfile.write( '{0:d} = Total number of open boundary nodes\n'.format(totBndNodes) )
    #for i,bnd in enumerate(openBnds) :
      #tagStr = '' if bnd.tag.isdigit() else bnd.tag
      ## n nodes in bnd 1
      #outfile.write( '{0:d} = Number of nodes for open boundary {1:d} {2:s}\n'.format(len(bnd.nodes),i+1,tagStr) )
      ## list of nodes
      #outfile.write( string.join( bnd.nodes.astype('|S30'),'\n')+'\n' )
    ## n land bnds
    def writeBnds(bndType,bndList) :
      #if len(bndList)==0 : return
      #bndType = 'open' if bndList[0].type == 'open' else 'land'
      outfile.write( '{0:d} = Number of {1:s} boundaries\n'.format(len(bndList),bndType) )
      totBndNodes = 0
      for bnd in bndList :
        totBndNodes += len(bnd.nodes)
      # tot n land bnd nodes
      outfile.write( '{0:d} = Total number of {1:s} boundary nodes\n'.format(totBndNodes,bndType) )
      for i,bnd in enumerate(bndList) :
        tagStr = '' if bnd.tag.isdigit() else bnd.tag
        islandTag = ''
        if bnd.type == 'island' : islandTag = ' 1'
        if bnd.type == 'land' : islandTag = ' 0'
        # n nodes in bnd 1
        outfile.write( '{n:d} {iTag:s} = Number of nodes for {bTyp:s} boundary {iBnd:d} {tag:s}\n'.format(n=len(bnd.nodes),iTag=islandTag,bTyp=bnd.type,iBnd=i+1,tag=tagStr) )
        # list of nodes
        nodes = bnd.nodes + 1
        outfile.write( string.join( nodes.astype('|S30'),'\n')+'\n' )
    writeBnds( 'open',openBnds )
    writeBnds( 'land',landBnds )
