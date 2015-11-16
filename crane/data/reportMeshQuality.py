"""
Computes element quality measures and reports worst elements

Tuomas Karna 2014-06-18
"""
import numpy as np
import os

from crane.data import meshContainer
from crane.files import gmshInterface
from crane.files import gr3Interface

def readAnyMeshFile( inFile, dataFile=None ) :
  """Reads mesh data in a meshContainer. Supported formats are GR3 (SELFE), MSH (GMSH) and meshContainer netCDF."""
  fname,ext = os.path.splitext(inFile)
  if ext == '.gr3' :
    mc = gr3Interface.readGR3FileToMC(inFile)
  elif ext == '.msh' :
    gmsh = gmshInterface.gmshMesh.fromFile( inFile )
    if dataFile :
      mc = gmsh.getMeshContainer( fromFile=dataFile )
    else :
      mc = gmsh.getMeshContainer( fromFile=inFile )
  elif ext == '.nc' :
    mc = meshContainer.meshContainer.loadFromNetCDF(inFile)
  else :
    raise Exception('Unknown file extension: '+ext+'\n'+inFile)
  return mc

def reportMeshQuality(meshFile,maxAngle=None,minQuality=None, oneBasedIndex=False) :
  mc = readAnyMeshFile( meshFile )
  
  elemOffset = 1 if oneBasedIndex else 0
  
  angles = mc.computeMaxAngle()
  quality = mc.computeTriangleQuality()
  
  if maxAngle != None :
    ix = np.argsort(angles)[::-1]
    print_ix = np.nonzero(angles[ix] >= maxAngle)[0]
    print '---'
    print 'iElem    angle  quality'
    for i in print_ix :
      print '{i:8d} {a:6.2f} {q:6.4f}'.format(i=ix[i]+elemOffset, a=angles[ix][i], q=quality[ix][i])

  if minQuality != None :
    ix = np.argsort(quality)
    print_ix = np.nonzero(quality[ix] <= minQuality)[0]
    print '---'
    print 'iElem    angle  quality'
    for i in print_ix :
      print '{i:8d} {a:6.2f} {q:6.4f}'.format(i=ix[i]+elemOffset, a=angles[ix][i], q=quality[ix][i])

#-------------------------------------------------------------------------------
# Command line interface
#-------------------------------------------------------------------------------
def parseCommandLine() :

  from optparse import OptionParser
  usage = 'Usage: %prog [options] meshFile'
  parser = OptionParser(usage=usage)
  parser.add_option('-a', '--maxAngle', action='store', type='float',
                      dest='maxAngle', help='Define the maximum angle in triangle to classify element as acceptable (e.g. 160)')
  parser.add_option('-q', '--minQuality', action='store', type='float',
                      dest='minQuality', help='Define minimum element quality to classify element as acceptable (e.g. 0.05)')
  parser.add_option('-O', '--one-based-index', action='store_true',default=False,
                      dest='oneBasedIndex', help='Report element indices starting from 1 (e.g. in gr3 file) instead of 0 (Default: %default)')

  (options, args) = parser.parse_args()
  
  if len(args) == 0 :
    parser.print_help()
    parser.error('mesh file undefined')

  meshFile = args[0]

  if options.maxAngle is None and options.minQuality is None :
    parser.print_help()
    parser.error('either maxAngle or minQuality must be given')

  print 'Parsed options:'
  print ' - mesh file',meshFile
  if options.maxAngle :
    print ' - max angle', options.maxAngle
  if options.minQuality :
    print ' - min quality', options.minQuality
  print ' - use 1 as first element number:', options.oneBasedIndex

  reportMeshQuality( meshFile, options.maxAngle, options.minQuality, options.oneBasedIndex )
  
if __name__=='__main__' :
  parseCommandLine()
