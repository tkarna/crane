#!/usr/bin/python
"""
Extract transect data using the efficient SELFE extract_mod python module.

Examples:

# extract vdff, -d data dir, -t defines transect bp file, -n transect name string, -o output dir, -s -e time range
python extractTransect.py -d ~pturner/db29/run29/outputs/ -v vdff -t nchannel_fine.bp -n nchannel -o run29/transect -s 2012-5-1 -e 2012-5-17

Tuomas Karna 2012-11-16
"""

import numpy as np
import time as timeMod
import os
import sys
import datetime
import subprocess as sub

from crane.data import dataContainer
from crane.data import timeArray
from data.loadHindcastStations import excludeNaNs,VALID_MIN
from data.extractStation import *

from files.buildPoints import BuildPoint

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def extractTransectForCoords( x, y, dataDir, varList, startTime, endTime, name, modelCoordSys='spcs' ) :
  """Extracts a transect defined by x,y coordinates and stores to disk in netCDF format."""
  dcs = []
  for var in varList :
    try :
      ee = extractTransect(dataDir, var, firstStack=1, modelCoordSys=modelCoordSys)
      ee.setTransect( name, x, y )
      dc = ee.extractDates( startTime, endTime )
      dcs.append( dc )
    except Exception as e :
      print 'Extraction failed'
      print e
  return dcs


def extractTransectForBPFile( bpFile, dataDir, varList, startTime, endTime, name, modelCoordSys='spcs' ) :
  """Extracts a transect defined in the bpFile and stores to disk in netCDF format."""
  bpObj = BuildPoint()
  bpObj.readFileFromDisk(bpFile)
  x = bpObj.points[:,1]
  y = bpObj.points[:,2]

  return extractTransectForCoords( x, y, dataDir, varList, startTime, endTime,
                                   name, modelCoordSys )
  
#-------------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------
class extractTransect(extractBase) :
  """A higher lever extraction object for transects"""
  def __init__(self, dataDir, fieldName, firstStack=1, modelCoordSys='spcs') :
    extractBase.__init__(self,dataDir,fieldName,firstStack,modelCoordSys)
    self.extractor.setTransectMode()
    self.name = None
    
  def setTransect( self, name, X, Y ) :
    '''Transect coordinates. All points outside the grid are removed.'''
    self.name = name
    if X.shape != Y.shape :
      raise Exception( 'X and Y must have same shape' )
    # This automatically discards outside coordinates in the extractor
    pointsInGrid = self.setCoordinates( X,Y )
    badPoints = [ s for i, s in enumerate(zip(X,Y)) if not pointsInGrid[i] ]
    for pointXY in badPoints:
      print 'Point out of grid, removing: ', pointXY
    if (pointsInGrid == 0).all() :
      print 'Warning: no points found in grid.'
    mask = pointsInGrid.copy().astype(bool)
    self.X = X[mask] # always in spcs coords
    self.Y = Y[mask]

  def setCoordinates( self, X, Y ) :
    x = X.copy()
    y = Y.copy()
    if self.modelCoordSys == 'utm' :
      # convert stations to utm for extraction
      for i in range(len(X)) :
        x[i],y[i] = spcs2utm( x[i], y[i] )
    return self.extractor.setCoordinates( x,y )

  def extract(self, stacks) :
    """Extract data for given stacks. Returns the transect data in dataContainer."""
    t,d = self.extractor.extract(stacks)
    if t == []:
      return []
    var = self.fieldName
    xUnique = self.X.flatten() # array where each x appears only once
    yUnique = self.Y.flatten()
    x = []
    y = []
    z = []
    data = []
    nComponents = self.extractor.getNumberOfComponents()
    for i in range(len(xUnique)) :
      # all zCoords for this x,y point, versus time
      zPoint = d[2,:,i,:] # (dim1, dim2, z), nvrt, nxy, ntime
      # remove bad values (below bottom)
      #goodIxZ = np.logical_and( np.logical_not( np.isnan( np.sum(zPoint,axis=1) ) ), 
                                #np.logical_not( np.isnan( np.sum(d[0,:,i,:],axis=1) ) ) )
      goodIxZ = np.logical_not( np.isnan( np.sum(zPoint,axis=1) ) )
      if np.all( goodIxZ == False ) :
        continue
      zPoint = zPoint[ goodIxZ, : ]
      # all data values for this x,y point, versus time
      dPoint = d[:nComponents,goodIxZ,i,:] # (dim1, dim2), nvrt, ntime
      dPoint = dPoint.swapaxes(0,1) # from (dim,nvrt,time) to (nvrt,dim,time)
      tmp = xUnique[i]*np.ones( (zPoint.shape[0],) )
      x.append( tmp )
      tmp = yUnique[i]*np.ones( (zPoint.shape[0],) )
      y.append( tmp )
      z.append( zPoint )
      data.append( dPoint )
    # TODO is this efficient?? len(x)~10000
    x = np.hstack( tuple(x) ) # shape (npt,)
    y = np.hstack( tuple(y) ) # shape (npt,)
    z = np.vstack( tuple(z) ) # shape (npt,ntime)
    data = np.vstack( tuple(data) )
    # create dataContainer
    ta = timeArray(t, 'simulation', self.extractor.startTime).asEpoch()
    # if suspected bad values, print warning
    hasBadValues = np.isnan(data).any() or np.isinf(data).any() or np.any( data < VALID_MIN )
    if hasBadValues :
      print 'Warning: bad values in', self.name
    meta = {}
    meta['location'] = self.name
    meta['instrument'] = 'model'
    meta['bracket'] = 'A'
    meta['variable'] = var
    meta['dataType'] = 'transect'
    dc = dataContainer('', ta, x,y,z, data, fieldNameList[var],
                       coordSys='spcs', metaData=meta)
    return dc


#-------------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------------
def parseCommandLine() :
  from optparse import OptionParser

  parser = OptionParser()
  parser.add_option('-r', '--runTag', action='store', type='string',
                      dest='runTag', help='Run tag, used as a label in post-proc.')
  parser.add_option('-d', '--dataDirectory', action='store', type='string',
                      dest='dataDir', help='directory where model outputs are stored')
  parser.add_option('-C', '--read-netcdf', action='store_true',
                    dest='readNetcdf',
                    help='Extract from SELFE netcdf output files instead of SELFE binary files (default %default)',default=False)
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
  parser.add_option('', '--stacks', action='store', type='string',
                    dest='stackStr', help='range of output files to read '
                    '(e.g 1,14) if start,end not given')
  parser.add_option('-v', '--variable', action='store', type='string',
                      dest='varList', help='variable(s) to extract: elev,temp,salt, ...\nTo use specific output file define extrension, e.g. salt.70')
  parser.add_option('-t', '--buildPointFile', action='store', type='string',
                      dest='bpFile', help='text file (*.bp) containing the transect (x,y) coordinates',
                             default=None)
  parser.add_option('-n', '--name', action='store', type='string',
                      dest='name', help='name of the transect (Default: first line in *.bp file)')
  parser.add_option('-c', '--modelCoordSys', action='store', type='string',
                      dest='modelCoordSys', default='spcs',
                      help='horizontal coordinate system used in model: '
                           'spcs or utm (Default: %default)')
  parser.add_option('-o', '--outDirectory', action='store', type='string',
                      dest='outDir', help='base directory for netCDF file tree '
                                                '(optional)')
  parser.add_option('-T', '--tracerModel', action='store', type='string', dest='tracerModel',
                    help='Enable extraction of tracers: sed, oxy, generic. Must '
                         'supply number of tracers for \'sed\' and \'generic\' '
                         'models via the -N switch. \'oxy\' model provides tracers: '
                         '\'NO3\',\'NH4\',\'phy\',\'zoo\',\'det\' and \'oxy\'.', default=None)
  parser.add_option('-N', '--numTracers', action='store', type='int', dest='numTracers',
                    help='Number of tracers to support for \'sed\' and \'generic\' models',
                    default=None)
  parser.add_option('', '--decimals', action='store', type='int', dest='digits',
                    help='Round extracted data to given decimal precision to save disk space', default=None)

  (options, args) = parser.parse_args()

  dataDir       = options.dataDir
  varList       = options.varList.split(',') if options.varList else None
  bpFile        = options.bpFile
  name          = options.name
  modelCoordSys = options.modelCoordSys
  outDir        = options.outDir
  runTag        = options.runTag
  startStr      = options.startStr
  endStr        = options.endStr
  stackStr      = options.stackStr
  readNetcdf = options.readNetcdf
  runTag        = options.runTag
  tracerModel   = options.tracerModel
  numTracers    = options.numTracers
  digits        = options.digits

  if tracerModel :
    if not numTracers and tracerModel.split('.')[0] in ['sed','generic']:
      parser.print_help()
      parser.error('numTracers must be provided if sed or generic tracer models are used.')
    addTracers( tracerModel, varList, numTracers=numTracers)

  if not dataDir :
    parser.print_help()
    parser.error('dataDir  undefined')
  if not varList and tracerModel is None :
    parser.print_help()
    parser.error('variable undefined')
  #if not outDir :
    #parser.print_help()
    #parser.error('outDir   undefined')
  if startStr is None and stackStr is None:
    parser.print_help()
    parser.error('startStr undefined')
  if endStr is None and stackStr is None:
    parser.print_help()
    parser.error('endStr   undefined')
  if not runTag :
    parser.print_help()
    parser.error('runTag  undefined')
  if not bpFile :
    parser.print_help()
    parser.error('buildPointFile  undefined')

  if stackStr is not None:
    limits = [int(v) for v in stackStr.split(',')]
    stacks = np.arange(limits[0], limits[1]+1)
  else:
    stacks = None

  if startStr is not None:
    startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
    endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')
  else:
    startTime = endTime = None

  print 'Parsed options:'
  if name :
    print ' - name ', name
  if stackStr is None:
    print ' - time range:',str(startTime),'->', str(endTime)
  else:
    print ' - stacks:',str(stacks[0]),'->', str(stacks[-1])
  print ' - dataDir',dataDir
  print ' - SELFE output format:','netCDF' if readNetcdf else 'binary'
  print ' - runTag',runTag
  print ' - transect file',bpFile
  if outDir :
    print ' - output dir',outDir
  print ' - variables ', varList
  print ' - model coord system', modelCoordSys
  sys.stdout.flush()

  dcs = []
  if readNetcdf :
    # use ncExtract instead of this module
    from data.ncExtract import extractTransectForBPFile as extractTransectForBPFileNC
    dcs = extractTransectForBPFileNC(bpFile, dataDir, varList,
                                     startTime, endTime, name=name,
                                     modelCoordSys=modelCoordSys, stacks=stacks)
  else :
    dcs = extractTransectForBPFile( bpFile, dataDir, varList,
                                    startTime, endTime, name=name,
                                    modelCoordSys=modelCoordSys )

  for dc in dcs :
    dc.setMetaData( 'tag', runTag )
  import data.dirTreeManager as dtm
  rule = dtm.oldTreeRule()
  #rule = dtm.defaultTreeRule()
  dtm.saveDataContainerInTree( dcs, path=outDir, rule=rule, dtype=np.float32,
                               overwrite=True, compress=True, digits=digits )

if __name__=='__main__' :
  parseCommandLine()
