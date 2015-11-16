"""
Computes plume volume metrics from SELFE outputs.

Reads *_salt.63.nc files. Only netcdf output files are supported.

Tuomas Karna 2013-11-08
"""
import sys
import numpy as np
import datetime
import traceback
import time as timeMod
from netCDF4 import Dataset as NetCDFFile

from crane.data import meshContainer
from crane.data import dataContainer
from crane.data import timeArray
from crane.files import gr3Interface
from crane.data import ncExtract
from crane.data import dirTreeManager


class plumeStatsComputer(object) :
  """Class that computes the default plume statistics from SELFE netcdf
  outputs.
  """
  def __init__(self, path, plumeRegionGR3File, saltThreshold) :
    """Initializes data structures for the given plume region file."""
    # read salt field from netcdf file
    self.extractor = ncExtract.selfeExtractBase(path, var='salt')
    self.extractor.initialize()
    self.ncReader = self.extractor.dataFile
    self.saltThreshold = saltThreshold

    # load plume region
    plumeRegionMC = gr3Interface.readGR3FileToMC(plumeRegionGR3File)
    # check that gr3 file matches the mesh in netCDF file
    if ( self.ncReader.nodeX.shape != plumeRegionMC.x.shape or
         not np.allclose(self.ncReader.nodeX,plumeRegionMC.x) or
         self.ncReader.nodeY.shape != plumeRegionMC.y.shape or
         not np.allclose(self.ncReader.nodeY,plumeRegionMC.y) or
         not np.array_equal(self.ncReader.faceNodes,
                            plumeRegionMC.connectivity) ) :
        raise Exception('Given GR3 file does not match the mesh in output files')
    plume_mask = plumeRegionMC.data[:,0,0] == 1

    # create minimesh only for plume region, save 50% cpu time
    goodElems = plume_mask[self.ncReader.faceNodes].max(axis=1)
    self.plume_mask = np.sort(np.unique(self.ncReader.faceNodes[goodElems,:]))
    self.nodeX = self.ncReader.nodeX[self.plume_mask]
    self.nodeY = self.ncReader.nodeY[self.plume_mask]
    #print self.nodeX.shape
    self.nodesToMiniMesh = -1*np.ones_like(self.ncReader.nodeX,dtype=int)
    #print self.nodeX.shape
    self.nodesToMiniMesh[self.plume_mask] = np.arange(len(self.nodeX))
    self.faceNodes = self.nodesToMiniMesh[ self.ncReader.faceNodes ]
    self.faceNodes = self.faceNodes[goodElems]

    self.areas = meshContainer.computeAreas(self.faceNodes,self.nodeX,self.nodeY)
    self.elem_center_x = self.nodeX[self.faceNodes].mean(axis=1)
    self.elem_center_y = self.nodeY[self.faceNodes].mean(axis=1)
    self.bath = self.ncReader.bath[self.plume_mask]

    ## full mesh
    #areas = computeAreas(ee.faceNodes,ee.nodeX,ee.nodeY)
    #elem_center_x = ee.nodeX[ee.faceNodes].mean(axis=1)
    #elem_center_y =ee.nodeY[ee.faceNodes].mean(axis=1)
    #bath = ee.bath
    #faceNodes = ee.faceNodes

  def processStack( self, stack ):
    """Computes plume stats for all time steps in the netcdf file
    
    Returns
    -------
    times : np.ndarray (nTime,)
      time stamps in epoch format
    values : np.ndarray (nTime,nStats)
      plume statistics: area, centroid_x, centroid_y, volume, thickness
    """
    ncfile = self.extractor.getNCFile(stack)
    nTime = len(ncfile.dimensions['time'])
    time = ncfile.getTime()
    # area, centroid_x, centroid_y, volume, thickness
    plume_stats = np.zeros((nTime,5))
    for iTime in xrange(nTime) :

      # minimesh
      salt = ncfile.variables['salt'][iTime,:,:][:,self.plume_mask] #(nTime,nVert,nNodes)
      eta = ncfile.variables['elev'][iTime,:][self.plume_mask]
      ## full mesh
      #salt = ncfile.variables['salt'][iTime,:,:] #(nTime,nVert,nNodes)
      #eta = ncfile.variables['elev'][iTime,:]

      Z, kbp2, iwet = self.ncReader.vCoords.computeVerticalCoordinates(eta,self.bath)
      salt_in_range = salt <= self.saltThreshold

      # minimesh
      plume_nodes = np.logical_and( salt_in_range.max(axis=0), iwet )
      ##full mesh
      #plume_nodes = np.logical_and( salt_in_range.max(axis=0), plume_mask )
      #plume_nodes = np.logical_and( plume_nodes, iwet )

      plume_triangles = plume_nodes[self.faceNodes].max(axis=1)
      plume_area = np.sum(self.areas[plume_triangles])
      plume_x = np.sum((self.areas*self.elem_center_x)[plume_triangles])/plume_area
      plume_y = np.sum((self.areas*self.elem_center_y)[plume_triangles])/plume_area

      # alpha in [0,1]: fraction of vertical edge in the plume
      salt_top = salt[1:,:] # top of each vertical 1d element
      salt_bot = salt[:-1,:]
      salt_min=np.minimum(salt_top,salt_bot)
      a=(self.saltThreshold-salt_min)
      b=np.abs(salt_top - salt_bot)
      alpha=np.zeros_like(a)
      ix = np.logical_and( b > 1e-10, ~b.mask )
      alpha[ix]=a[ix]/b[ix]
      alpha[alpha<0] = 0
      alpha[alpha>1] = 1

      # compute height in plume for each vertical line
      z_top = Z[1:,:]
      z_bot = Z[:-1,:]
      plume_height_nodes = (alpha*(z_top-z_bot)).sum(axis=0)
      plume_height_tri = plume_height_nodes[self.faceNodes].mean(axis=1)

      # et finalement le plume volume
      plume_volume = self.areas*plume_height_tri
      plume_volume = plume_volume[plume_triangles].sum()

      # plume sickness
      plume_thickness = plume_volume/plume_area
      plume_stats[iTime,:] = [plume_area,plume_x,plume_y,plume_volume,plume_thickness]

      #print 'plume area',plume_area
      #print 'plume centroid:',plume_x,plume_y
      #print 'plume volume:',plume_volume
      #print 'plume thickness:',plume_thickness

    return time, plume_stats

  def processStacks( self, stacks ) :
    """Computes plume statistics for all the given stacks
    
    Returns
    -------
    times : np.ndarray (nTime,)
      time stamps in epoch format
    values : np.ndarray (nTime,nStats)
      plume statistics: area, centroid_x, centroid_y, volume, thickness
    """
    times=[]
    values=[]
    cputime0 = timeMod.clock()
    
    sys.stdout.write(' * Processing stacks '+str(stacks[0])+' - ' +str(stacks[-1]))
    sys.stdout.flush()
    for i in range(len(stacks)) :
      try :
        stack = stacks[i]
        # (nTime,) (nTime,nStats)
        ti, vi = self.processStack( stack )
        times.append(ti)
        values.append(vi)
      except Exception as e:
        print 'computing plume stats failed'
        traceback.print_exc(file=sys.stdout)

    times = np.ma.concatenate( tuple(times), axis=0 )
    values = np.ma.concatenate( tuple(values), axis=0 )
    sys.stdout.write( ' duration %.2f s\n'%(timeMod.clock()-cputime0) )

    return times, values
  
  def processDates( self, startTime, endTime ) :
    """Computes plume volume for the given date range
    
    Returns
    -------
    times : np.ndarray (nTime,)
      time stamps in epoch format
    values : np.ndarray (nTime,nStats)
      plume statistics: area, centroid_x, centroid_y, volume, thickness
    """
    stacks = self.ncReader.getStacks(startTime, endTime, wholeDays=True)
    return self.processStacks(stacks)
  
  def getDataContainer( self, startTime, endTime ) :
    """Computes plume stats for the given time period. Returns a dataContainer
    for each metric.

    Returns
    -------
    times : np.ndarray (nTime,)
      time stamps in epoch format
    values : np.ndarray (nTime,nStats)
      plume statistics: area, centroid_x, centroid_y, volume, thickness
    """
    time,values = self.processDates(startTime,endTime)
    # make dataContainer
    goodIx = np.isfinite(np.sum(values,axis=1))
    time = time[goodIx]
    values = values[goodIx,:]
    varNames = ['plume_area','plume_center',
                'plume_volume','plume_thickness']
    fieldIndices = [[0],[1,2], [3],[4]]
    fieldNames = {'plume_center':['plume_center_x','plume_center_y']}
    dcs=[]
    for i,var in enumerate(varNames):
      sthSuffix = '_{0:d}'.format(int(self.saltThreshold))
      data = np.swapaxes(values[:,fieldIndices[i]],0,1)[None,:,:] # (1,nStats,nTime)
      ta = timeArray.timeArray( time, 'epoch' )
      meta = {}
      meta['location'] = 'plume'
      meta['instrument'] = 'model'
      meta['variable'] = var+sthSuffix
      meta['dataType'] = 'plumemetrics'
      meta['saltThreshold'] = str(self.saltThreshold)
      x = y = z = 0
      fNames = [ fn+sthSuffix for fn in fieldNames.get(var,[var]) ]
      dc = dataContainer.dataContainer('', ta, x,y,z, data, fNames,
                          coordSys='spcs',metaData=meta)
      dcs.append(dc)
    return dcs

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
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
  parser.add_option('-v', '--saltThreshold', action='store', type='float',
                    default=28.0, dest='saltThreshold',
                    help='salinity value that defines the plume, S<=threshold (default %default)')
  parser.add_option('-t', '--plumeRegionFile', action='store', type='string',
                      dest='plumeRegionGR3File', help='a gr3 file indicating the plume region (depth==1)',
                             default=None)

  (options, args) = parser.parse_args()

  runTag        = options.runTag
  startStr      = options.startStr
  endStr        = options.endStr
  dataDir       = options.dataDir
  plumeRegionGR3File = options.plumeRegionGR3File
  saltThreshold = options.saltThreshold

  if not dataDir :
    parser.print_help()
    parser.error('dataDir  undefined')
  if not startStr :
    parser.print_help()
    parser.error('startStr undefined')
  if not endStr :
    parser.print_help()
    parser.error('endStr   undefined')
  if not runTag :
    parser.print_help()
    parser.error('runTag  undefined')
  if not plumeRegionGR3File :
    parser.print_help()
    parser.error('plume region file undefined undefined')

  startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
  endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')

  print 'Parsed options:'
  print ' - time range:',str(startTime),'->', str(endTime)
  print ' - salinity threshold:',saltThreshold
  print ' - dataDir',dataDir
  print ' - runTag',runTag
  print ' - plume region file',plumeRegionGR3File

  # Extract 
  psc = plumeStatsComputer( dataDir, plumeRegionGR3File, saltThreshold )
  dcs = psc.getDataContainer( startTime, endTime )
  for dc in dcs :
    dc.setMetaData('tag',runTag)
  
  rule = 'monthlyFile'
  dirTreeManager.saveDataContainerInTree( dcs, rule=rule, dtype=np.float32,
                               overwrite=True )

if __name__=='__main__' :
  parseCommandLine()
