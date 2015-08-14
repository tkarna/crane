"""
Computes regional volume metrics from SELFE *.63.nc outputs.

Generalized version of computePlumeMetrics.py
"""
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from data.meshContainer import *
import files.gr3Interface as gr3Interface
from netCDF4 import Dataset as NetCDFFile
from data.ncExtract import *
import time as timeMod


#-----------------------------------------------------------------------------
# Classes
#-----------------------------------------------------------------------------
class volumeComputer(ncExtractBase):
  """Volumetric computations from SELFE netCDF outputs"""
  def __init__(self, path, volumeFile, volumeVar,
               volumeThresholdLow=None, volumeThresholdHi=None,
               criteriaFile=None, criteriaVar=None,
               criteriaThresholdLow=None, criteriaThresholdHi=None,
               regionFile=None, regionName='region'):
    """Initialize data structures for region based volumetric computations.

    Params
    ------
    path : Path to input netCDF files
    volumeFile : File name for variable to compute volume (*.63.nc)
    volumeVar : Variable to compute volume
    volumeThresholdLow : Compute volume above this threshold (inclusive)
    volumeThresholdHi : Compute volume below this threshold (inclusive)
    criteriaFile : Variable to constrain volume computation (*.63.nc)
    criteriaVar : Variable to constrain volume computation (*.63.nc)
    criteriaThresholdLow : Compute volume above this threshold (inclusive)
    criteriaThresholdHi : Compute volume below this threshold (inclusive)
    region : Path to file that defines region
    """
    self.volumeFile = volumeFile
    self.volumeVar = volumeVar
    self.volumeFileType = self.IdFileType(volumeFile)
    self.volumeThresholdLow = volumeThresholdLow
    self.volumeThresholdHi = volumeThresholdHi
    self.criteriaFile = criteriaFile
    self.criteriaVar = criteriaVar
    if criteriaFile is not None:
      self.volumeFileType = self.IdFileType(criteriaFile)
    self.criteriaThresholdLow = criteriaThresholdLow
    self.criteriaThresholdHi = criteriaThresholdHi
    self.regionFile = regionFile
    self.regionName = regionName

    # Read headers from netcdf files
    self.volumeNCReader = ncExtractBase(path, volumeFile)
    self.volumeNCReader.readHeader()

    if criteriaFile is not None:
      self.criteriaNCReader = ncExtractBase(path, criteriaFile)
      self.criteriaNCReader.readHeader()
    else:
      self.criteriaNCReader = None

    # Create mask from region file and create mini-mesh
    # Or, create a fake mask with every node included
    if regionFile is not None:
      self.region_mask = self.createRegionNodeMask(regionFile)
      self.createRegionMesh(self.region_mask)
    else:
      whole_mesh = np.ones_like(self.volumeNCReader.nodeX, dtype=bool)
      self.createRegionMesh(whole_mesh)

  def IdFileType(self, filename):
    """Identify and return file type."""
    # variable.type.nc
    type = filename.split('.')[1] 
    return type

  def createRegionNodeMask(self, region_file):
    """Load region file and check that it matches netcdf file

    Params
    ------
    region_file : Path to file that defines the region

    Returns
    -------
    mask : Masked defining region
    """
    region = gr3Interface.readGR3FileToMC(region_file)
    self.region = region
    # check that gr3 file matches the mesh in netCDF file
    if ( self.volumeNCReader.nodeX.shape != region.x.shape or
         not np.allclose(self.volumeNCReader.nodeX, region.x) or
         self.volumeNCReader.nodeY.shape != region.y.shape or
         not np.allclose(self.volumeNCReader.nodeY, region.y) or
         not np.array_equal(self.volumeNCReader.faceNodes,
                            region.connectivity) ) :
        raise Exception('Given GR3 file does not match the mesh in output files')
    mask = region.data[:, 0, 0] == 1
    print ' - %d nodes in region.' % np.sum(mask)

    return mask

  def createRegionElemMask(self, region_node_mask):
    """Load region file and check that it matches netcdf file

    Params
    ------
    region_mask : Nodal region mask 

    Returns
    -------
    mask : Masked defining region by elements
    """
    elem_nodes = self.volumeNCReader.faceNodes 
    elem_mask = np.alltrue(region_node_mask[elem_nodes], axis=1)
    print ' - %d elements in region.' % np.sum(elem_mask)

    return elem_mask

  def createRegionMesh(self, regionMask):
    """Create sub-mesh defined by region to speed calculations."""
    # create minimesh only for region, save 50% cpu time
    # Includes any element with a node in region
    #goodElems = regionMask[self.volumeNCReader.faceNodes].max(axis=1)
    # Includes only elements with all nodes in region
    goodElems = np.alltrue(regionMask[self.volumeNCReader.faceNodes], axis=1)
    self.region_mask = np.sort(np.unique(self.volumeNCReader.faceNodes[goodElems,:]))
    self.nodeX = self.volumeNCReader.nodeX[self.region_mask]
    self.nodeY = self.volumeNCReader.nodeY[self.region_mask]
    #print self.nodeX.shape
    self.nodesToMiniMesh = -1*np.ones_like(self.volumeNCReader.nodeX,dtype=int)
    #print self.nodeX.shape
    self.nodesToMiniMesh[self.region_mask] = np.arange(len(self.nodeX))
    self.faceNodes = self.nodesToMiniMesh[ self.volumeNCReader.faceNodes ]
    self.faceNodes = self.faceNodes[goodElems]

    self.areas = computeAreas(self.faceNodes,self.nodeX,self.nodeY)
    self.elem_center_x = self.nodeX[self.faceNodes].mean(axis=1)
    self.elem_center_y = self.nodeY[self.faceNodes].mean(axis=1)
    self.bath = self.volumeNCReader.bath[self.region_mask]

  def calculateVolume(self, timeStep, volFile, criteriaFile=None):
    """ Calculate the volume of a water mass based on thresholds from 1 or 2 variables

    Params
    ------
    timeStep : Time step of array
    volFile : NCFile of variable to constrain volume calculation
    criteriaFile : NCFile of second variable to constrain volume calculation

    Returns
    -------
    region_stats : [region_area, region_x, region_y, region_volume, region_conc, region_thickness]
    """
    if self.volumeFileType == '63':
      region_stats = self.calculateNodalVolume(timeStep, volFile, criteriaFile)
    elif self.volumeFileType == '70':
      region_stats = self.calculateElementVolume(timeStep, volFile, criteriaFile)
    else:
      raise Exception('Given GR3 file does not match the mesh in output files')

    return region_stats

  def calculateNodalVolume(self, timeStep, volFile, criteriaFile):
    """ Calculate the volume of a water mass based on thresholds from 1 or 2 nodal variables

    Params
    ------
    timeStep : Time step of array
    volFile : NCFile of variable to constrain volume calculation
    criteriaFile : NCFile of second variable to constrain volume calculation

    Returns
    -------
    region_stats : [region_area, region_x, region_y, region_volume, region_conc, region_thickness]
    """
    # Load vars
    if self.regionFile is not None:
      volVar = volFile.variables[self.volumeVar][timeStep, :, :][:, self.region_mask]
      eta = volFile.variables['elev'][timeStep, :][self.region_mask]
      if criteriaFile is not None:
        criteria_var = criteriaFile.variables[self.criteriaVar][timeStep, :, :][:, self.region_mask]
    else:
      volVar = volFile.variables[self.volumeVar][timeStep, :, :][:, :]
      eta = volFile.variables['elev'][timeStep, :][:]
      if criteriaFile is not None:
        criteria_var = criteriaFile.variables[self.criteriaVar][timeStep, :, :][:, :]

    # Prepare vertical
    Z, kbp2, iwet = self.volumeNCReader.vCoords.computeVerticalCoordinates(eta, self.bath)

    # Apply thresholds
    if self.volumeThresholdLow is not None or self.volumeThresholdHi is not None:
      vol_low = np.ones_like(volVar, dtype=bool)
      vol_hi = np.ones_like(volVar, dtype=bool)
      if self.volumeThresholdLow is not None:
        print '- filter out %s below %f ' % (self.volumeVar, self.volumeThresholdLow)
        vol_low = volVar > self.volumeThresholdLow
      if self.volumeThresholdHi is not None:
        print '- filter out %s above %f ' % (self.volumeVar, self.volumeThresholdHi)
        vol_hi = volVar < self.volumeThresholdHi
      vol_in_range = np.logical_and(vol_low, vol_hi)
    else:
      vol_in_range = np.ones_like(volVar, dtype=bool)

    if criteriaFile is not None:
      criteria_in_range = criteria_var.copy()
      cri_lo = np.ones_like(criteria_in_range, dtype=bool)
      cri_hi = np.ones_like(criteria_in_range, dtype=bool)
      if self.criteriaThresholdLow is not None:
        print '- filter out %s below %f ' % (self.criteriaVar, self.criteriaThresholdLow)
        cri_lo = criteria_var > self.criteriaThresholdLow
      if self.criteriaThresholdHi is not None:
        print '- filter out %s above %f ' % (self.criteriaVar, self.criteriaThresholdHi)
        cri_hi = criteria_var < self.criteriaThresholdHi
      criteria_in_range = np.logical_and(cri_lo, cri_hi)
      vol_in_range = np.logical_and(criteria_in_range, vol_in_range)

    # Remove "dry" values identified with value of -99.0
    wet_values = volVar > -99.0
    vol_in_range = np.logical_and(wet_values, vol_in_range)

    # Calculate region values
    region_nodes = np.logical_and(vol_in_range.max(axis=0), iwet)
    region_triangles = region_nodes[self.faceNodes].max(axis=1)
    region_area = np.sum(self.areas[region_triangles])
    region_x = np.sum((self.areas*self.elem_center_x)[region_triangles])/region_area
    region_y = np.sum((self.areas*self.elem_center_y)[region_triangles])/region_area
    region_var_nodes = np.ma.array(volVar, mask=(np.tile(~region_nodes, (volVar.shape[0],1))))

    # Calculate element concentrations - mean of 6 element nodes
    nLevels = region_var_nodes.shape[0]
    nElems = len(region_triangles)
    region_var_tris = np.zeros((nLevels, nElems))
    for i in np.arange(nLevels):
      region_var_tris[i, :] = region_var_nodes[i, :][self.faceNodes].mean(axis=1)
    region_var_tris = np.ma.masked_where(region_var_tris == 1.0000002e+20, region_var_tris)
    region_var_tris_top = region_var_tris[1:, :]
    region_var_tris_bot = region_var_tris[:-1, :]
    region_var_elems = (region_var_tris_top+region_var_tris_bot)/2.0

    # Vertical fractions of criteria
    alpha = np.ones_like(np.minimum(volVar[1:, :], volVar[:-1,:]))
    if self.volumeThresholdLow or self.volumeThresholdHi:
      vol_var_top = volVar[1:, :]
      vol_var_bot = volVar[:-1, :]
      vol_min = np.minimum(vol_var_top, vol_var_bot)
      vol_max = np.maximum(vol_var_top, vol_var_bot)
      dvar_dz = np.abs(vol_var_top - vol_var_bot)
      #ix = np.logical_and(dvar_dz > 1e-10, ~dvar_dz.mask)
      ix = dvar_dz > 1e-10
      alpha_hi = np.ones_like(alpha)
      alpha_low = np.ones_like(alpha)

      # Criteria hi (Filter out values above)
      if self.volumeThresholdHi:
        a = self.volumeThresholdHi - vol_min
        alpha_hi[ix] = a[ix]/dvar_dz[ix]
        alpha_hi[alpha_hi < 0] = 0
        alpha_hi[alpha_hi > 1] = 1

      # Criteria lo (Filter out values below)
      if self.volumeThresholdLow:
        a = vol_max - self.volumeThresholdLow
        alpha_low[ix] = a[ix]/dvar_dz[ix]
        alpha_low[alpha_low < 0] = 0
        alpha_low[alpha_low > 1] = 1

      alpha = np.minimum(alpha_low, alpha_hi)

    # Filter by secondary variable
    if criteriaFile is not None:
      var_top = criteria_var[1:, :]
      var_bot = criteria_var[:-1, :]
      var_min = np.minimum(var_top, var_bot)
      var_max = np.maximum(var_top, var_bot)
      dvar_dz = np.abs(var_top - var_bot)
      #ix = np.logical_and(dvar_dz > 1e-10, ~dvar_dz.mask)
      ix = dvar_dz > 1e-10
      alpha_var = np.ones_like(var_min)
      alpha_hi = np.ones_like(alpha)
      alpha_low = np.ones_like(alpha)

      # Criteria hi (Filter out values above)
      if self.criteriaThresholdHi:
        a = self.criteriaThresholdHi - var_min
        alpha_hi[ix] = a[ix]/dvar_dz[ix]
        alpha_hi[alpha_hi < 0] = 0
        alpha_hi[alpha_hi > 1] = 1

      # Criteria lo (Filter out values below)
      if self.criteriaThresholdLow:
        a = var_max - self.criteriaThresholdLow
        alpha_low[ix] = a[ix]/dvar_dz[ix]
        alpha_low[alpha_low < 0] = 0
        alpha_low[alpha_low > 1] = 1

      alpha_var = np.minimum(alpha_low, alpha_hi)
      alpha = np.minimum(alpha, alpha_var)

    # Compute height in region for each prism and vertical line
    z_top = Z[1:, :]
    z_bot = Z[:-1, :]
    region_height_nodes = alpha*(z_top - z_bot)
    region_height_nodes_total = region_height_nodes.sum(axis=0)
    region_height_tris = np.zeros((nLevels-1, nElems))
    for i in np.arange(nLevels-1):
#      region_height_tris[i, :] = region_height_nodes[i, :][self.faceNodes].mean(axis=1)
      region_height_tris[i, :] = region_height_nodes[i, :][self.faceNodes].mean(axis=1)
    region_height_tris_top = region_height_tris[1:, :]
    region_height_tris_bot = region_height_tris[:-1, :]
    region_height_elems = (region_height_tris_top+region_height_tris_bot)/2.0
    region_height_tri_total = region_height_nodes_total[self.faceNodes].mean(axis=1)

    # Region volume
    region_conc = self.areas*region_height_tris*region_var_elems
    region_conc_total = region_conc.sum()
    if np.sum(region_triangles) > 0:
      region_volume = self.areas*region_height_tris
      region_volume = region_volume[region_triangles]
      region_volume = region_volume.sum()
      region_thickness = region_volume/region_area
    else:
      region_volume = 0
      region_x = 0
      region_y = 0
      region_conc_total = 0
      region_thickness = 0

    region_stats = [region_area, region_x, region_y, region_volume, region_conc_total, region_thickness]
    return region_stats

  def calculateElementVolume(self, timeStep, volFile, criteriaFile):
    """ Calculate the volume of a water mass based on thresholds from 1 or 2 element variables

    Params
    ------
    timeStep : Time step of array
    volFile : NCFile of variable to constrain volume calculation
    criteriaFile : NCFile of second variable to constrain volume calculation

    Returns
    -------
    region_stats : [region_area, region_x, region_y, region_volume, region_conc, region_thickness]
    """
    eta = volFile.variables['elev'][timeStep, :][self.region_mask]
    Z, kbp2, iwet = self.volumeNCReader.vCoords.computeVerticalCoordinates(eta, self.bath)

    # Load vars
    if self.regionFile is not None:
      node_mask = self.createRegionNodeMask(self.regionFile)
      element_mask = self.createRegionElemMask(node_mask)
      volVar = volFile.variables[self.volumeVar][timeStep, :, :][:, element_mask]
      eta = volFile.variables['elev'][timeStep, :][element_mask]
      if criteriaFile is not None:
        criteria_var = criteriaFile.variables[self.criteriaVar][timeStep, :, :][:, element_mask]
    else:
      volVar = volFile.variables[self.volumeVar][timeStep, :, :][:, :]
      eta = volFile.variables['elev'][timeStep, :][:]
      if criteriaFile is not None:
        criteria_var = criteriaFile.variables[self.criteriaVar][timeStep, :, :][:, :]

    # Apply thresholds
    if self.volumeThresholdLow is not None or self.volumeThresholdHi is not None:
      vol_low = np.ones_like(volVar, dtype=bool)
      vol_hi = np.ones_like(volVar, dtype=bool)
      if self.volumeThresholdLow is not None:
        print '- filter %s below %f ' % (self.volumeVar, self.volumeThresholdLow)
        vol_low = volVar > self.volumeThresholdLow
      if self.volumeThresholdHi is not None:
        print '- filter %s above %f ' % (self.volumeVar, self.volumeThresholdHi)
        vol_hi = volVar < self.volumeThresholdHi
      vol_in_range = np.logical_and(vol_low, vol_hi)
    else:
      vol_in_range = np.ones_like(volVar, dtype=bool)

    if criteriaFile is not None:
      criteria_in_range = criteria_var.copy()
      cri_lo = np.ones_like(criteria_in_range, dtype=bool)
      cri_hi = np.ones_like(criteria_in_range, dtype=bool)
      if self.criteriaThresholdLow is not None:
        print '- filter %s below %f ' % (self.criteriaVar, self.criteriaThresholdLow)
        cri_lo = criteria_var > self.criteriaThresholdLow
      if self.criteriaThresholdHi is not None:
        print '- filter %s above %f ' % (self.criteriaVar, self.criteriaThresholdHi)
        cri_hi = criteria_var < self.criteriaThresholdHi
      criteria_in_range = np.logical_and(cri_lo, cri_hi)
      vol_in_range = np.logical_and(criteria_in_range, vol_in_range)

    # Remove "dry" values identified with value of -99.0
    wet_values = volVar > -99.0
    vol_in_range = np.logical_and(wet_values, vol_in_range)


    # Calculate region values
# - All nodes must be wet to include element
#    wet_elems = iwet[self.faceNodes].min(axis=1)
# - One node can be wet to include element
    wet_elems = iwet[self.faceNodes].max(axis=1)
    region_elems = np.logical_and(vol_in_range.max(axis=0), wet_elems)
    region_area = np.sum(self.areas[region_elems])
    region_x = np.sum(self.areas*self.elem_center_x)/region_area
    region_y = np.sum(self.areas*self.elem_center_y)/region_area

    # Compute height in region for each prism and vertical line
    nLevels = volVar.shape[0] 
    nElems = volVar.shape[1] 
    z_top = Z[1:, :]
    z_bot = Z[:-1, :]
    region_height_nodes = z_top - z_bot
    region_height_nodes_total = region_height_nodes.sum(axis=0)
    region_height_tris = np.zeros((nLevels-1, nElems))
    for i in np.arange(nLevels-1):
      region_height_tris[i, :] = region_height_nodes[i, :][self.faceNodes].mean(axis=1)
    region_height_tri_total = region_height_nodes_total[self.faceNodes].mean(axis=1)

    # Region volume - area*height*volume*binaryMask
    region_conc = self.areas*region_height_tris*volVar[1:, :]*vol_in_range[1:, :]
    region_conc_total = region_conc.sum()
    if region_area > 0:
      region_volume = self.areas*region_height_tris*vol_in_range[1:, :]
      #region_volume = region_volume[region_triangles]
      region_volume = region_volume.sum()
      region_thickness = region_volume/region_area
    else:
      region_volume = 0
      region_x = 0
      region_y = 0
      region_conc_total = 0
      region_thickness = 0

    region_stats = [region_area, region_x, region_y, region_volume, region_conc_total, region_thickness]
    return region_stats
  
  def processStack( self, stack ):
    """Computes region stats for all time steps in the netcdf file

    Returns
    -------
    times : np.ndarray (nTime,)
      time stamps in epoch format
    values : np.ndarray (nTime,nStats)
      region statistics: area, centroid_x, centroid_y, volume, thickness
    """
    print '*** Processing stack %s ***' % stack
    ncfile = self.volumeNCReader.getNCFile(stack)
    if self.criteriaNCReader is not None:
      criteria_file = self.criteriaNCReader.getNCFile(stack)
    else:
      criteria_file = None
    nTime = len(ncfile.dimensions['time'])
    time = self.volumeNCReader.getTime(ncfile)
    # area, centroid_x, centroid_y, volume, integrated_concentration, thickness
    region_stats = np.zeros((nTime,6))
    for iTime in xrange(nTime) :
      region_stats[iTime, :] = self.calculateVolume(iTime, ncfile, criteria_file)
      print 'region area', region_stats[iTime, 0]
      print 'region centroid:', region_stats[iTime, 1], region_stats[iTime, 2]
      print 'region volume:', region_stats[iTime, 3]
      print 'region integrated concentration:', region_stats[iTime, 4]
      print 'region thickness:', region_stats[iTime, 5]

    return time, region_stats

  def processStacks( self, stacks ) :
    """Computes region statistics for all the given stacks

    Returns
    -------
    times : np.ndarray (nTime,)
      time stamps in epoch format
    values : np.ndarray (nTime,nStats)
      region statistics: area, centroid_x, centroid_y, volume, int_conc, thickness
    """
    times=[]
    values=[]
    cputime0 = timeMod.clock()

    sys.stdout.write(' * Processing stacks '+str(stacks[0])+' - '+str(stacks[-1])+'\n')
    sys.stdout.flush()
    for stack in stacks: 
      try :
        # (nTime,) (nTime,nStats)
        ti, vi = self.processStack( stack )
        times.append(ti)
        values.append(vi)
      except Exception as e:
        print 'computing region stats failed'
        traceback.print_exc(file=sys.stdout)

    times = np.ma.concatenate( tuple(times), axis=0 )
    values = np.ma.concatenate( tuple(values), axis=0 )
    sys.stdout.write( ' duration %.2f s\n'%(timeMod.clock()-cputime0) )

    return times, values

  def processDates( self, startTime, endTime ) :
    """Computes region volume for the given date range

    Returns
    -------
    times : np.ndarray (nTime,)
      time stamps in epoch format
    values : np.ndarray (nTime,nStats)
      region statistics: area, centroid_x, centroid_y, volume, thickness
    """
    stacks = self.volumeNCReader.getStacks( startTime, endTime )
    return self.processStacks( stacks )

  def getDataContainer( self, startTime, endTime ) :
    """Computes region stats for the given time period. Returns a dataContainer
    for each metric.

    Returns
    -------
    times : np.ndarray (nTime,)
      time stamps in epoch format
    values : np.ndarray (nTime,nStats)
      region statistics: area, centroid_x, centroid_y, volume, thickness
    """
    time,values = self.processDates(startTime,endTime)
    # make dataContainer
    goodIx = np.isfinite(np.sum(values,axis=1))
    time = time[goodIx]
    values = values[goodIx,:]
    varNames = ['region_area','region_center',
                'region_volume','region_conc','region_thickness']
    fieldIndices = [[0],[1,2], [3],[4],[5]]
    fieldNames = {'region_center':['region_center_x','region_center_y']}
    dcs=[]
    for i,var in enumerate(varNames):
      data = np.swapaxes(values[:,fieldIndices[i]],0,1)[None,:,:] # (1,nStats,nTime)
      ta = timeArray( time, 'epoch' )
      meta = {}
      # Add region name
      region_desc = [self.volumeVar, self.volumeThresholdLow, self.volumeThresholdHi,
                     self.criteriaVar, self.criteriaThresholdLow, self.criteriaThresholdHi]
      desc = '_'.join([str(s) for s in region_desc if s is not None])
      meta['location'] = self.regionName 
      meta['instrument'] = 'model'
      meta['variable'] = var+'_'+desc
      meta['dataType'] = 'regionmetrics'
      meta['varLowThreshold'] = str(self.volumeThresholdLow)
      meta['varHiThreshold'] = str(self.volumeThresholdHi)
      meta['criteriaVar'] = str(self.criteriaVar)
      meta['criteriaThresholdLow'] = str(self.criteriaThresholdLow)
      meta['criteriaThresholdHi'] = str(self.criteriaThresholdHi)
      x = y = z = 0
      fNames = [ fn for fn in fieldNames.get(var,[var]) ]
      dc = dataContainer('', ta, x,y,z, data, fNames,
                          coordSys='spcs',metaData=meta)
      dcs.append(dc)
    return dcs

#------------------------------------------------------------------------------
# Main: Commandline interface
#------------------------------------------------------------------------------
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
  parser.add_option('-f', '--volumeFile', action='store', type='string',
                      dest='volFile', help='File with variable to calculate volume')
  parser.add_option('-v', '--volumeVar', action='store', type='string',
                      dest='volVar', help='Variable in volume file to calculate volume')
  parser.add_option('-l', '--volThreshLow', action='store', type='float',
                      dest='volThreshLow', help='Lower bound of value to describe region to calculate volume',
                      default=None)
  parser.add_option('-u', '--volThreshHi', action='store', type='float',
                      dest='volThreshHi', help='Upper bound of value to describe region to calculate volume',
                      default=None)
  parser.add_option('-F', '--criteriaFile', action='store', type='string',
                      dest='criFile', help='File with variable to constrain calculation of volume',
                      default=None)
  parser.add_option('-V', '--criVar', action='store', type='string',
                      dest='criVar', help='Variable in criteria file to calculate volume',
                      default=None)
  parser.add_option('-L', '--criThreshLow', action='store', type='float',
                      dest='criThreshLow', help='Lower bound of criteria value to describe region to calculate volume',
                      default=None)
  parser.add_option('-U', '--criThreshHi', action='store', type='float',
                      dest='criThreshHi', help='Upper bound of criterita value to describe region to calculate volume',
                      default=None)
  parser.add_option('-t', '--regionRegionFile', action='store', type='string',
                      dest='regionRegionGR3File', help='a gr3 file indicating the region region (depth==1)',
                      default=None)
  parser.add_option('-n', '--regionName', action='store', type='string',
                      dest='regionName', help='Name of the region', default='region')

  (options, args) = parser.parse_args()
  dataDir       = options.dataDir
  startStr      = options.startStr
  endStr        = options.endStr
  runTag        = options.runTag
  volFile       = options.volFile
  volVar        = options.volVar
  volThreshHi   = options.volThreshHi
  volThreshLow  = options.volThreshLow
  criFile       = options.criFile
  criVar        = options.criVar
  criThreshHi   = options.criThreshHi
  criThreshLow  = options.criThreshLow
  regionName    = options.regionName
  regionRegionGR3File = options.regionRegionGR3File

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
  if not volFile :
    parser.print_help()
    parser.error('volVar undefined')
  if not volVar :
    parser.print_help()
    parser.error('volVar undefined')
  if criFile :
    if not criVar :
      parser.print_help()
      parser.error('criVar undefined')
    if criThreshHi is None and criThreshLow is None :
      parser.print_help()
      parser.error('criThreshLow or criThreshHi must be defined')

  startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
  endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')

  print 'Parsed options:'
  print ' - time range:',str(startTime),'->', str(endTime)
  print ' - dataDir',dataDir
  print ' - runTag',runTag
  print ' - volFile',volFile
  print ' - volVar',volVar
  print ' - volThreshLow',volThreshLow
  print ' - volThreshHi',volThreshHi
  print ' - criFile',criFile
  print ' - criVar',criVar
  print ' - criThreshLow',criThreshLow
  print ' - criThreshHi',criThreshHi
  print ' - region file',regionRegionGR3File
  print ' - region name',regionName

  # Extract
  vc = volumeComputer( dataDir, volFile, volVar, volThreshLow, volThreshHi,
                       criFile, criVar, criThreshLow, criThreshHi,
                       regionRegionGR3File, regionName)
  dcs = vc.getDataContainer( startTime, endTime )
  for dc in dcs :
    dc.setMetaData('tag',runTag)

  import data.dirTreeManager as dtm
  #rule = dtm.oldTreeRule()
  rule = dtm.defaultTreeRule()
  dtm.saveDataContainerInTree( dcs, rule=rule, dtype=np.float32,
                               overwrite=True )

if __name__=='__main__' :
  parseCommandLine()
