"""
A simple database implementation to handle all observation and model data for
skill assessments.

Tuomas Karna 2012-09-27
"""
#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import os
import numpy as np
import datetime
import sys
import traceback

from crane.data import dataContainer
from crane.data import timeArray
import crane.data.netcdfCacheInterface as netcdfDB
import crane.data.loadHindcastStations as loadHindcastStations
from crane.data import harmonicAnalysis
from crane.plotting.plotBase import createDirectory
from crane.data import timeSeriesFilters

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

EXTRACT_VARS = ['elev','temp','salt','hvel','turbidity','NO3','oxy','fluores']

#-------------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------

class unitConversion(object):
  """
  Converts units from observation unit to model units and vice versa.
  """
  def __init__(self, obsTag='obs'):
    """Initialize convert object. obsTag defines the tag used to identify
    observation data. Data with all other tags are assumed to be in model
    units (for now).
    """
    self.obsTag = obsTag
    # TODO generalize with unit conversion tool
    self.scalarToModelUnits = {'oxy':44.661, # from 1 ml/l to mmol/m3
                              }
    self.modelUnits = {'oxy':'mmol/m3',
                       }

  def convertToModelUnits(self, dc):
    """Converts given dataContainer from observation units to model units."""
    if dc.getMetaData('tag', suppressError=True) != self.obsTag:
      return dc
    var = dc.getMetaData('variable')
    dc2 = dc.copy()
    dc2.data *= self.scalarToModelUnits.get(var, 1.0)
    if var in self.modelUnits:
      dc2.setMetaData('unit', self.modelUnits[var])
    if self.scalarToModelUnits.get(var, 1.0) != 1.0:
        print '*** Converted'
        print dc
        print dc2
    return dc2

def filterTuples( allTuples, queries, exclude=False ) :
  """A generic routine for filtering lists of tuples.
  Returns all items in allTuples that match the queries.
  Query is a tuple indicating the desired attributes. For example,
  ('a',1,None) would match all cases where the first item is 'a' and the second is 1.
  None is used as a wildcard, so in this case all values for the last item are accepted.
  Setting exclude=False excludes the matching tuples from the list.
  """
  olist = set()
  if not exclude :
    for query in queries :
      for tup in allTuples :
        if tupleMatches(tup,query) :
          olist.add( tup )
  else :
    olist = set(allTuples)
    for query in queries :
      for tup in allTuples :
        if tupleMatches(tup,query) :
          olist = olist.difference([tup])
  return list(olist)

def tupleMatches( tup, query ) :
  """Compares a tuple to a query tuple. Query tuple has the same size as the tuple.
  None is used as a wildcard; Items with None are not tested and any value is accepted.
  """
  # find non-None keys for query
  nonZeroIx = [ i for i in range(len(query)) if query[i] ]
  # compare tuples, here's the one-liner!
  return all( v[0] == v[1] for i,v in enumerate(zip(tup,query)) if i in nonZeroIx )

def fetchAvailableObservations( startTime, endTime, obsTag='obs',
                                variables=['elev','temp','salt'] ) :
  """Fetches all available observations from CMOP's netCDF cache.
  Returns a StationCollection object.
  """
  offerings = netcdfDB.getAvailableOfferings( startTime, endTime, variables )
  if not offerings :
    raise Exception( 'No offerings could be retrieved' )
  sc = StationCollection( startTime, endTime, obsTag )
  conv = unitConversion(obsTag)
  for off in offerings :
    try :
      dc = netcdfDB.getDataContainerFromOffering(off, startTime, endTime)
      dc = conv.convertToModelUnits(dc)
      sc.addSample(dc)
    except Exception as e :
      #print 'Fetching data failed', off
      print e
  return sc

def fetchHindcastFromDatFiles( tag, hindcastStr, offerings, startTime, endTime,
  removeBadValues=False, stationFile=None ) :
  """Fetches model data from the hindcast database *.dat files.
  3D variables are interpolated to the correct depth.
  Returns a StationCollection object."""
  from files.stationFile import StationFile
  sc = StationCollection( startTime, endTime )
  baseDir =  '/home/workspace/ccalmr/hindcasts/'

  # group offerings based on station
  groupedOff = {}
  for off in offerings :
    sta = off['location']
    if off['variable'] not in EXTRACT_VARS :
      continue
    groupedOff.setdefault( sta , [] ).append( off )
  stations = groupedOff.keys()

  csvReader = csvStationFile()
  csvReader.readFromFile(stationFile)
  for sta in stations :
    try :
      x,y = csvReader.getLocation( sta )
      extractRequest = []
      for off in groupedOff[sta] :
        z = -float(off['msldepth'])/100.0
        extractRequest.append( (off['variable'],z,off['bracket']) )
      dcs = loadHindcastStations.readDateRange(baseDir, hindcastStr, sta, extractRequest, startTime, endTime, x, y, 'spcs', removeBadValues)
      for dc in dcs :
        dc.setMetaData( 'tag',tag )
        sc.addSample( dc )
    except Exception as e :
      print 'Extraction failed for station',sta
      print e
  return sc

def extractForOfferings( tag, dataDir, offerings, startTime, endTime,
                          modelCoordSys='spcs', stationFile=None, profile=False, netcdf=False ) :
  """Extracts station time series from the model output files stored in dataDir, based on the given StationCollection."""
  import data.extractStation as es
  sc = StationCollection( startTime, endTime )
  # sort offerings for each variable
  offForVar = dict()
  for off in offerings :
    var = off['variable']
    if var.split('.')[0] not in EXTRACT_VARS :
      continue
    if not var in offForVar :
      offForVar[var] = []
    offForVar[var].append( off )
  for var in offForVar :
    # extract
    try :
      if netcdf :
        import data.ncExtract as nce
        dataContainers = nce.extractForOfferings( dataDir, var, offForVar[var],
                                                  startTime, endTime,
                                                  stationFile=stationFile )
      else :
        dataContainers = es.extractForOfferings( dataDir, var, offForVar[var],
                                                startTime, endTime,
                                                profile=profile,
                                                modelCoordSys=modelCoordSys,
                                                stationFile=stationFile )
      # add to collection
      for dc in dataContainers :
        dc.setMetaData( 'tag', tag )
        sc.addSample( dc )
    except Exception as e :
      print 'Extraction failed for variable',var
      traceback.print_exc(file=sys.stdout)      
  return sc

def extractTrackFromModelOutputs( tag, dataDir, sColl, startTime, endTime,
                          modelCoordSys='spcs') :
  """Extracts track data from model outputs stored in dataDir, based on the track samples available in sColl (other stationCollection)."""
  import data.extractTrack as et
  import data.extractStation as es
  sc = StationCollection( startTime, endTime )
  # extract track data
  for key in sColl.getKeys() :
    trackDC = sColl.getSample(**key)
    if not trackDC.isTrack() :
      continue
    # TODO FIXME
    #dc = et.extractForDataContainer( dataDir, trackDC )
    var = key[-1]
    offerings = [trackDC.description]
    dcs = es.extractForOfferings( dataDir, var, offerings, startTime, endTime, profile=True )
    dc = dcs[0]
    dc.setMetaData( 'tag',tag )
    sc.addSample( dc )
  return sc

def extractProfiles( tag, dataDir, varList, startTime, endTime,
                     stationFile, stationNames=[], modelCoordSys='spcs', netcdf=False) :
  """Extracts vertical profiles for all stations in the stationFile for the given variable."""
  import data.extractStation as es
  from files.csvStationFile import csvStationFile
  sc = StationCollection( startTime, endTime )

  # all stations
  csvReader = csvStationFile()
  csvReader.readFromFile(stationFile)

  if len(stationNames)==0 :
    stationNames = list( set(csvReader.getStations()) )
  else :
    # remove stations that are missing in stationFile
    stationNames = [ s for s in stationNames if s in csvReader.getStations() ]
  if len(stationNames)==0 :
    print 'warning: list of stations is empty'
    return
  x = np.array( [ csvReader.getX(s) for s in stationNames ] )
  y = np.array( [ csvReader.getY(s) for s in stationNames ] )

  print ' *** extracting profiles for stations *** '
  for s in stationNames :
    print s

  # execute
  for var in varList :
    try :
      print ' *',var
      if netcdf :
        import data.ncExtract as nce
        pe = nce.selfeExtract(dataDir,var=var)
        dataContainers = pe.extractVerticalProfile(startTime,endTime,
                                                   var,x,y,stationNames)
      else :
        ee = es.extractStation(dataDir,var,profile=True,modelCoordSys=modelCoordSys)
        ee.setStations( stationNames, x,y )
        dataContainers = ee.extractDates( startTime, endTime )
      # add to collection
      for dc in dataContainers :
        dc.setMetaData( 'tag', tag )
        sc.addSample( dc )
    except Exception as e :
      print 'Extraction failed for variable',var
      print e

  return sc

def computeDerivedProducts( sColl,tagsToProcess ) :
  """Calculates all derived products that depend on observation/extracted data.
  Returns the derived products in a new stationCollection."""
  sc = StationCollection( sColl.startTime, sColl.endTime )
  dataContainers = []
  
  # Saturn01 profiler data
  for tag in tagsToProcess :
    for var in ['salt','temp'] : # TODO get the vars based on obs?
      dc = sColl.getSaturn01ProfilerData( tag, var, regenerate=True )
      if dc :
        dataContainers.append( dc )

  dcs = sColl.generateStratificationData( tag=tagsToProcess, variable='salt' )
  dataContainers.extend( dcs )
  
  # add to collection
  for dc in dataContainers :
    sc.addSample( dc )
  return sc


class tinyDB( object ) :
  """A implementation of a simple dictionary-based database.

  Data associated with a tuple of keywords is stored in a dictionary.
  A list of relevant keywords are given to the constructor. Samples with fever
  keywords can be added in the database; missing keywords are given value None.
  The class supports basic filtering operations as well as combinations
  (union, difference etc.)."""
  def __init__( self, keywords, source=None ) :
    """Create an empty database using the given keywords."""
    if not source :
      self.keywords = keywords
      self.data = {}
    else :
      self.keywords = list( source.keywords )
      self.data = dict( source.data )

  def __getitem__( self, keys ) :
    return self.data[keys]
  def __iter__( self ) :
    return self.data.__iter__()
  def __contains__( self,item ) :
    return self.data.__contains__(item)
  def __len__(self) :
    return self.data.__len__()
  def keys( self ) :
    return self.data.keys()
  def has_key( self, key ) :
    return self.data.has_key(key)

  def genKey( self, **kwargs ) :
    """Given the kwargs, generates a tuple of keywords.
    Missing keywords are replaced by None."""
    return tuple( [ kwargs.get(k,None) for k in self.keywords ] )

  def addSample( self, sample, **keywords ) :
    """Adds sample to the database, for the given keywords.
    If a sample already exists for the given keywords, it is silently replaced."""
    key = self.genKey( **keywords )
    self.data[key] = sample

  def getTuples( self, query=None, exclude=False, **kwargs ) :
    """Returns a list of all keys that match the query.
    query is a list of dictionaries, with appropriate keywords.
    Alternatively, keywords can be given as argument to this function.
    In this case one keyword can be a list.
    Setting exclude=True negates the query.
    """
    if query == None :
      # kwargs is a single dict, one value may be list
      nListArgs = sum( isinstance(v,list) for v in kwargs.values() )
      if nListArgs > 1 :
        raise Exception( 'Only one of kwargs can be a list')
      if nListArgs :
        # construct tuples
        for k in kwargs.keys() :
          if isinstance(kwargs[k],list) :
            break
        tuples = [ ]
        for i in range(len(kwargs[k])) :
          # construct dicts with one item of the list in k
          d = dict(kwargs)
          d[k] = kwargs[k][i]
          tuples.append( self.genKey( **d ) )
      else :
        tuples = [ self.genKey( **kwargs ) ]
    else :
      # filter is a list of dict
      # convert to list of tuples
      tuples = [ self.genKey( **d ) for d in query ]
    #print tuples
    return filterTuples( self.data.keys(), tuples, exclude )

  def getKeys( self, query=None, exclude=False, **kwargs ) :
    tuples = self.getTuples( query,exclude, **kwargs )
    return [ dict( (self.keywords[j],tuples[i][j]) for j in range(len(self.keywords)) ) for i in range(len(tuples))]

  def getAttributes( self, keyword, query=None, exclude=False, **kwargs ) :
    """Returns all occurences of a keyword.
    Example:
    sc.setAttributes( 'location',variable='salt')
    returns all locations that have salinity data, e.g. ['saturn01','jetta']
    """
    allVals = list( set( k[keyword] for k in self.getKeys(query, exclude, **kwargs) ) )
    return allVals
  
  def getSample( self, **kwargs ) :
    """Returns a sample corresponding to the given keywords.
    Returns None if no sample is found."""
    for k in kwargs :
      if k not in self.keywords :
        raise Exception('Unknown keyword argument: '+k)
    keys = self.getKeys( **kwargs )
    if keys :
      key = keys[0]
      k = self.genKey( **key )
      return self.data[ k ]
    else :
      return None

  def getSubset( self, query=None, exclude=False, **kwargs ) :
    subset = tinyDB( self.keywords )
    for k in self.getKeys(query, exclude, **kwargs ) :
      subset.addSample( self.getSample( **k ), **k )
    return subset

  def update( self, other ) :
    """Add samples of other to self."""
    for k in other.getKeys() :
      self.addSample( other.getSample(**k), **k )

  def union( self, *others ) :
    keywords = set(self.keywords)
    for o in others :
      keywords.update( o.keywords )
    keywords = list( keywords )
    union = tinyDB( keywords )
    for k in self.getKeys() :
      union.addSample( self.getSample(**k), **k )
    for o in others :
      for k in o.getKeys() :
        union.addSample( o.getSample(**k), **k )
    return union

class StationCollection(tinyDB) :
  """A simple "database" container for station time series data."""
  def __init__( self, startTime, endTime, obsTag=None, data=None,
                keywords=['dataType', 'tag','location','msldepth','variable'] ) :
    if not obsTag :
      obsTag = 'obs'
    self.obsTag = obsTag
    self.startTime = startTime
    self.endTime = endTime
    # condition samples TODO better implementation?
    self.haManager = HAManager( self )
    tinyDB.__init__(self,keywords,source=data)

  def getKeywordsFromDC( self, dc ) :
    """Returns keywords in a dict as for the given dc dataContainer."""
    keywords = dc.getMetaData()
    for k in self.keywords :
      if k not in ['msldepth'] and k not in keywords :
        print ' * warning: keyword '+k+' missing', keywords
    return dict( (k,keywords[k]) for k in self.keywords if k in keywords )

  def getOfferings( self ) :
    """Returns offerings for all timeseries data in the collection."""
    offerings = []
    for k in self.getKeys( dataType='timeseries' ) :
      dc = self.getSample( **k )
      off = {}
      for i in ['location','msldepth','bracket','instrument','variable'] :
        off[i] = dc.getMetaData(i)
      offerings.append( off )
    return offerings

  def addSample( self, dc, **kwargs ) :
    """Adds dataContainer in the structure."""
    try :
      truncated = dc.timeWindow( self.startTime, self.endTime, includeEnd=True )
      if ( truncated.getMetaData('dataType') in ['timeseries','profile'] and
           len(truncated.time) < 3 ) :
        # reject too short time series
        m = truncated.getMetaData()
        print 'time series too short',m['location'],m['variable'],m['msldepth']
        return
      if len(kwargs) == 0 :
        # try to get keywords from dataContainer (should be default usage)
        kwargs = self.getKeywordsFromDC( dc )
      tinyDB.addSample( self, truncated, **kwargs )
    except Exception as e :
      print e

  def getModelTags(self) :
    """Returns all model tags in the collection."""
    modTags = self.getAttributes( 'tag' )
    if self.getObsTag() in modTags :
      modTags.remove(self.getObsTag())
    return modTags

  def getObsTag(self) :
    """Returns the observation tag in the collection."""
    if self.obsTag == None :
      raise Exception( 'obs tag has not been set' )
    return self.obsTag

  def getComparableKeys( self, query=None, exclude=False,
                         requireObs=True, requireMod=True, **kwargs ) :
    """Returns a list of observation, model tuples of the database.
    Each entry is [ obsKeys, [modKey1,modKey2,...] ], where obsKey and modKey
    are the keys for observation and model data, respectively.
    See getKeys for filtering arguments."""
    # matching keys
    keys = self.getKeys(query, exclude, **kwargs)
    # store in a dict
    models = dict()
    observ = dict()
    for k in keys :
      sta = k['location']
      var = k['variable']
      dep = k['msldepth']
      tag = k['tag']
      if tag == self.getObsTag() :
        observ[(sta,var,dep)] = k
      else :
        models.setdefault( (sta,var,dep), [] ).append( k )
    # sort model keys for consistency
    for k in models :
      models[k] = sorted(models[k])
    # construct obs,mod pairs
    omPairs = []
    if requireObs :
      for k in observ :
        if not requireMod or k in models :
          entry = [ k, observ[k], models[k] ]
          omPairs.append( entry )
    else :
      if requireMod :
        allKeys = models.keys()
      else :
        allKeys = set(models.keys()).union(set(observ.keys()))
      for k in allKeys :
        entry = [ k, observ.get(k,None), models.get(k,[]) ]
        omPairs.append( entry )
    return omPairs

  def getSubset( self, query=None, exclude=False, **kwargs ) :
    """Returns a subset of current StationCollection defined by the query list.
    See getKeys for arguments.
    """
    data = tinyDB.getSubset(self,query,exclude,**kwargs)
    return StationCollection( self.startTime, self.endTime, self.obsTag, data=data )
  
  def copy( self ):
    """Returns a copy of this StationCollection."""
    return self.getSubset()

  def union( self, *others ) :
    """Return a new StationCollection with elements from the current object and all others.
    """
    data = tinyDB.union( self, *others )
    e = StationCollection(self.startTime, self.endTime, self.obsTag, data=data)
    return e

  def intersection( self,  *others ) :
    """Return a new StationCollection with elements common to the current object and all others.
    """
    raise NotImplementedError("This method has not yet been implemented")

  def difference( self,  *others ) :
    """Return a new StationCollection with elements in the current object that are not in the others.
    """
    raise NotImplementedError("This method has not yet been implemented")
  
  def saveAsNetCDFCollection( self, dataDir='data', treeType='default' ) :
    """Saves each time series in a netCDF file. The destination is
    dataDir/tag/station_variable_msldepth_startTime_endTime.nc
    All necessary directories are created as needed
    """
    from data.dirTreeManager import netcdfTree, oldTreeRule, defaultTreeRule
    for k in self.getKeys() :
      dc = self.getSample( **k )
      tag = dc.getMetaData('tag')
      sta = dc.getMetaData('location')
      var = dc.getMetaData('variable')
      dep = dc.getMetaData('msldepth',suppressError=True)
      bra = dc.getMetaData('bracket',suppressError=True)
      dat = dc.getMetaData('dataType')

      rule = oldTreeRule(dataDir=dataDir) if treeType=='old' else defaultTreeRule(dataDir=dataDir)
      tree = netcdfTree( dataType=dat, tag=tag, location=sta, variable=var, msldepth=dep, bracket=bra, rule=rule )
      tree.write( dc, overwrite=True )

  @classmethod
  def loadFromNetCDFCollection( cls, tag, startTime=None, endTime=None,
                                obsTag=None, dataDir='data', treeRule=None, dataType=None, variable=None, verbose=True ) :
    from data.dirTreeManager import netcdfTreeTraverser

    tra = netcdfTreeTraverser( rule=treeRule, verbose=verbose )
    dcs,st,et = tra.readFiles( tag=tag, dataType=dataType, variable=variable,
                         startTime=startTime, endTime=endTime, msldepth=None )
    if not startTime :
      startTime = st
    if not endTime :
      endTime = et
    sc = cls( startTime, endTime, obsTag )
    # add to collection
    for dc in dcs :
      sc.addSample( dc )
    return sc

  def fetchDischargeData( self, offering=None ) :
    """Fetches discharge data from the database and stores it in the collection."""
    # get discharge data
    if not offering :
      offering = {'location':'bono3','msldepth':'0','bracket':'A',
                  'instrument':'FLUX','variable':'flux'}
    try :
      dc = netcdfDB.getDataContainerFromOffering( offering,
                                                  self.startTime, self.endTime )
    except Exception as e :
      #print 'Fetching data failed', off
      print e
      return
    self.addSample( dc )

  def getDischargeData( self, locs=['bono3','bvao3'] ) :
    """Returns river discharge time series for the time period in question"""
    keys = self.getKeys( variable='flux' )
    if len(keys) == 0 :
      self.fetchDischargeData()
      keys = self.getKeys( variable='flux' )
    if not keys :
      return None
    else :
      # find key with desired location(s)
      key = keys[0] # default
      found=False
      if isinstance(locs,str) : locs=[locs]
      for loc in locs :
        for k in keys :
          if k['location']==loc :
            found=True
            key=k
            break
        if found : break
      return self.getSample( **key )

  def fetchTidalData( self, offering=None, dc=None ) :
    """Fetches water level data, computes tidal range, smooths it,
    and stores it in the collection.
    """
    if not offering :
      offering = {'location':'tpoin','msldepth':'0','bracket':'A',
                  'instrument':'External','variable':'elev'}
    # get elevation
    T = 44714. # M2 period in seconds
    if not dc :
      pad = datetime.timedelta(seconds=3*T)
      dc = netcdfDB.getDataContainerFromOffering( offering,
                                self.startTime-pad, self.endTime+pad )
    # find tidal range
    import time
    t0 = time.clock()
    gaps,ranges,t = dc.detectGaps()
    x = np.squeeze(dc.data)
    tRes = np.array([])
    xRes = np.array([])
    for i in range( ranges.shape[0] ) :
      tt = t[ ranges[i,0]:ranges[i,1] ]
      xx = x[ ranges[i,0]:ranges[i,1] ]
      tt, xx = timeSeriesFilters.computeRunningRange( tt,xx,2*T )
      if len(tt) == 0 :
        continue
      tt, xx = timeSeriesFilters.computeRunningMean( tt,xx,3*T )
      tRes = np.append( tRes, tt )
      xRes = np.append( xRes, xx )
    print 'tidal range cpu time',time.clock() - t0
    if len(tRes) == 0 :
      print 'tidal data could not be computed, skipping (time series too short?)'
      return
    ta = timeArray.timeArray(tRes, 'epoch' )
    data = xRes.reshape( (1,1,-1) )
    meta = dc.getMetaData()
    meta['dataType'] = 'timeseries'
    meta['tag'] = self.getObsTag()
    meta['variable'] = 'tidal_range'
    dc2 = dataContainer( '', ta, dc.x,dc.y,dc.z, data,
                              ['tidal_range'], coordSys='',metaData=meta)

    self.addSample( dc2 )

  def getTidalData( self ) :
    """Returns tidal range time series for the time period in question"""
    keys = self.getKeys( variable='tidal_range' )
    if len(keys) == 0 :
      self.fetchTidalData()
      keys = self.getKeys( variable='tidal_range' )
    if not keys :
      return None
    else :
      k = keys[0]
      return self.getSample( **k )

  def fetchCoastalUpwellingData( self, filename=None ) :
    """Fetches costal upwelling index from default file and stores it in the collection."""
    from files.cuiFileParser import cuiParser
    try :
      dc = cuiParser(filename).getDataContainer( self.startTime, self.endTime )
    except Exception as e :
      print 'CUI data could not be retrieved for period :',self.startTime, self.endTime
      print e
      return
    self.addSample( dc )
    
  def getCoastalUpwellingData( self ) :
    """Returns coastal upwelling index time series for the time period in question"""
    keys = self.getKeys( variable='cui' )
    if len(keys) == 0 :
      self.fetchCoastalUpwellingData()
      keys = self.getKeys( variable='cui' )
    if not keys :
      return None
    else :
      k = keys[0]
      return self.getSample( **k )

  def generateSaturn01ProfilerData( self, tag, var ) :
    """Uses saturn01 pressure data, and temp/salt data to generate a temp/salt dataContainer with correct depth information in z coordinate."""
    import data.wprofiler as wp
    import time as timeMod
    if tag == self.getObsTag() :
      # TODO what is the correct dataType for raw profiler data?
      depDC = self.getSample( tag=tag, location='saturn01',
                              variable='elev', dataType='profiler_raw' )
      varDC = self.getSample( tag=tag, location='saturn01',
                              variable=var, dataType='profiler_raw' )
      if not depDC or not varDC :
        print 'Could not find saturn01 profiler obs data',var, varDC==None, depDC==None
        return None
      cputime0 = time.clock()
      sys.stdout.write( 'calculating saturn01 profiler obs data ...' )
      sys.stdout.flush()
      try :
        dc = wp.generateSat01ProfilerObsData( depDC, varDC )
        sys.stdout.write( ' duration %.2f s\n'%(timeMod.clock()-cputime0) )
      except Exception as e :
        print 'Could not generate saturn01 profiler obs data'
        print e
        dc = None
    else : # model
      # reference obs
      obsDC = self.getSaturn01ProfilerData( self.getObsTag(), var )
      vprofDC = self.getSample( tag=tag, location='saturn01',
                              variable=var, dataType='profile' )
      if not obsDC or not vprofDC :
        rStr = ''
        if not obsDC : rStr+=' obs missing,'
        if not vprofDC : rStr+=' model profile missing,'
        print 'Could not generate saturn01 profiler model data:', rStr, var
        return None
      try :
        dc = wp.generateSat01ProfilerModData( obsDC, vprofDC )
      except Exception as e :
        print 'Could not generate saturn01 profiler model data'
        print e
        dc = None
    return dc

  def getSaturn01ProfilerData( self, tag, var, regenerate=False ) :
    """Returns saturn01 profiler data correctly annotated with depth,
    if availailable."""
    dc = self.getSample( tag=tag, variable=var,
                         location='saturn01', dataType='profiler' )
    if not dc or regenerate :
      dc = self.generateSaturn01ProfilerData( tag, var )
      if dc :
        # add to collection
        dc.setMetaData( 'tag', tag )
        self.addSample( dc )
    return dc

  def generateStratificationData( self, **kwargs ) :
    """Computes stratification for all stations with profile data."""

    # get all stations with profiler data
    kwargs['dataType'] = 'profile'
    kwargs['variable'] = 'salt'
    keys = self.getKeys( **kwargs )
    dcs = []
    for k in keys :
      s = self.getSample( **k )
      # compute stratification
      ixTop = np.argmax( s.z[:,0] )
      ixBot = np.argmin( s.z[:,0] )
      strat = s.data[ixBot,0,:]-s.data[ixTop,0,:]
      # create new dataContainer
      data = strat.reshape( (1,1,-1) )
      z = 0
      meta = s.getMetaData()
      #meta.pop('msldepth',None) # TODO
      meta['msldepth'] = '0'
      meta['dataType'] = 'timeseries'
      meta['variable'] = 'strat'
      dc = dataContainer( '', s.time, s.x[0],s.y[0],z, data,
                                ['strat'], coordSys=s.coordSys, metaData=meta )
      self.addSample( dc )
      dcs.append(dc)
    return dcs

  def getStratificationData( self, **keywords ) :
    """Returns stratification data for the requested tag and station"""
    dc = self.getSample( **keywords )
    if not dc or regenerate :
      dc = self.generateStratificationData( **keywords )[0]
    return dc

  def getHarmonicAnalysis( self, key, obsKey=None ) :
    return self.haManager.get( key, obsKey )

class HAManager( object ) :
  """A helper class to compute/save/load harmonic analysis results."""
  def __init__( self, stationColl, subdirectory='data/ha' ) :
    self.coll = stationColl
    self.subDir = subdirectory

  def generateFilename( self, dc ) :
    fname = '_'.join( [dc.getMetaData('location'),
                       dc.getMetaData('variable'),
                       dc.getMetaData('msldepth')  ] )
    st = str(dc.time.getDatetime( 0).date())
    et = str(dc.time.getDatetime(-1).date())
    fname = '_'.join(['ha',fname,st,et])+'.txt'
    return fname

  def save( self, key, dc, tc ) :
    """Stores tidal constituent data to disk"""
    fname = self.generateFilename( dc )
    path = os.path.join( key['tag'], self.subDir )
    createDirectory(path)
    tc.saveToASCII( fname, path=path )

  def compute( self, signal ) :
    t0 = time.clock()
    sys.stdout.write( 'computing harmonic analysis: %s... '%(signal.description) )
    sys.stdout.flush()
    tc =  harmonicAnalysis.tidalConstituents.computeFromData( signal )
    sys.stdout.write( '%.2f s\n'%(time.clock()-t0) )
    return tc

  def get( self, key, obsKey=None, useErrorSignal=False ) :
    """Returns tidal constituent data.
    If it has been computed previously, reads the file from disk, otherwise
    computes constituents and saves to disk. If both key and obsKey are given,
    operates on the error signal. If useErrorSignal is set, computes HA of
    the error signal, otherwise computes the difference of amplitudes and
    phases directly."""
    if obsKey and not useErrorSignal :
      # compute the difference of constituents
      tcmod = self.get( key )
      tcobs = self.get( obsKey )
      if not tcmod or not tcobs :
        return None
      commonConsts = list( set(tcmod.constituents).intersection( set(tcobs.constituents) ) )
      ampDict = {}
      phaDict = {}
      for const in commonConsts :
        # compute difference
        ampDict[const] = tcmod.getConstituent(const)[0] - tcobs.getConstituent(const)[0]
        phaDict[const] = tcmod.getConstituent(const)[1] - tcobs.getConstituent(const)[1]
      # fast to compute, no need to save
      tc = harmonicAnalysis.tidalConstituents( 'error_'+tcmod.description, ampDict, phaDict, tcmod.fieldName, tcmod.startTime, tcmod.endTime )
      return tc

    signal = self.coll.getSample( **key )
    if obsKey : # compute error
      try :
        o = self.coll.getSample( obsKey )
        signal = o.computeError( signal )
      except Exception as e:
        print 'Error could not be computed, skipping:', key
        print e
        return None
    fname = self.generateFilename( signal )
    fname = os.path.join(key['tag'],self.subDir,fname)
    tc = None
    if os.path.isfile( fname ) :
      tc = harmonicAnalysis.tidalConstituents.loadFromASCII( fname )
    else :
      try :
        tc = self.compute( signal )
        self.save( key, signal, tc )
      except Exception as e:
        print 'Harmonic analysis failed, skipping'
        print e
    return tc
