"""
Generic timeArray class.

Supported time formats:
epoch -- seconds since 1970-01-01 00:00:00 UTC
corie -- days since 1995-12-31 00:00:00 PST
simulation -- seconds since simulation started
datetime object -- always treated as PST

System time zone information is discarded.
Daylight saving time is not dealt with.

NOTE:
By default, time.mktime converts from system's local time zone to UTC epoch.
By default, datetime.fromtimestamp converts UTC epoch to system's time zone.
System time zone is given by time.timezone
Here the system's time zone is ingnored and PST is used instead.

Tuomas Karna 2012-08-16
"""

import time
import datetime
from calendar import timegm
import numpy as np
pstTZOffset = np.float64(28800.0)
# corie zero time as datetime object
corieZeroDate = datetime.datetime(1995,12,31,0,0,0)
# corie zero time as epoch time stamp
corieZeroEpoch = np.float64(timegm(corieZeroDate.timetuple())) + pstTZOffset

class timeArray(object) :
  """Generic time array class.
  timeFormat can be:
  'epoch' -- seconds since 1970-01-01 00:00:00 UTC
  'corie' -- days since 1995-12-31 00:00:00 PST
  'simulation' -- seconds since simulation started, startDate must be provided
  """
  supportedFormats = [ 'corie', 'epoch', 'simulation' ]
  dtype = float
  def __init__(self,array,timeFormat,startDate=None,acceptDuplicates=False) :
    """Create a new timeArray instance.
    If timeFormat=='simulation', startDate is required (as datetime object).
    Otherwise startDate is ignored.
    
    Time array inherits basic array functionality:
    timeArray[0], timeArray[:10], 14.0 in timeArray
    are all valid operations.
    """
    if len(array) == 0 :
      raise Exception( 'Given time array is empty' )
    if timeFormat not in timeArray.supportedFormats :
      raise Exception( 'Given time format not implemented: '+timeFormat )
    if acceptDuplicates :
      isIncreasing = np.all(np.diff(array) >= 0)
    else :
      isIncreasing = np.all(np.diff(array) > 0)
    if not isIncreasing :
      print array
      ix = np.nonzero(np.diff(array) <= 0)[0]
      print array[max(0,ix[0]-3):min(len(array),ix[0]+3)]
      raise Exception( 'Time sequence not monotonically increasing' )
    self.array = array.astype(self.dtype)
    self.timeFormat = timeFormat
    if timeFormat == 'simulation' :
      self.startDate = startDate
  
  # inherit basic array functionality
  def __len__(self) :
    return self.array.__len__()
  def __getitem__(self, key) :
    return self.array.__getitem__(key)
  def __setitem__(self, key, value) :
    return self.array.__setitem__(key, value)
  def __iter__(self) :
    return self.array.__iter__()
  def __contains__(self, item) :
    return self.array.__contains__(item)
  def max(self) :
    return max(self)
  def min(self) :
    return min(self)
  def __lt__(self, other) :
    """true if time.max() < other.min() in epoch format"""
    return self.asEpoch().max() <  other.asEpoch().min()
  def __le__(self, other) :
    """true if time.max() <= other.min() in epoch format"""
    return self.asEpoch().max() <= other.asEpoch().min()
  def __eq__(self, other) :
    """true if all values are the same in epoch format"""
    return np.array_equal( self.asEpoch().array,  other.asEpoch().array )
  def __ne__(self, other) :
    """true if any value is different in epoch format"""
    return not self.__eq__(other)
  def __gt__(self, other) :
    """true if time.min() > other.max() in epoch format"""
    return self.asEpoch().min() >  other.asEpoch().max()
  def __ge__(self, other) :
    """true if time.min() >= other.max() in epoch format"""
    return self.asEpoch().min() >=  other.asEpoch().max()
    
  def asEpoch(self) :
    """Return an equivalent timeArray instance in epoch time"""
    array = self.array
    if self.timeFormat == 'corie' :
      array = corieToEpochTime(self.array)
    elif self.timeFormat == 'simulation' :
      array = simulationToEpochTime(self.array,self.startDate)
    return timeArray( array, 'epoch', acceptDuplicates=True )
    
  def asCorie(self) :
    """Return an equivalent timeArray instance in corie time"""
    array = self.array
    if self.timeFormat == 'epoch' :
      array = epochToCorieTime(self.array)
    elif self.timeFormat == 'simulation' :
      tmparray = simulationToEpochTime(self.array,self.startDate)
      array = epochToCorieTime(tmparray)
    return timeArray( array, 'corie', acceptDuplicates=True )

  def asSimulation(self, startDate) :
    """Return an equivalent timeArray instance in simulation time.
    Simulation startDate must be provided (as datetime object)."""
    array = self.array
    if self.timeFormat == 'epoch' :
      array = epochToSimulationTime(self.array, startDate)
    elif self.timeFormat == 'corie' :
      tmparray = corieToEpochTime(self.array)
      array = epochToSimulationTime(tmparray, startDate)
    elif self.timeFormat == 'simulation' and self.startDate != startDate :
      # translate start time
      tmparray = simulationToEpochTime(self.array,self.startDate)
      array = epochToSimulationTime(tmparray, startDate)
    return timeArray( array, 'simulation', startDate=startDate, acceptDuplicates=True )
  
  def toFormat(self, timeFormat, startDate = None ) :
    """Return an equivalent timeArray instance in given format.
    timeFormat can be another timeArray instance or a string.
    startDate is needed for 'simulation' format (as datetime object)."""
    
    if isinstance(timeFormat, type(self)) : # timeArray instance
      ta = timeFormat
      if ta.timeFormat == 'simulation' :
        startDate = ta.startDate
      timeFormat = ta.timeFormat
    if timeFormat == 'epoch' :
      return self.asEpoch()
    elif timeFormat == 'corie' :
      return self.asCorie()
    elif timeFormat == 'simulation' :
      if not startDate :
        raise Exception( 'startDate must be provided' )
      return self.asSimulation(startDate)
    else :
      raise Exception( 'Unknown time format: '+timeFormat )

  def getDatetime(self, index) :
    """Returns datetime object for the given index in the array."""
    t = self.array[index]
    if self.timeFormat == 'corie' :
      t = corieToEpochTime(t)
    elif self.timeFormat == 'simulation' :
      t = simulationToEpochTime(t, self.startDate )
    return epochToDatetime( t )
    
  def merge(self, otherArray, acceptDuplicates=False) :
    """Appends otherArray to the end of current array.
    otherArray is first converted to the same timeFormat."""
    otherArray = otherArray.toFormat(self)
    if acceptDuplicates:
      cannotAppend = self.max() > otherArray.min()
    else:
      cannotAppend = self.max() >= otherArray.min()
    if cannotAppend:
      print self
      print otherArray
      raise Exception( 'Time ranges cannot be concatenated' )
    self.array = np.hstack( ( self.array, otherArray.array ) )
  
  def covers( self, date ) :
    """Tests whether this array covers given time instance, i.e. array[0] <= date <= array[-1].
    date is datetime object"""
    return self.getDatetime(0) <= date and date <= self.getDatetime(-1)

  def overlaps( self, startTime, endTime ) :
    """Tests whether this array overlaps with the given time range.
    startTime,endTime are datetime objects"""
    #return self.getDatetime(0) <= date and date <= self.getDatetime(-1)
    return startTime <= self.getDatetime(-1) and endTime >= self.getDatetime(0)

  def getRangeIndices( self, startTime, endTime, includeEnd=False ) :
    """Returns indices of time stamps that belong to the given interval (end points inclusive). startTime, endTime are datetime instances."""
    if startTime == None and endTime == None :
      return range(len(self.array))
    if startTime == None :
      startTime = self.getDatetime(0)
    if endTime == None :
      endTime = self.getDatetime(-1)
      includeEnd = True
    if self.timeFormat == 'epoch' :
      startTime = datetimeToEpochTime( startTime )
      endTime = datetimeToEpochTime( endTime )
    elif self.timeFormat == 'corie' :
      startTime = datetimeToCorieTime( startTime )
      endTime = datetimeToCorieTime( endTime )
    elif self.timeFormat == 'simulation' :
      startTime = epochToSimulationTime(datetimeToEpochTime(startTime), self.startDate)
      endTime = epochToSimulationTime(datetimeToEpochTime(endTime), self.startDate)
    if includeEnd :
      ix = np.logical_and( self.array >= startTime, self.array <= endTime )
    else :
      ix = np.logical_and( self.array >= startTime, self.array < endTime )
    return np.nonzero( ix )[0] 

  def window( self, startTime, endTime, includeEnd=False ) :
    """Returns a copy of the current array, restricted on the given interval."""
    ix = self.getRangeIndices(startTime,endTime,includeEnd)
    arr = self.array[ ix ].copy()
    sd = None if self.timeFormat != 'simulation' else self.startDate
    return timeArray( arr, self.timeFormat, sd, acceptDuplicates=True )

  def __str__(self) :
    """Return a string that summarises timeArray content"""
    outStr = 'timeArray, {0:d} values\n format: {1:s}\n start: {2:s}\n end  : {3:s}'
    staStr = '{0:.4f} => {1:s} '.format( float(self.array[0]), str(self.getDatetime(0))  )
    endStr = '{0:.4f} => {1:s} '.format( float(self.array[-1]), str(self.getDatetime(-1)) )
    return outStr.format( len(self.array), self.timeFormat,staStr , endStr )
    
  def copy(self) :
    """Deep copy, all numpy arrays are copied"""
    sd = None if self.timeFormat != 'simulation' else self.startDate
    return timeArray( self.array.copy(), self.timeFormat, sd, acceptDuplicates=True )

  def detectGaps( self, dt=None, gapFactor=5 ) :
    """Detects gaps in the time series.

    Args:
    dt        -- (float) data sampling period. If None, taken as a mean step between data points.
    gapFactor -- (float) A factor to determine a minimum gap: gapFactor*dt

    Returns:
    gaps      -- (array) Indices that mark the beginning of each gap. (ngaps,)
    ranges    -- (array) Start and end indices of each contiguous data block. (ngaps+1,2)
    t         -- (array) time array in epoch format. (ntime,)
    """
    t = self.asEpoch().array
    steps = np.diff( t )
    if not dt :
      dt = np.mean(steps)
    gaps = np.nonzero( steps >= gapFactor*dt )[0]
    # convert to data ranges
    ranges = np.zeros( (len(gaps)+1,2),dtype=int )
    ranges[1:,0] = gaps + 1  # start indices
    ranges[:-1,1] = gaps # end indices
    ranges[-1,1] = len(self.array)-1
    return gaps, ranges, t

  def getAlignedTimeIndices(self,other,selfDt=None,otherDt=None) :
    """Returns time indices in other that fit in mutually overlapping time 
    ranges
    """
    # detect gaps
    sGaps,sRanges, st = self.detectGaps(dt=selfDt)
    oGaps,oRanges, ot = other.detectGaps(dt=otherDt)
    # merge ranges
    oTimeStamps = []
    for sIx in range( sRanges.shape[0] ) :
      sRangeStart = st[ sRanges[sIx,0] ]
      sRangeEnd = st[ sRanges[sIx,1] ]
      # check if any oRange fits in sRange
      for oIx in range( oRanges.shape[0] ) :
        oRangeStart = ot[ oRanges[oIx,0] ]
        oRangeEnd = ot[ oRanges[oIx,1] ]
        # test overlap
        if oRangeStart < sRangeEnd and oRangeEnd > sRangeStart :
          start = max( sRangeStart, oRangeStart )
          end = min( sRangeEnd, oRangeEnd )
          # get mathcing time stamps of other
          rangeIx = np.nonzero( np.logical_and( ot >= start, ot <= end ) )[0]
          oTimeStamps.extend( rangeIx )
    if not oTimeStamps :
      raise Exception( 'No overlapping intervals found in the given time series' )
    return oTimeStamps

# ------ helper functions ------
    
def generateSimulationTimeArray(startTime,endTime,dt) :
  """
  Generates a timeArray object for simulation time steps.
  
  startTime, endTime -- datetime objects
  dt                 -- time step in seconds
  returns:
  tsim               -- timeArray object
  """
  simDuration = datetimeToEpochTime(endTime)-datetimeToEpochTime(startTime)
  timeArr = np.arange(0,simDuration+dt,dt) # simulation time in sec
  time = timeArray( timeArr, 'simulation', startTime )
  return time
  
def corieToEpochTime(t) :
  """Convert corie to epoch time. t can be an array."""
  return np.float64(t)*86400.0 + corieZeroEpoch

def epochToCorieTime(t) :
  """Convert epoch to corie time. t can be an array."""
  return (t - corieZeroEpoch)/86400.0

def datetimeToEpochTime(t) :
  """Convert python datetime object to epoch time stamp.
  By default, time.mktime converts from system's local time zone to UTC epoch.
  Here the input is treated as PST and local time zone information is discarded.
  """
  return np.float64(timegm(t.timetuple())+t.microsecond*1e-6) + pstTZOffset

def epochToDatetime(t) :
  """Convert epoch time to python datetime object.
  By default, datetime.fromtimestamp converts UTC epoch to system's time zone.
  Here the local timezone information is discarded so that the output is always in PST.
  """
  return datetime.datetime.utcfromtimestamp(t) - datetime.timedelta(seconds=float(pstTZOffset)) 

def corieToDatetime(corieTime) :
  """Convert corie date to python datetime object"""
  return corieZeroDate + datetime.timedelta(float(corieTime))

def datetimeToCorieTime(t) :
  """Convert python datetime to corie time"""
  return epochToCorieTime( datetimeToEpochTime( t ) )

def simulationToEpochTime(t, simStartTime) :
  """Converts simulation time (in seconds from the start) to epoch time.
  simStartTime is datetime object. t can be an array."""
  if hasattr( t, 'astype' ) :
    return t.astype(np.float64) + datetimeToEpochTime(simStartTime)
  else :
    return np.float64(t) + datetimeToEpochTime(simStartTime)

def epochToSimulationTime(t, simStartTime) :
  """Converts epoch time stamps to simulation time (in seconds) since
  the simulation started.
  simStartTime is datetime object. t can be an array."""
  return t - datetimeToEpochTime(simStartTime)

def epochToMatlabDatenum(t):
  """Converts epoch time to matlab datenum in PST time zone."""
  return 719529.0 + (t - 8*3600) / (24*3600.0)

def matlabDatenumToEpochTime(t):
  """Converts matlab datenum to epoch time. Matlab time is assumed to be in PST which is converted to UTF."""
  return 24*3600*(t - 719529.0) + 8*3600

#def datetimeToSimulationTime(t,startDate) :
  #"""Converts datetime object to seconds versus simulation
  #start time (also a datetime object)"""
  #d = t - startDate
  #return d.total_seconds()
  
#def corieToSimulationTime(corieTime,startDate) :
  #"""Converts corie time to simulation time versus
  #start time (also in corie format)"""
  #return (corieTime-startDate)*86400.0

#class FixedOffsetTZ(datetime.tzinfo):
  #"""Fixed offset in hours east from UTC."""
  #def __init__(self, hours, name):
    #self.offset = datetime.timedelta(hours=hours)
    #self.name = name

  #def __repr__(self):
    #return self.name
  #def utcoffset(self, dt):
    #return self.offset

  #def tzname(self, dt):
    #return self.name

  #def dst(self, dt):
    #return datetime.timedelta(0)
        
#pstTZ = FixedOffsetTZ(-8,'PST')
#utcTZ = FixedOffsetTZ( 0,'UTC')
