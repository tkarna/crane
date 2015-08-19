"""
Reads Coastal Upwelling Index ACSII files.

Tuomas Karna 2012-10-18
"""

#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import datetime
from crane.data import timeArray
from data.dataContainer import dataContainer

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
DEFAULT_FILE = '/home/workspace/ccalmr/data/external/index/CUI/ERD_upwelling_125_45_1967_2012.txt'

#-------------------------------------------------------------------------------
# Classes and functions
#-------------------------------------------------------------------------------
class cuiParser(object) :
  def __init__( self, filename=None ) :
    if not filename :
      filename = DEFAULT_FILE
    self.filename = filename
    self.time = None
    self.data = None
    self.coordinates = None
    self.lat = None
    self.lon = None
    self.badValue = None
    self.fileIsRead = False

  def readData( self ) :
    """Reads data from file to internal variables"""
    f = open(self.filename)
    self.time = []
    self.data = []
    for line in f :
      if line.find('BAD FLAG') > 0 :
        self.badValue = float( line.split(':')[1].strip() )
      if line.find('LONGITUDE') > 0 :
        self.lon = line.split(':')[1].strip()
      if line.find('LATITUDE') > 0 :
        self.lat = line.split(':')[1].strip()
      if len(line) >6 and line[2] == '-' and line[6] == '-' :
        parts = line.rsplit(None,1)
        # data line
        timeStamp = datetime.datetime.strptime( parts[0], '%d-%b-%Y %H' )
        t = timeArray.datetimeToEpochTime( timeStamp )
        self.time.append( t )
        val = float( parts[1] )
        self.data.append( val )

    self.time = np.array( self.time )
    self.data = np.array( self.data )
    # remove bad values
    if self.badValue :
      goodIx = self.data != self.badValue
      self.time = self.time[goodIx]
      self.data = self.data[goodIx]
    self.fileIsRead = True

  def getDataContainer( self,startTime=None, endTime=None ) :
    """Returns dataContainer with cui data. Optional startTime, endTime are used
    to crop the time series.
    """
    if not self.fileIsRead :
      self.readData()
    station = 'cui'
    x = y = z = 0
    coordSys = ''
    if self.lat and self.lon :
      station = self.lat+self.lon
      # coordinates in deg north, east
      sign = 1 if self.lat[-1] == 'N' else -1
      lat = sign*float( self.lat[:-1] )
      sign = 1 if self.lon[-1] == 'E' else -1
      lon = sign*float( self.lon[:-1] )
      self.coordinates = [ lon,lat ]
      x,y = self.coordinates
      coordSys = 'lonlat'
    meta = {}
    meta['dataType'] = 'timeseries'
    meta['location'] = station
    meta['msldepth'] = '0'
    meta['variable'] = 'cui'
    #meta['instrument'] = 'CUI'
    meta['bracket'] = 'A'
    meta['tag'] = 'obs'
    dc = dataContainer.fromTimeSeries( '', self.time, self.data, ['cui'],
                                       x,y,z, 'epoch', coordSys, metaData=meta )
    if startTime and endTime :
      dc = dc.timeWindow( startTime, endTime, includeEnd=True )
    return dc

if __name__=='__main__' :
  st = datetime.datetime(2000,1,1)
  et = datetime.datetime(2012,8,31)
  dc = cuiParser().getDataContainer(st,et)
  #dc = cuiParser().getDataContainer()
  print dc
