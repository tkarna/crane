"""
Methods for reading a station.sta file and getting station locations.
"""
#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import os

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
FILE = '/home/workspace/ccalmr/hindcasts/reference/real_stations.sta'

#-------------------------------------------------------------------------------
# Classes and functions
#-------------------------------------------------------------------------------
class StationFile(object):
  """Handles reading of station.sta files."""
  def __init__(self, path=FILE, description=None, coordSys='spcs'):
    """Just returns an empty object"""
    self.path = path 
    self.description = description 
    self.coordSys = coordSys 
    self.nStations = ''
    self.stations = {}

  def __iter__( self ) :
    return self.stations.__iter__()

  def readFileFromDisk(self, file=None):
    """Reads a station.sta file.

    Args:
        file -- String of path to station.sta file
    """ 
    if file is None: file = self.path
    if not os.path.isfile(file):
      raise Exception('File does not exist: '+file)
    self.file = file

    f = open(file, 'r')
    self.nStations = int(f.readline().rstrip('\n'))
    for p in xrange(self.nStations):
      name = f.readline().rstrip() 
      coords = f.readline().rstrip().split()
      x = coords[0]
      y = coords[1]
      self.stations[name] = (float(x),float(y))

  def getLocation(self, station):
    """Returns xy location of a station. Returns None if station not found.

    Args:
      station -- String of the name of the station
    Returns
      xy -- Tuple (x,y) of station location
    """
    if not self.stations :
      self.readFileFromDisk()
    return self.stations.get(station)

  def addStation(self, name, xy):
    """Add a station to this collection.

    Args:
      name -- String of station name
      xy -- Tuple of station location
    """
    if type(xy) != tuple:
      print ('xy must be an xy tuple: (x,y)')
    if len(xy) != 2:
      print ('xy must have be a tupe of len 2: (x,y)')
    if type(name) != string:
      print ('name must be a string')
    self.stations[name] = xy

  def deleteStation(self, station):
    """Removes a station from the object.

    Args:
      station -- String of name of station to remove
    """
    if not station in self.stations:
      print ('Station not removed, not in collection')
    del self.stations[station]

  def writeToFile(self, path):
    """Write object to station.sta file"""
    nf = open(path, 'w')
    nf.write('%d\n' % self.nStations)
    for station in self.stations:
      nf.write(station+'\n')
      xy = self.stations[station]
      nf.write('%s %s\n' % (xy[0], xy[1]))
    nf.close()

#-------------------------------------------------------------------------------
# Test
#-------------------------------------------------------------------------------
def test_StationFile():
  """Test class """
  print 'Read file'
  sf = StationFile()
  sf.readFileFromDisk()

  print 'Write copy of file'
  testFile = 'test_file'
  sf.writeToFile(testFile)

  print 'Test copy'
  nf = StationFile(testFile)
  nf.readFileFromDisk()
  if nf.nStations != sf.nStations:
    print 'Test failed: Number of stations not equal'
    print 'Original: '+str(sf.nStations)
    print 'Test: '+str(nf.nStations)
    return
  for station in sf.stations:
    if sf.stations[station] != nf.stations[station]:
      print 'Test failed: Station locations differ.'
      og = sf.stations[station]
      print 'Original: %f %f' % (og[0], og[1])
      bg = nf.stations[station]
      print 'Test copy: %f %f' % (bg[0], bg[1])
      return
  os.system('rm -f test_file')
  print 'Test passed' 

#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
if __name__ == '__main__':
  test_StationFile()
