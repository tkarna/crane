import os
import sys
import datetime

import numpy as np
from scipy.interpolate import interp1d

from crane.data import dataContainer
from crane.data import timeArray

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
VALID_MIN = -89

#-------------------------------------------------------------------------------
# Functions and classes
#-------------------------------------------------------------------------------

def excludeNaNs( time, data ) :
  """Removes NaNs and Infs from the time and data arrays.

  Args:
    time, data  --  singleton arrays

  Returns:
    time, data  --  cleaned arrays
 """
  goodIx = np.logical_not( np.isnan(time) )
  goodIx = np.logical_and( goodIx, np.logical_not( np.isinf(time) ) )
  goodIx = np.logical_and( goodIx, np.logical_not( np.isinf(data) ) )
  goodIx = np.logical_and( goodIx, np.logical_not( np.isnan(data) ) )
  return time[goodIx], data[goodIx]

def getWeeklyDataDir(baseDir,year,week,tagStr,subdir='run') :
  """Generate a directory where hindcast data is stored.

  Args:
    baseDir -- (str) directory of all hidcasts, e.g. '/home/workspace/ccalmr/hindcasts/'
    year   -- (int) year of interest
    week   -- (int) year of interest
    tagStr -- (str) hindcast tag, e.g. '22'
    subdir -- (str) sub-directory where output files are, e.g. 'run' [default]
  
  Returns:
    dir    -- (str)  directory of the form: baseDir/YYYY-WW-DB/subdir
  """
  yearStr = str(int(year))
  weekStr = '{0:02d}'.format( int(week) )
  return os.path.join( baseDir, '-'.join([yearStr,weekStr,tagStr]), subdir )

class datFileReader(object) :
  def __init__( self, filename, removeBadValues ) :
    self.filename = filename
    self.removeBadValues = removeBadValues
    self.badValue = VALID_MIN
    self.data = dict()
    self.var2d = ['elev','depth']

  def readData( self ) :
    data = np.genfromtxt(self.filename)
    # read nvrt from first header : corieTime, depth, elev, nvrt
    nvrt = int(data[0,3])
    # total lines
    nLines = data.shape[0]
    # number of time steps
    nRec = nLines/(1+2*nvrt)
    # allocate arrays
    self.data['time'] = np.zeros(nRec)
    self.data['depth'] = np.zeros(nRec)
    self.data['elev'] = np.zeros(nRec)
    self.data['temp'] = np.zeros([nRec,nvrt])
    self.data['salt'] = np.zeros([nRec,nvrt])
    #self.data['diff'] = np.zeros([nRec,nvrt])
    #self.data['u'] = np.zeros([nRec,nvrt])
    #self.data['v'] = np.zeros([nRec,nvrt])
    #self.data['w'] = np.zeros([nRec,nvrt])
    self.data['z'] = np.zeros([nRec,nvrt])
    self.nRec = nRec
    self.nvrt = nvrt
    # parse data into separate arrays
    for i in range(nRec) :
      # header
      iLine = i*(1+2*nvrt)
      self.data['time'][i] = data[iLine,0]
      self.data['depth'][i] = data[iLine,1]
      self.data['elev'][i] = data[iLine,2]
      # first data chunk
      iLine = i*(1+2*nvrt) + 1
      chunk = data[iLine:nvrt+iLine]
      self.data['z'][i,:]    = chunk[:,0]
      self.data['salt'][i,:] = chunk[:,1]
      self.data['temp'][i,:] = chunk[:,2]
      #self.data['diff'][i,:] = chunk[:,3]
      # second data chunk
      #iLine = i*(1+2*nvrt) + 1 + nvrt
      #chunk = data[iLine:nvrt+iLine]
      #self.data['z'][i,:] = chunk[:,0]
      #self.data['u'][i,:] = chunk[:,1]
      #self.data['v'][i,:] = chunk[:,2]
      #self.data['w'][i,:] = chunk[:,3]

  def getVariable( self, variable ) :
    if not self.data :
      self.readData()

    outData = self.data[variable]
    if self.removeBadValues and variable != 'time' :
      outData[ outData < self.badValue ] = np.NaN
    return outData
    
  def interpolateVariables( self, variables, zCoord, zRelativeToSurf ) :
    if not self.data :
      self.readData()

    z = self.data['z']
    elev = self.data['elev']
    outData = dict()
    vars3d = []
    for var in variables :
      if var in self.var2d :
        outData[var] = self.data[var]
      else :
        vars3d.append( var )
        outData[var] = np.ones(self.nRec)*-99.0
    if vars3d :
      for i in range(self.nRec) :
        # find indices of strictly inreasing z coords
        incr_ix = np.nonzero(np.diff(z[i,:])>0)[0]
        incr_ix = np.unique(np.hstack(( incr_ix, incr_ix+1 )))
        zz = z[i,incr_ix]
        if len(zz)>0 :
          # bound interpolation point
          # TODO proper depth versus NGDV29 or such
          zInterp = elev[i] + zCoord if zRelativeToSurf else zCoord
          zInterp = min( zz.max(), zInterp )
          zInterp = max( zz.min(), zInterp )
          # interpolate
          for var in vars3d :
            outData[var][i] = interp1d( zz, self.data[var][i,incr_ix] )(zInterp)

    # remove bad values
    if self.removeBadValues :
      for var in variables :
        outData[var][ outData[var] < self.badValue ] = np.NaN
    return outData
    
def getWeeksAndDays(startTime, endTime) :
  """Deduce the weeks and week days for given time period. May encompass several years.
  
  Args:
    startTime -- (datetime) first day
    endTime   -- (datetime) last day
  
  Returns:
    ordinalDays -- (array) all days in the range in ordinal day format (unique int)
    years       -- (array) year for each day in range
    yearDays    -- (array) year day for each day. Jan 1st is day 1.
    weeks       -- (array) "week" for each day in range. 1st week starts Jan 1st.
    weekDays    -- (array) corresponding week day (int) for each day
  """
  # deduce correct years, weeks and week days for the given time range
  ordinalDays = np.arange( startTime.toordinal(), endTime.toordinal()+1 ) # unique int
  years = np.array( [ datetime.datetime.fromordinal( d ).year for d in ordinalDays ] )
  # day of year
  yearDays = ordinalDays - np.array( [ datetime.datetime(year,1,1).toordinal() for year in years ] ) + 1
  weeks = np.floor((yearDays-1)/7)+1 # corresponding week (1st week starts Jan 1st)
  weekDays = yearDays - (weeks-1)*7  # corresponding week day
  return ordinalDays, years, yearDays, weeks, weekDays

# TODO change to use list of offerings
def readDateRange(baseDir, tagStr, station, keys, startTime, endTime, stationX=0, stationY=0, coordSys=None, removeBadValues=False) :
  """Given a range of days, extracts time series from hindcasts.
  Bad values, Nans and Infs are removed from the data.
  Note: If hindcast files are not found, they are skipped silently. Sometimes all data may be bad.
  Thus the requested time period may not be covered and some variables may be entirely missing.
  If no data is retrieved for any variable, None is returned.
  
  Args:
    baseDir -- (str) directory of all hidcasts, e.g. '/home/workspace/ccalmr/hindcasts/'
    tagStr  -- (str) tag for a specific hindcast, e.g. '22' for the databases
    station -- (str) station name
    keys  -- (list of tuples) list of extraction requests ( variable, zCoord, bracket )
             variable (str) e.g. 'elev', 'temp'
             zCoord (float) e.g. -3.30
             bracket (str) e.g. 'A', 'F'
    startTime -- (datetime) first day to extract
    endTime -- (datetime) last day to extract
    stationX,stationY -- (float) horiz. coordinates of station [default 0]
    coordSys -- (str) coordinate system of horiz. coordinates
    removeBadValues -- (bool) If True, values below -90 are replaced by NaNs
    
  Returns:
    dcList -- (list of dataContainers) collection of all the time series associated with station.
              None empty list if extraction fails.
  """

  print 'Extracting station',station
  ordinalDays, years, yearDays, weeks, weekDays = getWeeksAndDays(startTime, endTime)
  print 'yearDays: ', yearDays[0], '...', yearDays[-1]
  
  # low-level extraction without dataContainer for efficiency
  time = [] # list of arrays, one array per day
  data = dict() # dictionary of lists
  # keys are ( var, z, bracket )
  for key in keys :
    data[key] = [] # list of arrays, one array per day
  twoDimVariables = ['elev','depth']
  # all 2d var keys
  keys2d = [ k for k in keys if k[0] in twoDimVariables ]
  # all 3d var keys
  keys3d = [ k for k in keys if k[0] not in twoDimVariables ]
  # 3d var keys sorted by (z,bracket)
  keys3dByDepth = dict()
  for var,z,bra in keys3d :
    if not (z,bra) in keys3dByDepth :
      keys3dByDepth[(z,bra)] = []
    keys3dByDepth[(z,bra)].append( var )
  for i in range(len(ordinalDays)) :
    #/home/workspace/ccalmr/hindcasts/2002-18-16/process/tpoin_elcirc_n003.dat
    dataDir = getWeeklyDataDir(baseDir,years[i],weeks[i],tagStr,subdir='process')
    filename = '{s:s}_elcirc_n{n:03d}.dat'.format(s=station, n=int(weekDays[i]))
    filename = os.path.join(dataDir,filename)
    if not os.path.isfile( filename ) :
      print 'File not found: '+filename+', skipping ...'
      continue
    if i%10 == 0 :
      sys.stdout.write( '\nday ' )
    sys.stdout.write( '{0:3d} '.format(yearDays[i]) )
    sys.stdout.flush()
    rr = datFileReader( filename, removeBadValues )
    time.append( rr.getVariable( 'time' ) )
    for var,z,bra in keys2d :
      data[(var,z,bra)].append( rr.getVariable( var ) )
    for z,bra in keys3dByDepth :
      vars3d = keys3dByDepth[(z,bra)]
      zRelToSurf = True if bra == 'F' else False
      tmp = rr.interpolateVariables( vars3d, z, zRelToSurf )
      for var in vars3d :
        data[(var,z,bra)].append( tmp[var] )
  sys.stdout.write( '\n' )
  
  if len(time) == 0 :
    print 'Extraction failed: no data was retrieved'
    return []
  # time sanity check
  for i in range(len(time)-1) :
    if time[i][-1] > time[i+1][0] :
      print ' overlapping days'
      print '  ',yearDays[i],'ends   at',str(timeArray.epochToDatetime( time[i][-1] ))
      print '  ',yearDays[i+1],'starts at',str(timeArray.epochToDatetime( time[i+1][0] ))
      raise Exception( 'overlapping days... something is wrong with the time stamps' )

  # merge arrays, one array for whole time period
  time = np.hstack( tuple(time) )
  for k in data :
    data[k] = np.hstack( tuple(data[k]) )
  dcList = []
  for var,zCoord,bracket in data :
    ti,di = excludeNaNs( time, data[(var,zCoord,bracket)] )
    if len(ti) == 0 :
      continue
    ta = timeArray(ti,'corie').asEpoch()
    di = di.reshape( (1,1,-1) )

    meta = {}
    meta['location'] = station
    meta['instrument'] = 'model'
    meta['bracket'] = bracket
    meta['msldepth'] = str(int(round(-zCoord*100)))
    meta['dataType'] = 'timeseries'
    meta['variable'] = var
    dc = dataContainer('',ta,stationX,stationY,zCoord,di,[var],coordSys,meta)
    dcList.append( dc )
    # if suspected bad values, print warning
    badValues = np.any( di < VALID_MIN )
    if badValues :
      print 'Warning: bad values in', tagStr, station,zCoord,bracket,var
  return dcList

# ----------- main ------------

def example_usage() :
  """Example usage"""
  # hard-coded stations
  stations = {
              #'sat01':[ 349653, 290854, -7.4 ],
              #'sat03':[ 344241, 287160, -2.4 ],
              #'sat04':[ 359101, 285529, -8.6 ],
              #'hmndb':[ 343541.1, 292041.2, -4.07 ],
              #'tpoin':[ 357302.28000000, 287741.89000000, -2.1 ],
              #'cbnc3':[ 361994.20000000, 287309.11000000, -6.5 ],
              #'sveni':[ 366850.5, 284192.5, -10.8 ],
              'marsh':[ 368985.0, 287756.0, -5.4 ],
              'woody':[ 375759.13000000, 291725.65000000, -2.4 ],
              #'eliot':[ 369566.17000000, 292736.66000000, -13.9 ],
              #'grays':[ 357926.00, 294761.00,  -1.6 ],
              #'tansy':[ 345784.0, 285882.0, -8.4 ]
              }
  # user defined input
  keys = [ ('elev',-1.2,'A'), ('temp',-1.2,'A'), ('temp',-3.0,'F'), ('salt',-1.2,'A'), ('salt',-3.0,'F')]
  startTime = datetime.datetime(2002, 4, 1)
  endTime = datetime.datetime(2002,4,12)
  baseDir = '/home/workspace/ccalmr/hindcasts/'
  
  # process 3 hindcasts
  dbs = [14,16,22]
  for db in dbs :
    # create output dir
    dbStr = 'tmp_db'+str(db)
    outDir = dbStr
    if not os.path.isdir(outDir) :
      os.mkdir(outDir)

    # extract each station
    for sta in stations.keys() :
      x,y,z = stations[sta]
      res = readDateRange(baseDir, str(db), sta, keys, startTime, endTime, stationX=x, stationY=y,coordSys='spcs',removeBadValues=True)
      for dc in res :
        dc.saveAsNetCDF( path=outDir )

if __name__ == '__main__' :
  example_usage()
