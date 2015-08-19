#!/usr/bin/python
"""
Collection of classes and methods to fetch, store and manipulate ADP data.

Tuomas Karna 2012-11-29 / (Modified by Alex Jaramillo 2013-02-04)
"""

import numpy as np
import os
import sys
import datetime

import data.netcdfCacheInterface as netcdfCacheInterface
from data.netcdfCacheInterface import netcdfCacheReader as ncCacheReader
from files.stationFile import StationFile
from crane.data import dataContainer
from crane.data import timeArray
#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
import pdb
#-------------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------
class netcdfCacheADPReader(ncCacheReader) :
  """ Child object overrides getData method to read ADP data"""

  def __init__(self, o) :
     ncCacheReader.__init__(self, o)

  def getData(self, starttime, endtime) :
     if not self.db:
       self.connectToDB()
     if not self.detailsFetched :
       self.getOfferingDetails()
    
     variables = ['height', 'bindepth', self.var]
     t, v, u = netcdfCacheInterface.getncdatastation(self.station, self.nc_offering, starttime, endtime, variables, quality='PD0') 

     if t.shape[0]==0:
       raise Exception("Empty files found for %s - %s time range" % (time.gmtime(starttime), time.gmtime(endtime)))

     # ambigous about following line, but is faster than a query to DB
     binsize = np.max(v['height'])
     t = t[::binsize]
     return t, v, u, binsize

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
def getADPData(offering, sT, eT, var) :
  """Fetches ADP data from the netcdf database, returns a dataContainer.
  
  Arguments:
  station --   name of the station where data is collected
  offering --  code for data offering
  starttime -- time of the first value in datetime
  endtime --   time of the last value in datetime
  var  --      name of the variable to be extracted from cache i.e.: 'bindepth', 
               'alongvel', 'crossvel', 'vel_n', 'height', 'vel_vert', 'vel_mag', 
               'time', 'bindist', 'vel_e', 'bs_avg'
  """
  s, d, b, i = offering.split('.')
  # convert datetime to epochtime with timeArray methods (safe timezone convert)
  sT = datetimeToEpochTime( sT )
  eT = datetimeToEpochTime( eT )
  off = {'location':s, 'msldepth':d, 'bracket':b, 'instrument':i, 'variable': var}
  ncreader = netcdfCacheADPReader( off )
  # Get data as 1D arrays 
  t, v, u, binsize = ncreader.getData(sT, eT)
  sta = StationFile()
  x,y = sta.getLocation( ncreader.station )
  z = v['bindepth']
  # Reshape time, x, y, and data
  zv = -np.reshape(z, (t.size, binsize)).T
  xv = np.zeros((binsize,))
  yv = np.zeros((binsize,))
  xv.fill(x)
  yv.fill(y)
  # next line is a hack because somebody thought it was brilliant to use "uppercase" for DB though
  # every other variable is lowercase and netcdf accepts lowercase only.
  var_key = var.replace("_", " ").lower()
  data = np.reshape(v[var_key], (t.size,binsize)).T[:,None,:]
  # remove time stamps that contain only bad values
  goodIx = np.logical_not( np.all( np.isnan( data[:,0,:] ), axis=0 ) )
  data = data[:,:,goodIx]
  zv = zv[:,goodIx]
  ta = timeArray.timeArray(t[goodIx], 'epoch')
  # Prepare metadata
  meta = {}
  meta['dataType'] =  'profile'
  meta['location'] = ncreader.station
  meta['msldepth'] = str(int(round(int(ncreader.msldepth))))
  meta['bracket'] = ncreader.bracket
  meta['instrument'] = ncreader.instrumenttype
  meta['variable'] = ncreader.var
  meta['tag'] = 'obs' 

  return dataContainer( '', ta, xv, yv, zv, data, [var], coordSys=sta.coordSys, acceptNaNs=True, metaData=meta)


#-------------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------------
def parseCommandLine() :
  from optparse import OptionParser

  parser = OptionParser()
  parser.add_option('-v', '--variable', action='store', type='string',
                      dest='var', help='variable to extract: alongvel, crossvel, vel_n, vel_vert, vel_mag,, bindist, vel_e, bs_avg')
  parser.add_option('-o', '--outDirectory', action='store', type='string',
                      dest='outDir', help='base directory for netCDF file tree '
                                                '(optional)')
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')

  (options, args) = parser.parse_args()

  var           = options.var
  outDir        = options.outDir
  startStr      = options.startStr
  endStr        = options.endStr

  if not var :
    parser.print_help()
    parser.error('variable undefined')
  if not startStr :
    parser.print_help()
    parser.error('start undefined')
  if not endStr :
    parser.print_help()
    parser.error('end undefined')

  sT = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
  eT = datetime.datetime.strptime( endStr ,'%Y-%m-%d')

  dc = getADPData( offering, sT, eT, var)
  dc.setMetaData( 'tag', runTag )
  import data.dirTreeManager as dtm
  rule = dtm.oldTreeRule()
  #rule = dtm.defaultTreeRule()
  dtm.saveDataContainerInTree( dc, path=outDir, rule=rule, dtype=np.float32,
                               overwrite=True  )

if __name__=='__main__' :

  Vs = 'vel_N'
  sT = datetime.datetime(2012, 11, 1)
  eT = datetime.datetime(2012, 11, 30)
  try:
    dc = getADPData('saturn01.1950.A.ADP', sT, eT, Vs)
  except Exception as e:
    print e
  #parseCommandLine()
