#!/usr/bin/python
"""
Collection of classes and methods to fetch, store and manipulate AUV data.

Tuomas Karna 2012-11-29
"""

import numpy as np
import os
import sys
import datetime

import data.netcdfCacheInterface as netcdfCacheInterface
from crane.data import dataContainer
from crane.data import timeArray
#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

transectFilesBaseDir = '/home/users/karnat/projects/AUV_comp'
alongTrBPFile = 'nc_auv_along.bp'
eastTrBPFile = 'nc_auv_east.bp'
westTrBPFile = 'nc_auv_west.bp'
ETMStationFile = 'nc_auv.sta'

# hard-coded dictionary for mission number to missionID
missionID = {
     9  : 805,
     10 : 804,
     12 : 806,
     13 : 807,
     14 : 808,
     15 : 809,
     16 : 810,
     17 : 811,
     20 :  12,
     22 :  11,
     23 :  21,
     24 :  22,
     25 :  23,
     26 :  27,
     27 :  28,
     28 :  31,
     29 :  32,
     30 :  34,
     40 : 821,
     41 : 822,
     44 : 828,
     45 : 838,
     46 : 839,
     47 : 834,
     48 : 835,
     49 : 836,
     50 : 845,
     52 : 841,
     53 : 842,
     65 : 853,
     66 : 856,
     67 : 852,
     68 : 854,
     69 : 855,
     70 : 857,
     71 : 859,
     72 : 858,
     73 : 861,
     74 : 860,
     }

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def getAUVData( missionNB, var ) :
  """Fetches AUV data from the netcdf database, returns a dataContainer.
  MissionNB is the number Craig uses for remus data.
  """
  x,y,dep,t,v = netcdfCacheInterface.getAUVData( missionID[missionNB], var )
  # reshape for dataContainer
  x = x[None,:]
  y = y[None,:]
  z = dep[None,:] # z coordinate versus free surface
  data = np.reshape( v, (1,1,-1) )
 
  ta = timeArray.timeArray( t, 'epoch' )
  
  print 'AUV mission',missionNB, ta.getDatetime(0),'->', ta.getDatetime(-1)

  meta = {}
  meta['dataType'] = 'track'
  meta['location'] = 'AUV'+str(missionNB)
  meta['instrument'] = 'AUV'
  meta['bracket'] = 'F'
  meta['variable'] = var
  dc = dataContainer.dataContainer( '', ta,x,y,z,data,[var],coordSys='spcs', metaData=meta )

  return dc

#-------------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------------
def parseCommandLine() :
  from optparse import OptionParser

  parser = OptionParser()
  parser.add_option('-v', '--variable', action='store', type='string',
                      dest='var', help='variable to extract: elev, temp, or salt')
  parser.add_option('-o', '--outDirectory', action='store', type='string',
                      dest='outDir', help='base directory for netCDF file tree '
                                                '(optional)')
  #parser.add_option('-s', '--start', action='store', type='string',
                      #dest='startStr', help='Date to start processing')
  #parser.add_option('-e', '--end', action='store', type='string',
                      #dest='endStr', help='Date to end processing')
  parser.add_option('-i', '--missionNB', action='store', type='int',
                      dest='missionNB', help='mission ID number, e.g. 46')

  (options, args) = parser.parse_args()

  var           = options.var
  outDir        = options.outDir
  #startStr      = options.startStr
  #endStr        = options.endStr
  missionNB     = options.missionNB
  runTag = 'AUV'

  if not var :
    parser.print_help()
    parser.error('variable undefined')
  #if not outDir :
    #parser.print_help()
    #parser.error('outDir   undefined')
  if not missionNB :
    parser.print_help()
    parser.error('missionNB undefined')
  #if not startStr :
    #parser.print_help()
    #parser.error('startStr undefined')
  #if not endStr :
    #parser.print_help()
    #parser.error('endStr   undefined')

  #startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
  #endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')

  dc = getAUVData( missionNB, var )
  dc.setMetaData( 'tag', runTag )
  import data.dirTreeManager as dtm
  rule = 'singleFile'
  dtm.saveDataContainerInTree( dc, rootPath=outDir, rule=rule, dtype=np.float32,
                               overwrite=True  )

if __name__=='__main__' :
  parseCommandLine()
