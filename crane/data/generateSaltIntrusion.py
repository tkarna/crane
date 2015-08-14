#!/usr/bin/env python
"""
Script for generating bottom friction Cd plots.

Tuomas Karna 2013-03-01
"""

import os
import sys
import numpy as np
from scipy.interpolate import interp1d
import time as timeMod

from data.timeArray import *
from data.dataContainer import dataContainer
from data.dirTreeManager import netcdfTree, oldTreeRule, defaultTreeRule
from plotting.transectPlot import generateTransectFromDataContainer

TRANSECT_NAME = 'mainChannel'
TRANSECT_BPFILE = '/home/workspace/users/pturner/db29/processing.dev/scripts.working/intrusion_length.bp'

def computeSaltIntrusion( transectDC, salt_threshold_list ) :
  print 'Generating salt intrusion length... ',
  sys.stdout.flush()
  t0 = timeMod.clock()
  nTime = len(transectDC.time)
  nThresholds = len( salt_threshold_list )
  sil = [ np.zeros((nTime,)) for i in xrange(nThresholds) ]
  # for each time step
  for it in range(nTime) :
    # convert transect to array
    Xalong, Z, salt, time, uniqueXYCoords = generateTransectFromDataContainer(transectDC,it)
    Xalong = Xalong[0,:]
    # compute max salt in each column
    maxSalt = salt.max(axis=0)
    for i,salt_threshold in enumerate(salt_threshold_list) :
      binaryIx = maxSalt > salt_threshold
      x_threshold = 0
      if binaryIx[-1] : # whole transect is within threshold
        x_threshold = Xalong[-1]
      elif not binaryIx.any() :
        x_threshold = Xalong[0]
      else :
        # find max S > salt_threshold
        thIx = np.nonzero( maxSalt > salt_threshold )[0]
        thIx = thIx.max()
        # interpolate distance
        #print thIx, maxSalt[thIx]
        x_threshold = interp1d( maxSalt[[thIx+1,thIx]], Xalong[[thIx+1,thIx]] ) (salt_threshold)
      sil[i][it] = x_threshold/1000. # convert to km
  print 'duration ',timeMod.clock()-t0,'s'

  output_dc_list = []
  for i,salt_threshold in enumerate(salt_threshold_list) :
    # generate dataContainer
    meta = {}
    meta['location'] = transectDC.getMetaData('location')
    meta['instrument'] = 'model'
    meta['variable'] = '''sil_%d''' % (salt_threshold,)
    meta['dataType'] = 'sil'
    meta['tag'] = transectDC.getMetaData('tag')
    data = sil[i][None,None,:]
    silDC = dataContainer( '', transectDC.time, 0,0,0, data,
                            [meta['variable']], coordSys='',metaData=meta)
    # compute daily max
    dailyMax = []
    dailyTime = []
    dayBegin = silDC.time.getDatetime(0).replace(minute=0, hour=0, second=0, microsecond=0)
    while dayBegin < silDC.time.getDatetime(-1) :
      dayEnd = dayBegin + datetime.timedelta(hours=24.8)
      try :
        dailySIL = silDC.timeWindow(dayBegin,dayEnd)
        m = dailySIL.data.max()
        dailyMax.append( m )
        dailyMax.append( m )
        dailyTime.append( dailySIL.time.array[0] )
        dailyTime.append( dailySIL.time.array[-1] )
      except Exception as e :
        print 'Cannot compute SIL: day missing, skipping'
        print e
      dayBegin = dayEnd
    data = np.array( dailyMax )[None,None,:]
    ta = timeArray( np.array(dailyTime), 'epoch' )
    meta = {}
    meta['location'] = transectDC.getMetaData('location')
    meta['instrument'] = 'model'
    meta['variable'] = '''max_sil_%d''' % (salt_threshold,)
    meta['dataType'] = 'sil'
    meta['tag'] = transectDC.getMetaData('tag')
    dailyMaxDC = dataContainer( '', ta, 0,0,0, data,
                            [meta['variable']], coordSys='',metaData=meta)
    output_dc_list.append(silDC)
    output_dc_list.append(dailyMaxDC)
  return output_dc_list

def readSILTransect( tag, startTime, endTime ) :
  # load transect dataContainer
  rule = oldTreeRule()
  #rule = defaultTreeRule()
  tree = netcdfTree( dataType='transect', tag=tag, location=TRANSECT_NAME, variable='salt', rule=rule )
  dc = tree.read(startTime,endTime)
  if dc == None:
    raise Exception('salt intrusion_length transect could no be read')
  return dc

#-------------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------------
def parseCommandLine() :
  from optparse import OptionParser

  parser = OptionParser()
  parser.add_option('-r', '--runTag', action='store', type='string',
                      dest='runTag', help='Run tag, used as a label in post-proc.')
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
  parser.add_option('-d', '--dataDirectory', action='store', type='string',
                      dest='dataDir', help='directory where model outputs are stored')
  parser.add_option('-v', '--salinityThreshold', action='store', type='string',
                      dest='salinityThreshold', help='PSU value for determining presence of salinity along the transect line (default %default). Multiple thresholds can be given in a list e.g. 1.0,5.0', default='1.0')
  parser.add_option('-t', '--buildPointFile', action='store', type='string',
                      dest='bpFile', help='text file (*.bp) containing the transect (x,y) coordinates (optional)',
                             default=None)
  parser.add_option('', '--variable', action='store', type='string',
                      dest='var', help='Define alternative output file process, e.g. salt.70 (default %default)', default='salt')
  parser.add_option('-c', '--modelCoordSys', action='store', type='string',
                      dest='modelCoordSys', default='spcs',
                      help='horizontal coordinate system used in model: '
                           'spcs or utm (Default: %default)')
  parser.add_option('-C', '--read-netcdf', action='store_true',
                    dest='readNetcdf',
                    help='Extract from SELFE netcdf output files instead of SELFE binary files (default %default)',default=False)

  (options, args) = parser.parse_args()

  runTag        = options.runTag
  startStr      = options.startStr
  endStr        = options.endStr
  dataDir       = options.dataDir
  bpFile        = options.bpFile
  modelCoordSys = options.modelCoordSys
  salinityThreshold = options.salinityThreshold
  readNetcdf    = options.readNetcdf
  var           = options.var

  salt_threshold_list = [ float(vStr) for vStr in salinityThreshold.split(',')]

  if not bpFile :
    # TODO move to shared location
    bpFile = TRANSECT_BPFILE

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

  startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
  endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')

  print 'Parsed options:'
  print ' - time range:',str(startTime),'->', str(endTime)
  print ' - salinity threshold(s):',salt_threshold_list
  print ' - dataDir',dataDir
  print ' - variable',var
  print ' - SELFE output format:','netCDF' if readNetcdf else 'binary'
  print ' - runTag',runTag
  print ' - transect file',bpFile
  print ' - model coord system', modelCoordSys

  # Extract 
  varList = [var]
  name = TRANSECT_NAME
  dcs = []
  if readNetcdf :
    from data.ncExtract import extractTransectForBPFile
    dcs = extractTransectForBPFile( bpFile, dataDir, varList,
                                    startTime, endTime, name )
  else :
    from data.extractTransect import extractTransectForBPFile
    dcs = extractTransectForBPFile( bpFile, dataDir, varList,
                                    startTime, endTime, name=name,
                                    modelCoordSys=modelCoordSys )
  for dc in dcs :
    dc.setMetaData( 'tag', runTag )

  # compute SIL
  silDCs = computeSaltIntrusion( dcs[0], salt_threshold_list )
  for dc in silDCs :
    print dc

  import data.dirTreeManager as dtm
  #rule = dtm.oldTreeRule()
  rule = dtm.defaultTreeRule()
  dtm.saveDataContainerInTree( silDCs, rule=rule, dtype=np.float32,
                               overwrite=True )

if __name__=='__main__' :
  parseCommandLine()
