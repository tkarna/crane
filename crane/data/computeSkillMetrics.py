#!/usr/bin/env python
"""
Computes various skill metrics comparing model time series against reference.

Tuomas Karna 2013-11-07
"""
import os
import numpy as np
import datetime
import sys
from optparse import OptionParser
import data.statistics as stat
from data.stationCollection import tinyDB
from data.collection import *
from data.timeSeriesFilters import *
import data.dirTreeManager as dtm
from files.csvStationFile import *
import traceback

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def getSurrogateModelDiff(surTag, modTag, location, z=None, k=None) :
  mod = getTSFromProfile(modTag, location, z, k)
  sur = getTSFromProfile(surTag, location, z, k)
  # compute difference
  diff = mod.computeError( sur )
  tagStr = sur.getMetaData('tag').split('-')[0]+'-'+mod.getMetaData('tag').split('-')[0]
  diff.setMetaData('tag',tagStr)
  return diff

def getDataArray( signal ) :
  """Returns a singleton data array from dataContainer, returns numpy arrays
  directly.
  """
  if isinstance( signal, np.ndarray ):
    return signal
  return signal.data.ravel()

def computeStats( reference, signal ):
  """Returns statistics comparing signal to reference in a dictionary"""
  r = getDataArray(reference)
  s = getDataArray(signal)
  d = {}
  d['bias'] = stat.bias( r, s )
  d['rmse'] = stat.rootMeanSquareError( r, s, centered=False )
  d['nmse'] = stat.normalizedMeanSquareError( r, s )
  d['stde'] = stat.standardDeviation( s-r )
  d['std'] = stat.standardDeviation( s )
  d['ioa'] = stat.indexOfAgreement( r,s )
  d['corr'] = stat.correlationCoefficient( r,s )
  d['murphy'] = stat.murphySkillScore( r,s )
  return d

def getStats( reference, other ) :
  """Computes statistics comparing reference and other signals"""
  r,o = reference.alignTimes( other )
  stats = computeStats( r, o )

  lp_r = removeTides( r )
  lp_o = removeTides( o )
  lp_stats = computeStats( lp_r, lp_o )

  hp_r = r.computeError(lp_r)
  hp_o = o.computeError(lp_o)
  hp_r2, hp_o2 = hp_r.alignTimes(hp_o)
  hp_stats = computeStats(hp_r2, hp_o2)
  return {'original':stats,'lowpass':lp_stats,'highpass':hp_stats}

def printStats( fid, label, stats, metrics ) :
  """Print one line of statistics summary"""
  fmt = '%7.3f'
  content = []
  content.append( '%20s'%label )
  if stats != None :
    for dataSetName,measure in metrics :
      content.append( fmt%stats[dataSetName][measure] )
  else :
    for dataSetName,measure in metrics :
      content.append( '%7s'%'--' )
  contentStr = ' '.join(content)+'\n'
  fid.write(contentStr)

def printStatsLegend( fid, metrics ) :
  """Print legend line of statistics summary"""
  fmt = '%7.3f'
  content = []
  content.append( '%20s'%'name' )
  for dataSetName,measure in metrics :
    setAbbrev = '' if dataSetName in ['original', 'all'] else dataSetName[:2]+'_'
    content.append( '%3s'%setAbbrev+measure )
  contentStr = ' '.join(content)+'\n'
  fid.write(contentStr)

def printStationStats( fid, allStats, refTag, runTags ) :
  """Print summary of station statistics to given open file"""
  metrics = [('original','rmse'),
             ('original','bias'),
             ('original','nmse'),
             ('original','stde'),
             ('lowpass','stde'),
             ('highpass','stde'),
             ('highpass','std'),
             ]
  fid.write('\nStation data sets\n\n')

  # find all unique station,variable,depth combinations
  entries = [ (loc,var,msld) for tag,loc,var,msld in allStats.getTuples() ]
  entries = uniqueList(entries)
  # sort by var,loc,msldepth
  entries = sorted(entries, key=lambda s: s[2])
  entries = sorted(entries, key=lambda s: s[0])
  entries = sorted(entries, key=lambda s: s[1])
  for loc,var,msld in entries:
    fid.write( '{0:s} {1:s} {2:s}\n'.format(var,loc,msld) )
    printStatsLegend(fid,metrics)
    stats = allStats.getSample(tag=refTag,location=loc,
                               variable=var,msldepth=msld)
    printStats(fid,refTag,stats,metrics)
    for runTag in runTags :
      stats = allStats.getSample(tag=runTag,location=loc,
                                variable=var,msldepth=msld)
      printStats(fid,runTag,stats,metrics)

def printMurphyStats( fid, allStats, refTag, runTags ) :
  """Print summary of station statistics to given open file"""
  metrics = [('original','MS'),
             ]
  fid.write('\nMurphy scores\n\n')
  # find all unique station,variable,depth combinations
  entries = [ (loc,var,msld) for tag,loc,var,msld in allStats.getTuples() ]
  tag = allStats.getTuples()[0][0]
  entries = uniqueList(entries)
  # sort by var,loc,msldepth
  entries = sorted(entries, key=lambda s: s[2])
  entries = sorted(entries, key=lambda s: s[0])
  entries = sorted(entries, key=lambda s: s[1])
  for loc,var,msld in entries:
    fid.write( '{0:s} {1:s} {2:s}\n'.format(var,loc,msld) )
    #printStatsLegend(fid,metrics)
    stats = allStats.getSample(tag=tag,location=loc,
                               variable=var,msldepth=msld)
    printStats(fid,tag,stats,metrics)

def printGlobalStats( fid, allStats, refTag, runTags ) :
  """Print summary of global statistics to given open file"""
  metrics = [('all','rmse'),
             ('all','bias'),
             ('all','nmse'),
             ('all','corr'),
             ('all','ioa'),
             ]
  fid.write('\nCombined data sets\n\n')

  # find all unique station,variable,depth combinations
  entries = [ var for tag,var in allStats.getTuples() ]
  entries = uniqueList(entries)
  # sort by var
  entries = sorted(entries)
  for var in entries:
    fid.write( '{0:s}\n'.format(var) )
    printStatsLegend(fid,metrics)
    stats = allStats.getSample(tag=refTag, variable=var)
    printStats(fid,refTag,stats,metrics)
    for runTag in runTags :
      stats = allStats.getSample(tag=runTag, variable=var)
      printStats(fid,runTag,stats,metrics)

#-------------------------------------------------------------------------------
# Main routine
#-------------------------------------------------------------------------------
def processStats(refTag, runTags, stationFile, startTime=None, endTime=None,
                 outFile=None):

  if outFile == None :
    fid = sys.stdout
  else :
    path = os.path.split(outFile)[0]
    if path: # create dir
      if os.path.exists(path) :
        if not os.path.isdir(path) :
          raise Exception( 'file with same name exists',path )
      else :
        os.makedirs(path)
    fid = open(outFile,'w')

  defRule = dtm.defaultTreeRule()

  if isinstance(stationFile, str) and stationFile[-4:]=='.csv' :
    # read station file and convert to a filter (list of dict)
    stations = csvStationFileWithDepth()
    stations.readFromFile(stationFile)
    stationFilter = []
    for loc,x,y,z,zType,var in stations.getTuples():
      msldepth = str(int(round(abs(z)*100)))
      stationFilter.append({'location':loc, 'msldepth':msldepth,
                            'variable':var})
  elif isinstance(stationFile, list):
    # assume that user provided the filter directly
    stationFilter = stationFile
  else:
    raise Exception('Unknown stationFile format')

  # loop over data sets and store statistics in structure:
  # stationStats[('db31','saturn01','salt')]['lowpass']['rmse']
  stationStats=tinyDB(['tag','location','variable','msldepth'])
  # concatenate all time series to one array for global metrics
  # allData[('db31','salt','sig')] # aligned signal
  # allData[('db31','salt','ref')] # aligned reference signal
  allData = {}
  murphyScores = tinyDB(['tag','location','variable','msldepth'])
  #for loc,x,y,z,zType,var in stations.getTuples() :
    #print ' * ',loc,z,zType,var
  for filt in stationFilter :
    loc = filt['location']
    var = filt['variable']
    msldepth = filt['msldepth']
    print ' * ',loc,msldepth,var
    try :
      #msldepth = str(int(round(abs(z)*100)))
      # load reference and other time series
      ref = dtm.getDataContainer(tag=refTag,location=loc,dataType='timeseries',
                                  msldepth=msldepth,variable=var,rule=defRule,
                                  startTime=startTime,endTime=endTime)
      stationStats.addSample(getStats(ref,ref), tag=refTag,
                      location=loc, variable=var, msldepth=msldepth)
      modSignals = {}
      for runTag in runTags :
        try :
          sig = dtm.getDataContainer(tag=runTag,location=loc,
                                    dataType='timeseries',
                                    msldepth=msldepth,variable=var,rule=defRule,
                                    startTime=startTime,endTime=endTime)
          stationStats.addSample(getStats(ref,sig), tag=runTag,
                          location=loc, variable=var, msldepth=msldepth)
          s,r = sig.alignTimes( ref )
          allData.setdefault((runTag,var,'ref'),[]).append(r.data.flatten())
          allData.setdefault((runTag,var,'sig'),[]).append(s.data.flatten())
          modSignals[runTag] = sig
        except Exception as e:
          print ' computing model {0:s} stats failed, skipping...'.format(runTag)
          traceback.print_exc(file=sys.stdout)

      if len(modSignals)==2:
          # compute Murphy score, do 3-way merge of time steps
          sig1 = modSignals[runTags[0]]
          sig2 = modSignals[runTags[1]]
          st = max(ref.time.getDatetime(0), sig1.time.getDatetime(0), sig2.time.getDatetime(0))
          et = min(ref.time.getDatetime(-1), sig1.time.getDatetime(-1), sig2.time.getDatetime(-1))
          if st > et:
              raise Exception('could not merge 3 data sets')
          r = ref.timeWindow(st, et, includeEnd=True)
          s1 = sig1.timeWindow(st, et, includeEnd=True)
          s2 = sig2.timeWindow(st, et, includeEnd=True)
          s1, s2 = s1.alignTimes(s2)
          s1, r = s1.alignTimes(r)
          st2 = max(r.time.getDatetime(0), s1.time.getDatetime(0), s2.time.getDatetime(0))
          et2 = min(r.time.getDatetime(-1), s1.time.getDatetime(-1), s2.time.getDatetime(-1))
          newTime = r.time.window(st2, et2)
          r = r.interpolateInTime(newTime, True)
          s1 = s1.interpolateInTime(newTime, True)
          s2 = s2.interpolateInTime(newTime, True)
          #s1, r = s1.alignTimes(r)
          #s2 = s2.interpolateInTime(s1.time, True)
          r = r.data.ravel()
          s1 = s1.data.ravel()
          s2 = s2.data.ravel()
          mur = stat.murphySkillScore(r, s1, reference_model=s2)
          s = {'original':{}}
          s['original']['MS'] = mur
          murphyScores.addSample(s, tag=runTags[0]+' vs '+runTags[1],
                                 variable=var, location=loc, msldepth=msldepth)
    except Exception as e:
      print ' computing stats failed, skipping...'
      traceback.print_exc(file=sys.stdout)

  # concatenate allData
  for k in allData:
    allData[k] = np.hstack( tuple(allData[k]) )

  # compute global statistics
  globalStats = tinyDB(['tag','variable'])
  allVars = uniqueList([t[1] for t in allData ])
  for var in allVars:
    for tag in runTags :
      try :
        obs = allData[(tag,var,'ref')]
        mod = allData[(tag,var,'sig')]
        s = {}
        s['all']=computeStats(obs,mod)
        globalStats.addSample(s,tag=tag,variable=var)
      except Exception as e:
        print e

  # print statistics
  if outFile != None :
    print 'writing output to', outFile

  printGlobalStats( fid, globalStats, refTag, runTags )

  printStationStats( fid, stationStats, refTag, runTags )

  if len(murphyScores.getKeys()) > 0:
    printMurphyStats( fid, murphyScores, refTag, runTags )

#-------------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------------
def parseCommandLine() :

  usage = ('Usage: %prog -r refTag -t [csvStationFile] runID1 runID2 ...\n')

  parser = OptionParser(usage=usage)
  parser.add_option('-r', '--refTag', action='store', type='string',
                      dest='refTag', help='Name of data set to use as reference')
  parser.add_option('-t', '--csvStationFile', action='store', type='string',
                      dest='csvStationFile', help='file that defines station coordinates and\
                      variables to compare')
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing (optional)')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing (optional)')
  parser.add_option('-o', '--output-file', action='store', type='string',
                      dest='outFile', help='File where statistics are stored (optional, default:print to stdout)')

  (options, args) = parser.parse_args()

  refTag = options.refTag
  csvStationFile = options.csvStationFile
  startStr = options.startStr
  endStr = options.endStr
  outFile = options.outFile

  runTags = args

  if refTag == None:
    parser.print_help()
    parser.error('refTag undefined')
  if csvStationFile == None:
    parser.print_help()
    parser.error('csvStationFile undefined')
  
  startTime = None
  if startStr :
    startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
  endTime = None
  if endStr :
    endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')

  print 'Parsed options:'
  print ' - reference',refTag
  print ' - runTags',runTags
  print ' - using csvStationFile',csvStationFile
  if startTime or endTime:
    print ' - time range:',str(startTime),'->', str(endTime)
  if outFile :
    print ' - output file:',outFile
    

  processStats(refTag,runTags,csvStationFile,startTime,endTime,outFile)
  
if __name__=='__main__' :
  parseCommandLine()
