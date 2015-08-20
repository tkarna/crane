
import numpy as np
import datetime
import os
import sys
from optparse import OptionParser

import data.stationCollection as StationCollection
from data.loadHindcastStations import readDateRange
from files.stationFile import StationFile

from plotting.plot import Plots

def createDir(path) :
  if os.path.exists(path) :
    if not os.path.isdir(path) :
      raise Exception( 'file with same name exists',path )
  else :
    os.makedirs(path)
  return path

#-------------------------------------------------------------------------------
# Main routine
#-------------------------------------------------------------------------------
def skillPlots(runTag, imgDir, modelDir, hindcastDirs, startTime, endTime, stationFile=None) :

  imgDir = createDir(imgDir)

  baseDir = '/home/workspace/ccalmr/hindcasts/'
  removeBadValues = True

  # limits for time series plots
  ts_ylim = {'elev':[-2.0,10.0], 'temp':[5,20], 'salt':[-0.45,35]}
  ts_err_ylim = {'elev':[-0.6,0.6], 'temp':[], 'salt':[]}
  # limits for error histograms
  eh_lim = {'elev':[-1.0,1.0], 'temp':[-8,8], 'salt':[-32,32]}
  eh_bins = 60
  # limits for spectral plot
  esp_ylim = [0.0, 0.2]

  #-------------------------------------------------------------------------------
  # Part I  : set station information
  #-------------------------------------------------------------------------------

  # get "all" station names and coordinates
  sta = StationFile()
  sta.readFileFromDisk(stationFile)
  stationNames = sta.stations.keys()
  stationXY = sta.stations
  stationX = dict( (s,stationXY[s][0]) for s in stationXY )

  # all stations for elevation comparison
  elevStations = [ 'hmndb', 'tpoin', 'skaw1', 'bvao3', 'lonw1', 'stho3', 'vanw1', 'saturn06', 'prto3', 'wilo3', 'bono3' ]
  presStations = [ 'tansy', 'woody', 'cbnc3', 'marsh', 'saturn01', 'saturn04', 'grays','eliot' ]
  # ~distance from mouth
  elevStationX = dict( (s, stationX[s]) for s in elevStations if s in stationX )

  # filter for surface salinity stations
  surfSaltFilter = [
                    (None, 'saturn01', 'salt', '0'), # 0, 740, 1950
                    (None, 'saturn02', 'salt', '100'), # 100, 600, 1000, 1100, 1600, 2100, 3500
                    (None, 'saturn03', 'salt', '240'), # 240, 820, 1300
                    (None, 'saturn04', 'salt', '30'), # 30, 830, 860, 910
                    (None, 'am169', 'salt', '260'), # 260, 1100, 1130, 1430
                    (None, 'red26', 'salt', '330'), # 330, 750, 900
                    (None, 'ogi01', 'salt', '80'), # 80 500 1100 5000
                    (None, 'nh10', 'salt', '200'), # 200 1000 6000 7300
                  ]

  # filter for bottom salinity stations
  botSaltFilter = [
                    (None, 'saturn01', 'salt', '1950'), # 0, 740, 1950
                    (None, 'saturn02', 'salt', '3500'), # 100, 600, 1000, 1100, 1600, 2100, 3500
                    (None, 'saturn03', 'salt', '1300'), # 240, 820, 1300
                    (None, 'saturn04', 'salt', '910'), # 30, 830, 860, 910
                    (None, 'saturn05', 'salt', None), # 250
                    (None, 'saturn06', 'salt', None), # 50
                    (None, 'saturn07', 'salt', None), # 60
                    (None, 'am169', 'salt', '1430'), # 260, 1100, 1130, 1430
                    (None, 'red26', 'salt', '900'), # 330, 750, 900
                    (None, 'tansy', 'salt', None), # 840
                    (None, 'jetta', 'salt', None), # 640
                    (None, 'sandi', 'salt', None), # 790
                    (None, 'dsdma', 'salt', None), # 730 820
                    (None, 'coaof', 'salt', None), # 320 210
                    (None, 'grays', 'salt', None), # 160
                    (None, 'cbnc3', 'salt', None), # 650 670 900
                    (None, 'sveni', 'salt', None), # 1080
                    (None, 'marsh', 'salt', None), # 540
                    (None, 'eliot', 'salt', None), # 1390
                    (None, 'yacht', 'salt', None), # 560 590
                    (None, 'ogi01', 'salt', '5000'), # 80 500 1100 5000
                    (None, 'nh10', 'salt', '6000'), # 200 1000 6000 7300
                    (None, 'yb101', 'salt', None), # 410 320
                    (None, 'seahs', 'salt', None), # 100
                    (None, 'lwsck', 'salt', None), # 420
                    (None, 'lght6', 'salt', None), # 290
                    (None, 'lght2', 'salt', None), # 440
                    (None, 'coaww', 'salt', None), # 440
                    (None, 'chnkr', 'salt', None), # 0
                    (None, 'chnke', 'salt', None), # 260
                    (None, 'abpoa', 'salt', None), # 410
                  ]

  # ~distance from mouth
  surfSaltStationX = dict( (t[1], stationX[ t[1] ]) for t in surfSaltFilter if t[1] in stationX )
  botSaltStationX = dict( (t[1], stationX[ t[1] ]) for t in botSaltFilter if t[1] in stationX )

  surfTempFilter = [
                    (None, 'saturn01', 'temp', '0'), # 0, 740, 1950
                    (None, 'saturn02', 'temp', '100'), # 100, 600, 1000, 1100, 1600, 2100, 3500
                    (None, 'saturn03', 'temp', '240'), # 240, 820, 1300
                    (None, 'saturn04', 'temp', '30'), # 30, 830, 860, 910
                    (None, 'am169', 'temp', '260'), # 260, 1100, 1130, 1430
                    (None, 'red26', 'temp', '330'), # 330, 750, 900
                    (None, 'ogi01', 'temp', '80'), # 80 500 1100 5000
                    (None, 'nh10', 'temp', '200'), # 200 1000 6000 7300
                  ]

  botTempFilter = [
                    (None, 'saturn01', 'temp', '1950'), # 0, 740, 1950
                    (None, 'saturn02', 'temp', '3500'), # 100, 600, 1000, 1100, 1600, 2100, 3500
                    (None, 'saturn03', 'temp', '1300'), # 240, 820, 1300
                    (None, 'saturn04', 'temp', '910'), # 30, 830, 860, 910
                    (None, 'saturn05', 'temp', None), # 250
                    (None, 'saturn06', 'temp', None), # 50
                    (None, 'saturn07', 'temp', None), # 60
                    (None, 'am169', 'temp', '1430'), # 260, 1100, 1130, 1430
                    (None, 'red26', 'temp', '900'), # 330, 750, 900
                    (None, 'tansy', 'temp', None), # 840
                    (None, 'yacht', 'temp', None), # 560 590
                    (None, 'sveni', 'temp', None), # 1080
                    (None, 'sandi', 'temp', None), # 790
                    (None, 'ogi01', 'temp', '5000'), # 80 500 1100 5000
                    (None, 'nh10', 'temp', '6000'), # 200 1000 6000 7300
                    (None, 'ncbn1', 'temp', None), # 1200
                    (None, 'marsh', 'temp', None), # 540
                    (None, 'jetta', 'temp', None), # 640
                    (None, 'eliot', 'temp', None), # 1390
                    (None, 'dsdma', 'temp', None), # 820 790
                    (None, 'cbnc3', 'temp', None), # 900 670 650
                    (None, 'yb101', 'temp', None), # 410 320
                    (None, 'woody', 'temp', None), # 240
                    (None, 'tpoin', 'temp', None), # 60
                    (None, 'tnslh', 'temp', None), # 10
                    (None, 'seahs', 'temp', None), # 100
                    (None, 'nbd89', 'temp', None), # 60
                    (None, 'nbd41', 'temp', None), # 60
                    (None, 'nbd29', 'temp', None), # 60
                    (None, 'nb243', 'temp', None), # 60
                    (None, 'lwsck', 'temp', None), # 420
                    (None, 'lonw1', 'temp', None), # 0
                    (None, 'lght6', 'temp', None), # 290
                    (None, 'lght2', 'temp', None), # 440
                    (None, 'grays', 'temp', None), # 160
                    (None, 'coaww', 'temp', None), # 440
                    (None, 'coaof', 'temp', None), # 320 210
                    (None, 'chnkr', 'temp', None), # 0
                    (None, 'chnke', 'temp', None), # 260
                    (None, 'abpoa', 'temp', None), # 410
                    (None, 'bono3', 'temp', None), # 0
                  ]

  # ~distance from mouth
  surfTempStationX = dict( (t[1], stationX[ t[1] ]) for t in surfTempFilter if t[1] in stationX )
  botTempStationX = dict( (t[1], stationX[ t[1] ]) for t in botTempFilter if t[1] in stationX )

  #-------------------------------------------------------------------------------
  # Part II  : get observation data
  #-------------------------------------------------------------------------------

  # a collection for all the data
  dataDir = 'data'
  obsTag = 'obs'
  sc = StationCollection.StationCollection( startTime, endTime, obsTag )

  obsColl = None
  if os.path.isdir( os.path.join(dataDir,obsTag) ) :
    obsColl = StationCollection.StationCollection.loadFromNetCDFCollection( dataDir, obsTag, startTime, endTime )
  if not obsColl or len(obsColl) == 0 :
    obsColl = StationCollection.fetchAvailableObservations( startTime, endTime, obsTag=obsTag )
    obsColl.fetchDischargeData()
    obsColl.fetchTidalData()
    obsColl.fetchCoastalUpwellingData()
    obsColl.saveAsNetCDFCollection( dataDir )
  print ' *** obs data *** '
  for k in obsColl :
    print k

  # add to the main collection
  sc.update( obsColl )

  #-------------------------------------------------------------------------------
  # Part IIIa  : extract model output
  #-------------------------------------------------------------------------------

  offStrings = obsColl.getOfferingStrings()
  if modelDir :
    modColl = None
    if os.path.isdir( os.path.join(dataDir,runTag) ) :
      modColl = StationCollection.StationCollection.loadFromNetCDFCollection( dataDir, runTag, startTime, endTime )
    if not modColl or len(modColl) == 0 :
      modColl = StationCollection.extractFromModelOutputs( runTag, modelDir, offStrings, startTime, endTime,
                            modelCoordSys='spcs', staFile=stationFile )
      modColl.saveAsNetCDFCollection( dataDir )
    print ' *** extracted model data *** '
    for k in modColl :
      print k
    sc.update( modColl )

  #-------------------------------------------------------------------------------
  # Part IIIb  : get hindcast data
  #-------------------------------------------------------------------------------
  for dbTag in hindcastDirs :
    modColl = None
    if os.path.isdir( os.path.join(dataDir,dbTag) ) :
      modColl = StationCollection.StationCollection.loadFromNetCDFCollection( dataDir, dbTag, startTime, endTime )
    if not modColl or len(modColl) == 0 :
      modColl = StationCollection.fetchHindcastFromDatFiles( dbTag, hindcastDirs[dbTag],
                                  offStrings, startTime, endTime, removeBadValues )
      modColl.saveAsNetCDFCollection( dataDir )
    print ' *** hindcast data *** '
    for k in modColl :
      print k
    sc.update( modColl )

  print ' *** StationCollection *** '
  print 'samples',len(sc)
  print 'observation', sc.getObsTag()
  print 'models', sc.getModelTags()

  #-------------------------------------------------------------------------------
  # Part IV : plot results
  #-------------------------------------------------------------------------------
  # remove bad stations from all
  excludeStations = ['tnslh']
  sc = sc.getSubset( exclude=True, station=excludeStations )
  
  # conditions
  pl = Plots(imgDir, sc )
  pl.makeConditionPlot()

  # elevations for water level stations
  data = sc.getSubset(variable='elev',station=elevStations)
  pl = Plots(imgDir, data )
  # taylor
  pl.makeTaylorDiagramsVar()
  # minmax plot
  pl.makeStationExtremaPlot(stationCoords=elevStationX)
  #time series stack plot
  pl.makeTimeSeries( xlim=[startTime,endTime], ylim=ts_ylim, err_ylim=ts_err_ylim )
  pl.makeTimeSeriesStack( elevStationX, xlim=[startTime,endTime], ylim=ts_ylim['elev'], plotError=True )
  #error histogram stack plot
  #pl.makeErrorHistograms( plotMinMax=True, range=eh_lim['elev'], bins=eh_bins )
  pl.makeErrorHistogramStack( elevStationX, plotMinMax=True, range=eh_lim['elev'], bins=eh_bins )
  #error spectral plot
  pl.makeErrorSpectralPlots()
  pl.makeErrorSpectralStack( elevStationX, ylim=esp_ylim )
  # harmonic analysis NOTE slow!!
  pl.makeHarmonicAnalysisErrorPlots( elevStationX )
  
  # all salt
  data = sc.getSubset(variable='salt')
  pl = Plots(imgDir, data )
  pl.makeTimeSeries( xlim=[startTime,endTime], ylim=ts_ylim, err_ylim=ts_err_ylim )
  pl.makeErrorHistograms( plotMinMax=True, range=eh_lim['salt'], bins=eh_bins )

  # surface salt
  data = sc.getSubset(surfSaltFilter)
  pl = Plots(imgDir, data )
  pl.makeStationExtremaPlot(varPrefix='surface',stationCoords=surfSaltStationX)
  pl.makeTaylorDiagramsVar(varPrefix='surface')
  pl.makeTimeSeriesStack( surfSaltStationX, varPrefix='surface', xlim=[startTime,endTime], ylim=ts_ylim['salt'] )
  pl.makeErrorHistogramStack( surfSaltStationX, varPrefix='surface', plotMinMax=True, range=eh_lim['salt'], bins=eh_bins )

  # bottom salt
  data = sc.getSubset(botSaltFilter)
  pl = Plots(imgDir, data )
  pl.makeStationExtremaPlot(varPrefix='bottom',stationCoords=botSaltStationX)
  pl.makeTaylorDiagramsVar(varPrefix='bottom')
  pl.makeTimeSeriesStack( botSaltStationX, varPrefix='bottom', xlim=[startTime,endTime], ylim=ts_ylim['salt'] )
  pl.makeErrorHistogramStack( botSaltStationX, varPrefix='bottom', plotMinMax=True, range=eh_lim['salt'], bins=eh_bins )

  # all temp
  data = sc.getSubset(variable='temp')
  pl = Plots(imgDir, data )
  pl.makeTimeSeries( xlim=[startTime,endTime], ylim=ts_ylim, err_ylim=ts_err_ylim )
  pl.makeErrorHistograms( plotMinMax=True, range=eh_lim['temp'], bins=eh_bins )

  # surface temperature
  data = sc.getSubset(surfTempFilter)
  pl = Plots(imgDir, data )
  pl.makeStationExtremaPlot(varPrefix='surface',stationCoords=surfTempStationX)
  pl.makeTaylorDiagramsVar(varPrefix='surface')
  pl.makeTimeSeriesStack( surfTempStationX, varPrefix='surface', xlim=[startTime,endTime], ylim=ts_ylim['temp'] )
  pl.makeErrorHistogramStack( surfTempStationX, varPrefix='surface', plotMinMax=True, range=eh_lim['temp'], bins=eh_bins )

  # bottom temperature
  data = sc.getSubset(botTempFilter)
  pl = Plots(imgDir, data )
  pl.makeStationExtremaPlot(varPrefix='bottom',stationCoords=botTempStationX)
  pl.makeTaylorDiagramsVar(varPrefix='bottom')
  pl.makeTimeSeriesStack( botTempStationX, varPrefix='bottom', xlim=[startTime,endTime], ylim=ts_ylim['temp'] )
  pl.makeErrorHistogramStack( botTempStationX, varPrefix='bottom', plotMinMax=True, range=eh_lim['temp'], bins=eh_bins )

  ### salt stratification
  ##data = sc.getSubset(botSaltFilter).union( sc.getSubset(surfSaltFilter) ).computeStratification()
  ##for k in sorted(data.getKeys()) :
    ##print k
  ###

#-------------------------------------------------------------------------------
# Parse commandline arguments
#-------------------------------------------------------------------------------
if __name__=='__main__' :

  usage = ('Usage: %prog -r [run id] -s [start date YYYY-MM-DD] -e [end date YYYY-MM-DD]\n')

  parser = OptionParser(usage=usage)
  parser.add_option('-r', '--run', action='store', type='string',
                      dest='runTag', help='Run ID')
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
  parser.add_option('-d', '--modelDirectory', action='store', type='string',
                      dest='modelDir', help='model outputs directory', default=None)

  (options, args) = parser.parse_args()

  runTag = options.runTag
  startStr = options.startStr
  endStr = options.endStr
  modelDir = options.modelDir

  if runTag == None:
      parser.error('Run ID undefined')
  if startStr == None:
      parser.error('Start date undefined')
  if endStr == None:
      parser.error('End date undefined')

  startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
  endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')

  print runTag, startStr, endStr
  print modelDir

  imgDir = createDir(runTag+'/images/')
  
  # hindcasts to consider
  # key (e.g. 'db22') used for plotting, value (e.g. '22') used for the hindcast directory 'xxxx-xx-22'
  hindcastDirs = { }
  if not modelDir :
    hindcastDirs[runTag] = runTag

  skillPlots(runTag, imgDir, modelDir, hindcastDirs, startTime, endTime)
