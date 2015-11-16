#!/usr/bin/env python
"""
A script for generating time series plots.

Examples:

# plot all salt time series from saturn02
makeTimeSeriesPlots -s 2012-5-1 -e 2012-5-10 -o tmp/images */data/stations/saturn02/saturn02.*.salt/*.nc

# plot all salt time series from saturn02, show image do not store
makeTimeSeriesPlots -s 2012-5-1 -e 2012-5-10 */data/stations/saturn02/saturn02.*.salt/*.nc

Tuomas Karna 2013-01-07
"""

import os
import datetime
import sys
from optparse import OptionParser

import numpy as np
from crane import plt

from crane.data import dirTreeManager as dtm
from crane.data import stationCollection
from crane.files import csvStationFile
from crane.data import dataContainer
from crane.data import statistics

from crane.plotting import plot
from crane.plotting import plotBase
from crane.plotting import timeSeriesPlot


#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
def uniqueList(seq):
  "Return list with all duplicates removed. Preserves ordering."""
  seen = set()
  return [x for x in seq if x not in seen and not seen.add(x)]

def generateTSPlot(startTime, endTime, dcs, imgDir=None, ylim={},
                   obsTag='obs', **kwargs) :
  """Plots time series stack plot"""

  doMurphyScore = True

  # a collection for all the data
  dataDir = 'data'
  sc = stationCollection.StationCollection( startTime, endTime, obsTag )

  # add data from the files
  for dc in dcs :
    if dc.data.shape[1]>1 : # more than one component
      for i in range(dc.data.shape[1]) :
        # add each component as separate time series
        dci = dc.extractField(i)
        # change variable to component name: 'hvel' -> 'u','v'
        dci.setMetaData('variable',dci.fieldNames[0])
        sc.addSample(dci)
    else :
      sc.addSample(dc)

  # plot
  # master container to hold all the axes, x axis will be shared
  canvas = plotBase.stackPlotBase(rightPad=1.25)

  modelColors = plot.makeColorsForModels(sc)
  comKeys = sc.getComparableKeys(requireObs=False, requireMod=False)
  for c in comKeys :
    if c[0][2] is None :
      print c
      c[0] = (c[0][0],c[0][1],'0')
  # fancy sort: sort by 'entry' tuple with msldepth converted to float
  comKeys.sort(key=lambda tup: (tup[0][0],tup[0][1],float(tup[0][2])))
  for entry, obsKey, modKeys in comKeys :
    station,var,msldepth = entry
    if obsKey :
      o = sc.getSample( **obsKey )
      if len(o.time) < 3 :
        print 'observation time series too short:',obsKey
        continue
    else :
      o = None

    title = ' '.join((plot.VARS.get(var,var),'['+plot.UNITS.get(var,'-')+']'))
    xlabel = startTime.year if startTime.year==endTime.year else 'Date'
    dia = timeSeriesPlot.timeSeriesPlotDC2( title=title, ylabel=station.upper()+' '+msldepth, ylim=ylim.get(var,None), xlabel=xlabel )
    canvas.addPlot( dia, tag=plot.VARS.get(var,var)+station+msldepth )
    if obsKey :
      dia.addSample( o, color='r', label=obsKey['tag'] )
    for modKey in modKeys :
      m = sc.getSample( **modKey )
      l = modKey['tag'].split('-')[0]
      if doMurphyScore:
          o2, m2 = o.alignTimes(m)
          murphy = statistics.murphySkillScore(o2.data.ravel(), m2.data.ravel())
          l += ' MS={ms:.2f}'.format(ms=murphy)
      dia.addSample(m, color=modelColors[modKey['tag']], label=l)
    dia.showLegend()

  if 'title' in kwargs:
      canvas.addTitle(kwargs.pop('title'))
  if 'xaxis_tight' in kwargs and kwargs['xaxis_tight']==True:
      
      xLim = [canvas.axGrid[0].dataLim.xmin, canvas.axGrid[0].dataLim.xmax]
      for ax in canvas.axGrid[1:]:
          a, b = ax.dataLim.xmin, ax.dataLim.xmax
          xLim = [min(xLim[0], a), max(xLim[1], b)]
      canvas.axGrid[0].set_xlim(xLim)
  if imgDir :
    # ----- Save file -----
    sT = str(sc.startTime.date())
    eT = str(sc.endTime.date())
    imgDir = plotBase.createDirectory(imgDir)
    varList = uniqueList([ tup[0][1] for tup in comKeys ])
    #depSet = list(depSet)
    tagStr = '-'.join(uniqueList([ dc.getMetaData('tag') for dc in dcs ]))
    stationList = []
    for dc in dcs :
      meta = dc.getMetaData()
      entry = meta.get('location')
      if 'msldepth' in meta : entry+='-'+meta['msldepth']
      stationList.append( entry )
    stationStr = '-'.join(uniqueList(stationList))
    varStr = '-'.join(varList)
    fname = '_'.join(['ts',varStr,stationStr,sT,eT])
    fn = kwargs.get('filename', fname)
    saveFigure(imgDir, fn, 'png', verbose=True, bbox_tight=True)
    plt.close()
  else :
    # ----- show plot -----
    plt.show()

#-------------------------------------------------------------------------------
# Main routine
#-------------------------------------------------------------------------------
def makeTSPlotForStationFile(runTags, stationFile, startTime, endTime,
                              imgDir=None, ylim={}):
  # load data
  dcs = []
  stationsToExtract = csvStationFile.csvStationFileWithDepth()
  stationsToExtract.readFromFile(stationFile)
  for runTag in runTags :
    for key in stationsToExtract.getTuples() :
      loc,x,y,z,zType,var = key
      try :
        if var[:3] == 'sil' :
          dc = dtm.getDataContainer(tag=runTag,dataType='sil',
                                    location=loc, variable=var,
                                    startTime=startTime,endTime=endTime)
        elif loc == 'plume' :
          dc = dtm.getDataContainer(tag=runTag,dataType='plumemetrics',
                                    location=loc, variable=var,
                                    startTime=startTime,endTime=endTime)
        else :
          msldepth = str(int(round(abs(z)*100)))
          dc = dtm.getDataContainer(tag=runTag,dataType='timeseries',
                                    location=loc,
                                    variable=var,msldepth=msldepth,
                                    startTime=startTime,endTime=endTime)
        dcs.append(dc)
      except Exception as e:
        print 'reading failed'
        print e

  if dcs :
    generateTSPlot(startTime, endTime, dcs, imgDir, ylim)

#-------------------------------------------------------------------------------
# Parse commandline arguments
#-------------------------------------------------------------------------------
def parseCommandLine() :

  usage = ('Usage: %prog -s [start date YYYY-MM-DD] -e [end date YYYY-MM-DD] -o [path] -t [stationFile] runTag1 runTag2 ...\n')

  parser = OptionParser(usage=usage)
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
  parser.add_option('-o', '--imageDirectory', action='store', type='string',
                      dest='imgDir', help='(optional) directory where generated images are stored. If not specified, shows the image instead.')
  parser.add_option('-y', '--ylim', action='store', type='string',
                      dest='ylimStr', help='Custom limits for y axis, a string like salt:0:30,temp:4:12')
  parser.add_option('-t', '--csvStationFile', action='store', type='string',
                      dest='csvStationFile',
                      help='file that defines station coordinates and\
                      variables for time series to plot')

  (options, args) = parser.parse_args()
  
  runTags = args

  startStr = options.startStr
  endStr = options.endStr
  imgDir = options.imgDir
  ylimStr = options.ylimStr
  csvStationFile = options.csvStationFile

  if not runTags :
    parser.print_help()
    parser.error('RunTags are missing')
  if not csvStationFile :
    parser.print_help()
    parser.error('csvStationFile must be given')
  if startStr is None:
    parser.print_help()
    parser.error('Start date undefined')
  if endStr is None:
    parser.print_help()
    parser.error('End date undefined')

  startTime = plotBase.parseTimeStr( startStr )
  endTime = plotBase.parseTimeStr( endStr )

  ylim = {}
  if ylimStr :
    for entry in ylimStr.split(',') :
      var,vmin,vmax = entry.split(':')
      ylim[var] = [float(vmin),float(vmax)]

  print 'Parsed options:'
  print ' - time range:',str(startTime),'->', str(endTime)
  if imgDir :
    print ' - output dir',imgDir
  else :
    print ' - show image'
  if ylim :
    print ' - using y limits',ylim
  if csvStationFile :
    print ' - stations read from',csvStationFile

  makeTSPlotForStationFile(runTags, csvStationFile, startTime, endTime,
                           imgDir, ylim)

if __name__=='__main__' :
  parseCommandLine()
