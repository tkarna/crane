import numpy as np
import datetime
import sys
from optparse import OptionParser

from crane.data import dataContainer
from crane.data import timeArray
from crane.data import collection
from crane.files import csvStationFile

from crane import matplotlib
from crane import plt
from crane.physicalVariableDefs import VARS
from crane.physicalVariableDefs import UNITS
from crane.utility import saveFigure
from crane.utility import createDirectory
from crane.plotting import transectPlot

import multiprocessing
# NOTE this must not be a local function in runTasksInQueue
def _launch_job( task ) :
  """Excecutes a task"""
  function, args = task
  try :
    function( *args )
  except KeyboardInterrupt as e:
    raise e
  except Exception as e :
    print e

def _runTasksInQueue( num_threads, tasks ) :
  """Generic routine for processing tasks in parallel.

  Tasks are defined as a (function,args) tuple, whick will is executed as
  function(*args).

  num_threads - (int) number of threads to use
  tasks - list of (function,args)
  """
  pool = multiprocessing.Pool(num_threads)
  p = pool.map_async(_launch_job, tasks, chunksize=1)
  timeLimit = 24*3600
  try:
    result = p.get(timeLimit)
  except KeyboardInterrupt:
    print ' Ctrl-c received! Killing...'

def processFrame(dcs, logScaleVars, xlim, ylim, nCLines, defaultNCLines, clim, diffClim,
                 cmapDict, stationFile, diff, imgDir, fPrefix, filetype,
                 maxPlotSize=10.0):
    # make empty plot
    dia = transectPlot.stackTransectPlotDC(figwidth=maxPlotSize,
                              plotheight=2.0/10.0*maxPlotSize)
    it = 0
    # round time to nearest minute
    varList = []
    for dc in dcs :
      # add transect
      transectName = str(dc.getMetaData('location',suppressError=True))
      tag = str(dc.getMetaData('tag',suppressError=True))
      var = dc.getMetaData('variable')
      logScale = var in logScaleVars
      climIsLog = logScale
      pltTag = tag+transectName+dc.getMetaData('variable')
      if isinstance(cmapDict,dict) :
        cmap = plt.get_cmap(cmapDict.get(var,'Spectral_r'))
      else :
        cmap = plt.get_cmap(cmapDict)
      dateEpoch = round(dc.time.asEpoch()[it]/60)*60
      dateStr = timeArray.epochToDatetime(dateEpoch).strftime('%Y-%m-%d %H:%M')
      titleStr = tag+' '+transectName+' '+dateStr+' (PST)'
      if it==0:
        dia.addSample(pltTag, dc,it, N=nCLines.get(var,defaultNCLines),
                      clim=clim.get(var,None),
                      clabel=VARS.get(var,var),unit=UNITS.get(var,'-'),
                      logScale=logScale, climIsLog=climIsLog, cmap=cmap,
                      ylim=ylim, xlim=xlim)
        dia.showColorBar(tag=pltTag)
        dia.addTitle(titleStr,tag=pltTag)
      else:
        # the mesh changes in time, must render entire figure again
        # remove all previous contours TODO is there a higher level version?
        for p in dia.plots[pltTag].cax.collections:
            p.remove()
        dia.addSample(pltTag, dc,it, N=nCLines.get(var,defaultNCLines),
                       clim=clim.get(var,None),
                       clabel=VARS.get(var,var),unit=UNITS.get(var,'-'),
                       logScale=logScale, climIsLog=climIsLog, cmap=cmap,
                       ylim=ylim, xlim=xlim)
        dia.updateColorBar(tag=pltTag)
        dia.updateTitle(titleStr,tag=pltTag)
      varList.append(var)

    # plot difference of two transects
    if diff :
      diffDC = dcs[0].copy()
      diffDC.data = dcs[0].data - dcs[1].data
      tag = '('+str(dcs[0].getMetaData('tag',suppressError=True))+'-'+str(dcs[1].getMetaData('tag',suppressError=True))+')'
      pltTag = tag+transectName+var
      titleStr = tag+' '+transectName+' '+dateStr+' (PST)'
      if it==0:
        dia.addSample(pltTag, diffDC, it, N=nCLines.get(var,20),
                      clim=diffClim.get(var,None), clabel=VARS.get(var,var),
                      unit=UNITS.get(var,'-'), cmap='RdBu_r',
                      ylim=ylim, xlim=xlim)
        dia.addTitle(titleStr,tag=pltTag)
        dia.showColorBar(tag=pltTag)
      else:
        for p in dia.plots[pltTag].cax.collections:
            p.remove()
        # remove all previous contours TODO is there a higher level version?
        dia.addSample(pltTag, dc,it, N=nCLines.get(var,defaultNCLines),
                       clim=clim.get(var,None),
                       clabel=VARS.get(var,var),unit=UNITS.get(var,'-'),
                       logScale=logScale, climIsLog=climIsLog, cmap=cmap,
                       ylim=ylim, xlim=xlim)
        dia.updateColorBar(tag=pltTag)
        dia.updateTitle(titleStr,tag=pltTag)
      varList.append(var)

    # add station markers (if any)
    if it==0:
      if stationFile :
        for sta in stationFile :
          xSta,ySta = stationFile.getLocation( sta[0] )
          dia.addStationMarker( 'all', xSta, ySta,
                                label=sta[0].replace('saturn','sat'), color='k',
                                linewidth=1.5, linestyle='dashed' )

    # save to disk
    dateStr = dateStr.replace(' ','_').replace(':','-')
    file = '_'.join([fPrefix,transectName,'-'.join(collection.uniqueList(varList)),dateStr])
    saveFigure(imgDir,file,filetype,verbose=True, dpi=100, bbox_tight=True)
    plt.close(dia.fig)

#-------------------------------------------------------------------------------
# Main routine
#-------------------------------------------------------------------------------
def makeTransects(netCDFFiles, imgDir, startTime=None, endTime=None, skip=1,
                  stationFilePath=None, userClim={}, cmapDict=None, ylim=None,
                  xlim=None, diff=False, userDiffClim={}, num_threads=1,
                  maxPlotSize=10.0) :

  imgDir = createDirectory(imgDir)
  fPrefix = 'trans'
  filetype = 'png'
  
  # logScale flag
  logScaleVars = ['kine','vdff','tdff']
  # color range
  clim = { 'salt':[0,35],
           'temp':[5,20],
           'kine':[-6,-0.5],
           'vdff':[-6,-0.5],
           'tdff':[-6,-0.5],
           'hvel':[-3.0,3.0],
           'sed':[0,100],
           'sed_1':[0,0.01],
           'sed_2':[0,0.3],
           'sed_3':[0,0.4],
           'sed_4':[0,0.5],
           'NO3':[20,40],
           'NH4':[0,1],
           'PHYM':[0,5],
           'PHYF':[0,5],
           'SZOO':[0,1],
           'BZOO':[0,1],
           'DETN':[0,3],
           'DETC':[0,9],
           'BACT':[0,1],
           'DON':[0,3],
           'DOC':[0,9],
           'CHLM':[0,12],
           'CHLF':[0,12],
           'OXY':[0,300],
           'Diag':[0,1],
           'nem':[-10,10],
           'nem_DAVG':[-0.5,0.5],
           'nem_DI':[-5,5],
           'totalN':[0,40],
           'bnth_1':[0,25],
           'bnth_2':[0,15],
           'bnth_3':[0,150],
           'bnth_4':[0,300],
         }
  clim.update( userClim )
  
  # color range for differences 
  diffClim = { 'salt':[-5,5], 
             }
  diffClim.update( userDiffClim )

  # number of contour lines
  nCLines = { 'salt': 25,#71,
              'temp':16,
              'kine':23,
              'vdff':23,
              'hvel':31,
              'sed':11,
              'sed_1':11,
              'sed_2':11,
              'sed_3':11,
              'sed_4':11,
            }
  defaultNCLines = 30

  dcs = []
  for fn in netCDFFiles :
    dc = dataContainer.dataContainer.loadFromNetCDF(fn)
    if ylim is None:
        ylim = [np.nanmin(dc.z), np.nanmax(dc.z)]
    dcs.append( dc )

  if stationFilePath :
    stationFile = csvStationFile.csvStationFile()
    stationFile.readFromFile( stationFilePath )
  else :
    stationFile = None

  if startTime or endTime :
    if startTime == endTime :
      # interpolate to given time stamp
      t0 = timeArray.datetimeToEpochTime(startTime)
      newTime = timeArray.timeArray( np.array([t0]), 'epoch' )
      for i in range(len(dcs)) :
        dcs[i] = dcs[i].interpolateInTime( newTime )
    else :
      # restrict time window based on the given start/end time
      if not startTime :
        startTime = dcs[0].time.getDatetime(0)
      if not endTime :
        endTime = dcs[0].time.getDatetime(-1)
      for i in range(len(dcs)) :
        dcs[i] = dcs[i].timeWindow( startTime, endTime, includeEnd=True )

  # check that all time arrays match
  time = dcs[0].time.asEpoch().array
  for dc in dcs[1:] :
    if not np.array_equal(dcs[0].time.asEpoch().array, time) :
      raise Exception('time steps must match')

  # add all plotting tasks in queue and excecute with threads
  #import time as timeMod
  tasks = []
  #t0 = timeMod.time()
  for it in range(0,len(time),skip):
    dcs_single = []
    for dc in dcs:
        dcs_single.append(dc.subsample([it]))
    function = processFrame
    args = [dcs_single, logScaleVars, xlim, ylim, nCLines, defaultNCLines,
            clim, diffClim, cmapDict, stationFile, diff, imgDir, fPrefix, filetype,
            maxPlotSize]
    tasks.append((function, args))
  _runTasksInQueue(num_threads, tasks)
  #print 'duration', timeMod.time()-t0

#-------------------------------------------------------------------------------
# Parse commandline arguments
#-------------------------------------------------------------------------------

def parseCommandLine() :
  usage = ('Usage: %prog -s [start date YYYY-MM-DD] -e [end date YYYY-MM-DD] -o [path] -t [stationFile] transectFile1 transectFile2 ...\n')

  parser = OptionParser(usage=usage)
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startTime', help='Date to start processing (Default: start of data)')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endTime', help='Date to end processing (Default: end of data)')
  parser.add_option('-o', '--imageDirectory', action='store', type='string',
                      dest='imgDir', help='directory where generated images are stored')
  parser.add_option('-t', '--stationFile', action='store', type='string',
                      dest='stationFile', help='file that defines station coordinates to append in the plots')
  parser.add_option('-k', '--skip', action='store', type='int',
                      dest='skip', help='Generate every skip -th time step. (Default: %default)',default=1)
  parser.add_option('-c', '--clim', action='store', type='string',
                      dest='climStr', help='Custom limits for color bar, a string like salt:0:30,kine:-6:-2')
  parser.add_option('-y', '--ylim', action='store', type='string',
                      dest='ylimStr', help='Custom limits for the transect depth, a string like -400,5')
  parser.add_option('-x', '--xlim', action='store', type='string',
                      dest='xlimStr', help='Custom limits for the transect length, a string like 0,40 (in km)')
  parser.add_option('-M', '--colormap', action='store', type='string',
                      dest='cmapStr', help='name of matplotlib colormap to use. Can be specified for each variable: hvel:RdBu_r,salt:Spectral_r',
                      default='jet')
  parser.add_option('-D', '--diff', action='store_true', 
                      dest='diff', help='plot difference of two transects (Default: %default)', default=False)
  parser.add_option('-d', '--diffClim', action='store', type='string',
                      dest='diffClimStr', help='Custom limits for color bar if plotting differences, a string like salt:-1:1,kine:-0.5:-0.5')
  parser.add_option('-j', '--numThreads', action='store', type='int',
                      dest='num_threads', help='Number of concurrent threads to use (default %default).',default=1)
  parser.add_option('', '--maxPlotSize', action='store', type='float',
                      dest='maxPlotSize', help='The width of each plot in inches (default %default).',
                      default=10)
  parser.add_option('', '--fontSize', action='store', type='float',
                      dest='fontSize', help='Set the font size for all labels (default %default).',
                      default=12)

  (options, args) = parser.parse_args()

  startTime = options.startTime
  endTime = options.endTime
  imgDir = options.imgDir
  stationFile = options.stationFile
  skip = options.skip
  climStr = options.climStr
  ylimStr = options.ylimStr
  xlimStr = options.xlimStr
  cmapStr = options.cmapStr
  diff = options.diff
  diffClimStr = options.diffClimStr
  num_threads = options.num_threads
  maxPlotSize = options.maxPlotSize
  matplotlib.rcParams['font.size'] = options.fontSize

  if imgDir == None:
    parser.print_help()
    parser.error('imgDir undefined')
  if len(args) < 1 :
    parser.print_help()
    parser.error('transectFile missing')
  if len(args) != 2 and diff :
    parser.print_help()
    parser.error('can only take difference of two transects')
  if diffClimStr and not diff :
    parser.print_help()
    parser.error('can not specify colormap for differences. \nRequires -D to make difference transect plot.') 

  netCDFFiles = args

  if startTime :
    startTime = parseTimeStr( startTime )
  if endTime :
    endTime = parseTimeStr( endTime )

  clim = {}
  if climStr :
    for entry in climStr.split(',') :
      var,vmin,vmax = entry.split(':')
      clim[var] = [float(vmin),float(vmax)]

  diffClim = {}
  if diffClimStr :
    for entry in diffClimStr.split(',') :
      var,vmin,vmax = entry.split(':')
      diffClim[var] = [float(vmin),float(vmax)]

  ylim = [ float(f) for f in ylimStr.split(',') ] if ylimStr else None
  xlim = [ float(f) for f in xlimStr.split(',') ] if xlimStr else None

  if cmapStr :
    if len(cmapStr.split(','))>1 :
      cmap = {}
      for entry in cmapStr.split(',') :
        var,mapname = entry.split(':')
        cmap[var] = mapname
    else :
      cmap = str(cmapStr)

  print 'Parsed options:'
  if startTime :
    print ' - time range:',str(startTime),'->', str(endTime)
  print ' - output dir',imgDir
  print ' - transect files'
  for fn in netCDFFiles :
    print '   ',fn
  if stationFile :
    print ' - using stationFile',stationFile
  if clim :
    print ' - using color limits',clim
  if ylim :
    print ' - using depth limits',ylim
  if xlim :
    print ' - using x limits',xlim
  if cmapStr :
    print ' - using color map',cmap
  if diff :
    print ' - plotting difference of transects'
  if diffClim :
    print ' - using difference color limits',diffClim
  print ' - number of parallel threads', num_threads
  print ' - max plot size (in)', maxPlotSize
  print ' - font size (pt)', options.fontSize 

  makeTransects(netCDFFiles, imgDir, startTime, endTime, skip, stationFile,
                clim, cmap, ylim, xlim, diff, diffClim, num_threads,
                maxPlotSize)

if __name__=='__main__' :
  parseCommandLine()
