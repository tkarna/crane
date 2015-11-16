#!/usr/bin/env python
"""
Main script for generating transect images.

Examples:

makeSlabPlots -o run29/images/slab/ -r run29 -b [326000,352000,284000,302000] -t nc_tran_stations.sta -a nchannel_fine.bp run29/data/slab/slab_salt_slev1_2012-05-01_2012-05-11.nc

Tuomas Karna 2012-12-02
"""
import numpy as np
import datetime
import sys
from optparse import OptionParser

from crane.data import meshContainer
from crane.data import timeArray
from crane.data import collection
from crane.files import stationFile
from crane.files import buildPoints
from crane.data.selfeGridUtils import readAnyMeshFile

import traceback

from crane import plt
from crane import matplotlib
from crane.physicalVariableDefs import VARS
from crane.physicalVariableDefs import UNITS
from crane.plotting import slabPlot
from crane.utility import createDirectory
from crane.utility import saveFigure
from crane.utility import parseTimeStr

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
    traceback.print_exc(file=sys.stdout)
    raise e

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
  except Exception as e:
    print ' Plot generation failed'
    raise e

def processFrame(dcs, time, logScaleVars, aspect, clim, diffClim, cmap, bBox,
                 transectFile, stationFileObj, bathMC, isobaths, diff, imgDir, fPrefix, filetype,
                 maxPlotSize=6.0):
  """Plots only the first time step of the meshContainers."""
  it = 0

  # initialize plot
  height = maxPlotSize if aspect < 1.0 else maxPlotSize/aspect
  width = maxPlotSize if aspect > 1.0 else maxPlotSize*aspect
  #height = minsize if aspect > 1.0 else minsize/aspect
  #width = minsize*aspect if aspect > 1.0 else minsize
  width += 1.3 # make room for colorbar & labels
  dia = slabPlot.stackSlabPlotDC(figwidth=width, plotheight=height)

  varList = []
  for i,dc in enumerate(dcs) :
    dateStr = dc.time.getDatetime(it).strftime('%Y-%m-%d %H:%M')
    meta = dc.getMetaData()
    name = meta.get('location','')
    var = dc.fieldNames[0]
    logScale = var in logScaleVars
    climIsLog = logScale
    tag = str(dc.getMetaData('tag',suppressError=True))
    pltTag = tag+name+dc.getMetaData('variable')+'-'+str(i)
    titleStr = tag+' '+dateStr+' (PST)'

    dia.addPlot(pltTag, clabel=VARS.get(var,var), unit=UNITS.get(var,'-'), bbox=bBox)
    # add bathymetry contours (if any)
    if bathMC is not None and len(isobaths) > 0:
      dia.addSample(pltTag, bathMC, 0, plotType='contour', levels=isobaths, colors='k', bbox=bBox, zorder=1, draw_cbar=False)
    dia.addSample(pltTag, dc, it,
                  clabel=VARS.get(var,var), unit=UNITS.get(var,'-'),
                  clim=clim.get(var, None), cmap=cmap, bbox=bBox,
                  logScale=logScale, climIsLog=climIsLog, zorder=0)

    dia.addTitle(titleStr,tag=pltTag)

    varList.append(var)

  # plot differences between slabs
  if diff :
    # assuming sequential dataContainers are comparable (salt,salt, temp,temp)
    nPairs = int(np.floor(len(dcs)/2))
    pairs = []
    for i in range(nPairs):
      pairs.append((dcs[2*i], dcs[2*i+1]))
    for joe, amy in pairs:
      diffDC = joe.copy()
      diffDC.data = joe.data - amy.data
      meta = diffDC.getMetaData()
      name = meta.get('location','')
      var = diffDC.fieldNames[0]
      tag = '(' + str(joe.getMetaData('tag',suppressError=True)) + '-' +\
            str(amy.getMetaData('tag',suppressError=True)) + ')'
      pltTag = tag+name+var
      titleStr = tag+' '+dateStr+' (PST)'
      # add bathymetry contours (if any)
      if bathMC is not None and len(isobaths) > 0:
        dia.addSample(pltTag, bathMC, 0, plotType='contour', levels=isobaths, colors='k',
                      clabel=VARS.get('bath','Bathymetry'), unit=UNITS.get('bath','m'), bbox=bBox, zorder=1)
      dia.addSample(pltTag, diffDC, it,
                    clabel='diff '+VARS.get(var,var), unit=UNITS.get(var,'-'),
                    clim=diffClim.get(var, None), cmap=cmap, bbox=bBox,
                    logScale=logScale, climIsLog=climIsLog)
      dia.addTitle(titleStr,tag=pltTag)
      varList.append('diff_'+var)

  # add transect markers (if any)
  if transectFile is not None:
    for bp in transectFile :
        dia.addTransectMarker('all', bp.getX(), bp.getY(),
                              color='w', linewidth=2.0)
        dia.addTransectMarker('all', bp.getX(), bp.getY(),
                              color='k', linewidth=1.0)

  # add station markers (if any)
  if stationFileObj is not None:
    for sta in stationFileObj :
        xSta,ySta = stationFileObj.getLocation( sta )
        dia.addStationMarker('all', xSta, ySta, sta.replace('saturn','sat'),
                             printLabel=True, color='k')

  # save to disk
  dateStr = dateStr.replace(' ','_').replace(':','-')
  varStr = '-'.join(collection.uniqueList(varList))
  file = '_'.join([fPrefix,name,varStr,dateStr])
  saveFigure(imgDir, file, filetype, verbose=True, dpi=200, bbox_tight=True)
  plt.close(dia.fig)

#-------------------------------------------------------------------------------
# Main routine
#-------------------------------------------------------------------------------
def makeSlabPlots(netCDFFiles, imgDir, runTag=None, startTime=None,
                  endTime=None, skip=1, stationFilePath=None,
                  transectFilePath=None, bathFilePath=None, isobaths=[], userClim={}, cmapStr=None, bBox=None,
                  diff=False, userDiffClim={}, num_threads=1, maxPlotSize=6.0):

  imgDir = createDirectory(imgDir)
  fPrefix = 'slab'
  filetype = 'png'
  cmap = plt.get_cmap(cmapStr)

  # logScale flag
  logScaleVars = ['kine','vdff','bed_stress']
  # color range
  clim = {'salt':[0,35],
          'temp':[5,20],
          'kine':[-6,-0.5],
          'vdff':[-6,-0.5],
          'hvel':[-3.0,3.0]
          }
  clim.update(userClim)

  # color range for differences
  diffClim = {'salt':[-5,5],
              'temp':[-1,1],
              'hvel':[-1,1],
              }
  diffClim.update( userDiffClim )

  # number of contour lines
  nCLines = {'salt':71,
             'temp':16,
             'kine':23,
             'vdff':23,
             'hvel':31,
             }

  dcs = []
  for fn in netCDFFiles :
    if startTime and endTime and startTime != endTime :
      dc = meshContainer.meshContainer.loadFromNetCDF(fn, startTime, endTime)
    else:
      dc = meshContainer.meshContainer.loadFromNetCDF(fn)
    if bBox :
      dc = dc.cropGrid(bBox)
    dcs.append(dc)

  stationFileObj = transectFile = None
  if stationFilePath:
    stationFileObj = stationFile.StationFile( stationFilePath )
    stationFileObj.readFileFromDisk()

  if transectFilePath :
    transectFile = []
    for trf in transectFilePath :
      bp = buildPoints.BuildPoint()
      bp.readFileFromDisk(trf)
      transectFile.append(bp)
  
  if bathFilePath:
    bathMC = readAnyMeshFile(bathFilePath)
  else:
    bathMC = None

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

  # append each component separately (if more than one)
  dcs_comp = []
  for iComp in range(3):
    for dc in dcs:
      if iComp < dc.data.shape[1]:
        dcs_comp.append(dc.extractFields(iComp))
  dcs = dcs_comp

  aspect = 1.0 # figure aspect ratio
  if bBox :
    aspect = float(bBox[1]-bBox[0])/float(bBox[3]-bBox[2])
  else :
    aspect = (dcs[0].x.max()-dcs[0].x.min())/(dcs[0].y.max()-dcs[0].y.min())
    aspect = max(min(aspect,20.),1./20.)

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
    args = [dcs_single, time, logScaleVars, aspect, clim, diffClim, cmap, bBox,
            transectFile, stationFileObj, bathMC, isobaths, diff, imgDir, fPrefix, filetype,
            maxPlotSize]
    tasks.append((function, args))
  _runTasksInQueue(num_threads, tasks)
  #print 'duration', timeMod.time()-t0


#-------------------------------------------------------------------------------
# Parse commandline arguments
#-------------------------------------------------------------------------------
def parseCommandLine() :

  usage = ('Usage: %prog -s [start date YYYY-MM-DD] -e [end date YYYY-MM-DD] -o [path] -t [stationFile] slabFile\n')

  parser = OptionParser(usage=usage)
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startTime', help='Date to start processing (Default: start of data)')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endTime', help='Date to end processing (Default: end of data)')
  parser.add_option('-o', '--imageDirectory', action='store', type='string',
                      dest='imgDir', help='directory where generated images are stored')
  parser.add_option('-t', '--stationFile', action='store', type='string',
                      dest='stationFile', help='file that defines station coordinates to append in the plots')
  parser.add_option('-r', '--runID', action='store', type='string',
                      dest='runTag', help='Run ID, added in the plot title, e.g. run29 or "selfe k-e" (Default: empty string)')
  parser.add_option('-k', '--skip', action='store', type='int',
                      dest='skip', help='Generate every skip -th time step. (Default: %default)',default=1)
  parser.add_option('-c', '--clim', action='store', type='string',
                      dest='climStr', help='Custom limits for color bar, a string like salt:0:30,kine:-6:-2')
  parser.add_option('-b', '--boundingBox', action='store', type='string',
                      dest='bBox', help='Limits for x and y, [xmin,xmax,ymin,ymax]')
  parser.add_option('-a', '--transect', action='store', type='string',
                      dest='trFiles', help='comma separated list of transect bp files to draw in the plot')
  parser.add_option('-M', '--colormap', action='store', type='string',
                      dest='cmapStr', help='name of matplotlib colormap to use')
  parser.add_option('-D', '--diff', action='store_true',
                      dest='diff', help='plot difference of two transects (Default: %default)', default=False)
  parser.add_option('-d', '--diffClim', action='store', type='string',
                      dest='diffClimStr', help='Custom limits for color bar if plotting differences, a string like salt:-1:1,kine:-0.5:-0.5')
  parser.add_option('-j', '--numThreads', action='store', type='int',
                      dest='num_threads', help='Number of concurrent threads to use (default %default).',default=1)
  parser.add_option('', '--maxPlotSize', action='store', type='float',
                      dest='maxPlotSize', help='Max size of the longest edge of individual plot in inches (default %default).',
                      default=8)
  parser.add_option('', '--fontSize', action='store', type='float',
                      dest='fontSize', help='Set the font size for all labels (default %default).',
                      default=12)
  parser.add_option('', '--bathymetry', action='store', type='string',
                      dest='bathMeshFile', help='A mesh file that contains bathymetry data (for plotting isobath lines)')
  parser.add_option('', '--isobaths', action='store', type='string',
                      dest='isobathStr', help='list of isobaths to plot, e.g. \"50,80,100,400\". Units meters below datum.')

  (options, args) = parser.parse_args()

  startTime = options.startTime
  endTime = options.endTime
  imgDir = options.imgDir
  stationFile = options.stationFile
  runTag = options.runTag.strip('\"') if options.runTag else None
  skip = options.skip
  climStr = options.climStr
  bBox = options.bBox
  trFiles = options.trFiles.split(',') if options.trFiles else None
  cmapStr = options.cmapStr
  diff = options.diff
  diffClimStr = options.diffClimStr
  num_threads = options.num_threads
  maxPlotSize = options.maxPlotSize
  matplotlib.rcParams['font.size'] = options.fontSize
  bathMeshFile = options.bathMeshFile
  isobathStr = options.isobathStr

  if imgDir is None:
    parser.print_help()
    parser.error('imgDir undefined')
  if len(args) < 1 :
    parser.print_help()
    parser.error('slabFile missing')
  if len(args) != 2 and diff :
    parser.print_help()
    parser.error('can only take difference of two slabs')
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

  if bBox :
    bBox = [ float(f) for f in bBox.strip('[').strip(']').split(',') ]
    
  isobaths = []
  if isobathStr:
    isobaths = [float(s) for s in isobathStr.split(',')]

  print 'Parsed options:'
  if runTag :
    print ' - runID ', runTag
  if startTime :
    print ' - time range:',str(startTime),'->', str(endTime)
  print ' - output dir',imgDir
  print ' - slab files'
  for fn in netCDFFiles :
    print '    ', fn
  if stationFile :
    print ' - using stationFile',stationFile
  if clim :
    print ' - using color limits',clim
  if trFiles :
    print ' - plotting trasects', trFiles
  if cmapStr :
    print ' - using color map',cmapStr
  if bBox :
    print ' - using bounding box',bBox
  if diff :
    print ' - plotting difference of transects'
  if diffClim :
    print ' - using difference color limits',diffClim
  print ' - number of parallel threads', num_threads
  print ' - max plot size (in)', maxPlotSize
  print ' - font size (pt)', options.fontSize


  makeSlabPlots(netCDFFiles, imgDir, runTag, startTime, endTime, skip,
                stationFile, trFiles, bathMeshFile, isobaths, clim, cmapStr, bBox, diff, diffClim,
                num_threads, maxPlotSize, )

if __name__=='__main__' :
  parseCommandLine()
