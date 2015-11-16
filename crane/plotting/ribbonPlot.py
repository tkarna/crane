#!/usr/bin/env python
"""
Plots error in time series versus tidal range (x) and river discharge (y).


Tuomas Karna 2013-04-29
"""
import os
import numpy as np
import datetime
import sys
from optparse import OptionParser
import glob
from scipy.interpolate import interp1d, griddata, splev, splrep
import matplotlib.gridspec as gridspec

from crane import matplotlib
from crane import plt

from crane.data import dirTreeManager
from crane.data import dataContainer
from crane.data import timeSeriesFilters

from crane.plotting import trackPlot
from crane.physicalVariableDefs import VARS
from crane.physicalVariableDefs import UNITS
from crane.utility import createDirectory
from crane.utility import saveFigure
from crane.utility import parseTimeStr

fontsize=20
matplotlib.rcParams['font.size']=fontsize

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def cmap_discretize(cmap, N):
  """Return a discrete colormap from the continuous colormap cmap.

      cmap: colormap instance, eg. cm.jet.
      N: number of colors.

  Example
      x = resize(arange(100), (5,100))
      djet = cmap_discretize(cm.jet, 5)
      imshow(x, cmap=djet)
  """

  if type(cmap) == str:
    cmap = plt.get_cmap(cmap)
  colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
  colors_rgba = cmap(colors_i)
  indices = np.linspace(0, 1., N+1)
  cdict = {}
  for ki,key in enumerate(('red','green','blue')):
    cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
  # Return colormap object
  return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

def getDischarge(tag,loc,st,et) :
  """Loads discharge data from disk using given tag, station name and time
  range and applies low pass filter. Variable is assumed to be \'flux\'."""
  disch = dirTreeManager.getDataContainer(tag=tag,dataType='timeseries',
                               location=loc,variable='flux',
                               startTime=st,endTime=et)

  disch.data /= 1000.0 # change unit: 1000 m3 -s
  dischB = removeTides( disch ) #, 44714.0*2 )
  return dischB

def getTidalRange(tag,loc,st,et) :
  """Loads elevation data from disk using given tag, station name and time
  range and computes tidal range."""
  T = 44714. # M2 period in seconds
  pad = datetime.timedelta(seconds=3*T)
  tideDC = dirTreeManager.getDataContainer(tag=tag,dataType='timeseries',
                                location=loc,variable='elev',
                                startTime=st-pad,endTime=et+pad)
  tidalRange = timeSeriesFilters.runningRange(tideDC,T=2*T)
  tidalRange = timeSeriesFilters.removeTides(tidalRange) #,T=3*T)
  tidalRange.setMetaData('variable','tidal_range')
  return tidalRange

def getDiff( ref, sig, refTag=None, sigTag=None, absolute=False ) :
  # compute difference
  diff = ref.computeError( sig )
  if absolute :
    diff.data=np.abs(diff.data)
  if sigTag is None :
    sigTag = sig.getMetaData('tag').split('-')[0]
  if refTag is None :
    refTag = ref.getMetaData('tag').split('-')[0]
  tagStr = sigTag+'-'+refTag
  diff.setMetaData('tag',tagStr)
  return diff

def plotScatter( ax, xDC, yDC, color=None, typ='both', cbarArgs={}, **kwargs ) :
  year = xDC.time.getDatetime(0).year
  zorder = kwargs.pop('zorder',5)
  if typ in ['line','both'] :
    _,ranges,_ = xDC.detectGaps()
    for i in range(ranges.shape[0]) :
      ax.plot( np.squeeze(xDC.data[:,:,ranges[i,0]:ranges[i,1]]),
              np.squeeze(yDC.data[:,:,ranges[i,0]:ranges[i,1]]), 'k-',lw=0.7,zorder=zorder+1 )
  x = np.squeeze(xDC.data)
  y = np.squeeze(yDC.data)
  if typ in ['scatter', 'both'] :
    colorData = False
    if color is None :
      color=np.linspace(0,1,len(x))[::-1]
      colorData=True
    elif isinstance( color, dataContainer.dataContainer ) :
      color = np.squeeze( color.data )[::-1]
      colorData=True
    #plt.scatter( x, y, c=color, edgecolors='none',s=40 )
    kw = dict(kwargs)
    loc = yDC.getMetaData('location')
    kw.setdefault('ylabel','Discharge ('+loc+')')
    kw.setdefault('yunit','1000x m3 s-1')
    loc = xDC.getMetaData('location')
    kw.setdefault('xlabel','Tidal Range ('+loc+')')
    kw.setdefault('xunit','m')
    kw.setdefault('xIsTime',False)
    kw.setdefault('s',30)
    dia = trackPlot.trackTimeSeriesPlot(**kw)
    dia.setAxes( ax )
    dia.addSample( x[::-1], y[::-1], color, zorder=zorder )
    if colorData :
      dia.showColorBar( **cbarArgs )

def makePlot(ax,xDC,yDC,colorDC,xDCline,yDCline,clim=None, prefix=None, **kwargs):
  ax.set_axis_bgcolor('Gainsboro')
  var = colorDC.fieldNames[0]
  tag = colorDC.getMetaData('tag')
  loc = colorDC.getMetaData('location')
  if loc != 'plume' :
    #clabel = VARS[var]+' ('+tag+' '+loc+')'
    clabel = VARS.get(var,var)+' ('+loc.replace('saturn','saturn-').upper()+')'
  else :
    clabel = VARS.get(var,var)+' ('+tag+')'
  if prefix : clabel = prefix+' '+clabel
  unit = UNITS.get(var,'')
  st = xDC.time.getDatetime(0)
  et = xDC.time.getDatetime(-1)
  # three way alignment of data
  newTime = colorDC.time.copy()
  if len(xDC.time) > len(newTime) : newTime = xDC.time.copy()
  if len(yDC.time) > len(newTime) : newTime = yDC.time.copy()
  xTimeIx = xDC.time.getAlignedTimeIndices( newTime )
  yTimeIx = yDC.time.getAlignedTimeIndices( newTime )
  cTimeIx = colorDC.time.getAlignedTimeIndices( newTime )
  commonIx = np.intersect1d( np.intersect1d( xTimeIx, yTimeIx ), cTimeIx )
  newTime.array = newTime.array[commonIx]
  xDC_aligned = xDC.interpolateInTime( newTime )
  yDC_aligned = yDC.interpolateInTime( newTime )
  colorDC_aligned = colorDC.interpolateInTime( newTime )
  
  zorder = kwargs.pop('zorder',5)
  # entire data set with black line
  plotScatter(ax, xDCline,yDCline, typ='line',zorder=zorder)
  # add color line
  plotScatter(ax, xDC_aligned,yDC_aligned, color=colorDC_aligned,
              typ='scatter',s=40, clabel=clabel, clim=clim, unit=unit, zorder=zorder,**kwargs)

  #ax.text(0.92, 0.96, tag.upper(), weight='bold',fontsize=fontsize+2,
          #verticalalignment='top', horizontalalignment='right',
          #transform=ax.transAxes)
  ax.set_title( tag.upper(), weight='bold', fontsize=fontsize+2 )
  return ax

def saveRibbonPlot(xDC,yDC,colorDC, imgDir='plots',prefix='ribbon',fpref=''):
  st = xDC.time.getDatetime(0)
  et = xDC.time.getDatetime(-1)
  figExt = ['png']
  var = colorDC.fieldNames[0]
  tag = colorDC.getMetaData('tag')
  loc = colorDC.getMetaData('location')
  msldepth = colorDC.getMetaData('msldepth')
  xstr='-'.join([xDC.getMetaData('tag'),xDC.getMetaData('location')])
  ystr='-'.join([yDC.getMetaData('tag'),yDC.getMetaData('location')])
  dateStrShort = st.strftime('%Y-%b')+'-'+et.strftime('%b')
  locStr = '-'.join([loc,msldepth])
  fname = '_'.join([prefix,xstr,ystr,tag,locStr,var,fpref,dateStrShort])
  saveFigure( imgDir,fname,figExt,verbose=True, dpi=200, bbox_tight=True )
  #plt.close()

def loadAndPrepareData(tideStation,tideTag,dischargeStation,dischargeTag,
                 referenceTag,testTag,variable,station,msldepth,
                 startTime,endTime) :
  # load time series for x/y axis
  def printGaps(dc) :
    print '* Contiguous data intervals', dc.getMetaData('tag'), dc.getMetaData('location'),dc.getMetaData('variable')
    _,ranges,_ = dc.detectGaps()
    for i in range(ranges.shape[0]) :
      print i, dc.time.getDatetime(ranges[i,0]),dc.time.getDatetime(ranges[i,1])

  ydata = getDischarge(dischargeTag,dischargeStation,startTime,endTime)
  printGaps(ydata)
  xdata = getTidalRange(tideTag,tideStation,startTime,endTime)
  printGaps(xdata)
  xdata2,ydata2 = xdata.alignTimes( ydata )

  cmap = cmap_discretize(plt.get_cmap(),64)

  # load data for comparison
  rule='montlyFile'
  dataType = 'sil' if variable[:3]=='sil' else 'timeseries'
  refDC = dirTreeManager.getDataContainer(tag=referenceTag,dataType=dataType,
                              location=station,variable=variable,rule=rule,
                              msldepth=msldepth, startTime=startTime,endTime=endTime)
  testDC = dirTreeManager.getDataContainer(tag=testTag,dataType=dataType,
                              location=station,variable=variable,rule=rule,
                              msldepth=msldepth, startTime=startTime,endTime=endTime)
  diff = getDiff( refDC, testDC, absolute=True)
  
  return xdata,ydata,xdata2,ydata2,diff

#-------------------------------------------------------------------------------
# Main routine
#-------------------------------------------------------------------------------

def doRibbonPlot(tideStation,tideTag,dischargeStation,dischargeTag,
                 referenceTag,testTag,variable,station,msldepth,
                 startTime,endTime,imgDir,clim=None) :

  xdata,ydata,xdata2,ydata2,diff = loadAndPrepareData(
                            tideStation,tideTag,dischargeStation,dischargeTag,
                            referenceTag,testTag,variable,station,msldepth,
                            startTime,endTime)

  for i in range(len(diff.fieldNames)) :
    colorDC = diff.getFields(i)
    
    #fig = plt.figure(figsize=(8,8))
    #ax = fig.add_subplot(111)
    #makePlot(ax,xdata,ydata,colorDC,xdata2,ydata2,None, prefix='Error')
    #saveRibbonPlot(xdata,ydata,colorDC, imgDir,prefix='ribbon',fpref='err')

    colorDC = timeSeriesFilters.removeTides(colorDC) #,T=2*44714.)
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    makePlot(ax,xdata,ydata,colorDC,xdata2,ydata2,clim, prefix='LP Error')
    saveRibbonPlot(xdata2,ydata2,colorDC, imgDir,prefix='ribbon',fpref='err_lpf')

#------------------------------------------------------------------------------
#Command line interface
#------------------------------------------------------------------------------
def parseCommandLine() :

  from optparse import OptionParser
  usage = 'Usage: %prog [options] refTag testTag'
  parser = OptionParser(usage=usage)
  parser.add_option('-t', '--tide-station', action='store', type='string',
                      dest='tideStation', help='station for tidal data')
  parser.add_option('-d', '--discharge-station', action='store', type='string',
                      dest='dischargeStation', help='station for river discharge data')
  parser.add_option('', '--tide-tag', action='store', type='string',
                      dest='tideTag', default='obs',
                      help='tag defining the origin of tidal data (default: %default)')
  parser.add_option('', '--discharge-tag', action='store', type='string',
                      dest='dischargeTag', default='obs',
                      help='tag defining the origin of discharge data (default: %default)')
  parser.add_option('-o', '--output-directory', action='store', type='string',
                      dest='imgDir', help='base directory where generated figures are stored')
  parser.add_option('-v', '--variable', action='store', type='string',
                      dest='var', help='time series to compare: variable (elev,salt,hvel,...)')
  parser.add_option('-a', '--station', action='store', type='string',
                      dest='station', help='time series to compare: station (jetta,saturn01,...)')
  parser.add_option('-m', '--msldepth', action='store', type='string',
                      dest='msldepth', help='time series to compare: mean sea level depth (640,1950,...)')
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
  parser.add_option('-c', '--clim', action='store', type='string',
                      dest='clim', help='Limits for the color bar, e.g. [0.2,12]')

  (opt, args) = parser.parse_args()

  if len(args)<1 :
    parser.print_help()
    parser.error('refTag missing')
  if len(args)<2 :
    parser.print_help()
    parser.error('testTag missing')

  refTag = args[0]
  testTag = args[1]

  if not opt.tideStation :
    parser.print_help()
    parser.error('tide station undefined')
  if not opt.dischargeStation :
    parser.print_help()
    parser.error('discharge station undefined')
  if not opt.imgDir :
    parser.print_help()
    parser.error('output directory undefined')
  if not opt.var :
    parser.print_help()
    parser.error('variable undefined')
  if not opt.station :
    parser.print_help()
    parser.error('station undefined')
  #if not opt.msldepth :
    #parser.print_help()
    #parser.error('msldepth undefined')
  if not opt.startStr :
    parser.print_help()
    parser.error('start time undefined')
  if not opt.endStr :
    parser.print_help()
    parser.error('end time undefined')

  if opt.startStr :
    startTime = parseTimeStr( opt.startStr )
  if opt.endStr :
    endTime = parseTimeStr( opt.endStr )
  if opt.clim :
    clim = [ float(v) for v in opt.clim.strip('[').strip(']').split(',') ]

  print 'Parsed options:'
  print ' - tidal data:', opt.tideTag, opt.tideStation
  print ' - discharge data:', opt.dischargeTag, opt.dischargeStation
  print ' - time series to compare:', opt.var, opt.station, opt.msldepth
  print ' - reference data set:', refTag
  print ' - test data set:', testTag
  print ' - time range:', startTime, endTime
  print ' - output dir:',opt.imgDir
  
  if opt.clim :
    print ' - color limits:',clim
  else :
    clim = None

  if opt.imgDir :
    createDirectory(opt.imgDir)

  doRibbonPlot(opt.tideStation,opt.tideTag,opt.dischargeStation,
               opt.dischargeTag, refTag,testTag,opt.var,
               opt.station,opt.msldepth, startTime,endTime,opt.imgDir,clim)

if __name__=='__main__' :
  parseCommandLine()
