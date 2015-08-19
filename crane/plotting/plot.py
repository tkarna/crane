"""
High-level class to create timeseries and Taylor Diagrams.
"""
#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import os
import numpy as np
import matplotlib.pyplot as plt

from crane.data import timeArray
from crane.plotting import plotBase
from crane.plotting import timeSeriesPlot
from crane.plotting import taylorDiagram
from crane.plotting import stationExtremaPlot
from crane.plotting import errorHistogram
from crane.plotting import spectralPlot
from crane.plotting import tidalConstPlot
from crane.plotting import profilePlot
from crane.plotting import trackPlot
from crane.files import stationFile

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
# These strings will appear on plot axes etc
VARS = {'salt': 'Salinity',
        'temp': 'Temperature',
        'cond': 'Electrical conductivity',
        'elev': 'Elevation',
        'pres': 'Pressure',
        'fluores': 'Fluorescence',
        'kine': 'TKE',
        'vdff': 'Viscosity',
        'tdff': 'Diffusivity',
        'mixl': 'Mixing length',
        'hvel': 'Velocity',
        'vert': 'Vert velocity',
        'u': 'u-velocity',
        'v': 'v-velocity',
        'w': 'w-velocity',
        'alongvel': 'Along Vel.',
        'discharge': 'Discharge',
        'tidal_range': 'Tidal range',
        'cui': 'Upwelling index',
        'strat': 'Stratification',
        'sil': 'Salt intrusion length',
        'srad': 'Solar radiation',
        'sed': 'Sediment',
        'bed_stress': 'bed stress',
        'bed_depth': 'bed depth',
        'sed_1': 'S-class 1',
        'sed_2': 'S-class 2',
        'sed_3': 'S-class 3',
        'sed_4': 'S-class 4',
        'di_sed': 'DI sed',
        'NO3': 'NO3',
        'NH4': 'NH4',
        'PHYM': 'PHYM',
        'PHYF': 'PHYF',
        'SZOO': 'SZOO',
        'BZOO': 'BZOO',
        'DETN': 'DETN',
        'DETC': 'DETC',
        'BACT': 'BACT',
        'DON': 'DON',
        'DOC': 'DOC',
        'CHLM': 'CHLM',
        'CHLF': 'CHLF',
        'bnth_1':'NO3 benthos',
        'bnth_2':'NH4 benthos',
        'bnth_3':'DETN benthos',
        'bnth_4':'OXY benthos',
        'Diag':'Diagnostic array',
        'nem':'NEM',
        'nem_DAVG':'Averaged NEM',
        'nem_DI':'Integrated NEM',
        'nemb':'benthic NEM',
        'nemi':'vert. int. NEM',
        'totalN':'total Nitrogen',
        'NO3': 'NO3',
        'NH4': 'NH4',
        'phy': 'Phy',
        'zoo': 'Zoo',
        'det': 'Det',
        'oxy': 'Oxy',
        'plumevol': 'Plume volume',
        'plumethick': 'Plume thickness',
        'flntu':'FLNTU',
        'ecoblue':'ECO backscatter 470nm',
        'ecored':'ECO backscatter 700nm',
        'dens':'Density',
        'turbidity': 'Sediment',
        'bath': 'Bathymetry',
        'depth': 'Depth',
        'bed_depth': 'Depth',
        'bed_stress': 'bottom stress',
        'sil_1':'SIL 1',
        'sil_2':'SIL 2',
        'sil_3':'SIL 3',
        'sil_4':'SIL 4',
        'sil_5':'SIL 5',
        'plume_area_12':'Plume area 12',
        'plume_volume_12':'Plume volume 12',
        'plume_center_12':'Plume center 12',
        'plume_center_x_12':'Plume center x 12',
        'plume_center_y_12':'Plume center y 12',
        'plume_thickness_12':'Plume thickness 12',
        'plume_area_28':'Plume area 28',
        'plume_volume_28':'Plume volume 28',
        'plume_center_28':'Plume center 28',
        'plume_center_x_28':'Plume center x 28',
        'plume_center_y_28':'Plume center y 28',
        'plume_thickness_28':'Plume thickness 28',
        'plume_area_30':'Plume area 30',
        'plume_volume_30':'Plume volume 30',
        'plume_center_30':'Plume center 30',
        'plume_center_x_30':'Plume center x 30',
        'plume_center_y_30':'Plume center y 30',
        'plume_thickness_30':'Plume thickness 30',
        'volume': 'Volume',
        'flux': 'Flux',
        'bed_flux':'Sediment bed flux',
        'iRiv':'Riverine water',
        'iOce':'Oceanic water',
        'iPlu':'Plume water',
        'iRen':'Renewed water',
        'aRiv':'Riv. age conc.',
        'aOce':'Oce. age conc.',
        'aPlu':'Plu. age conc.',
        'aRen':'Ren. age conc.',
        'ARiv':'Riv. age',
        'AOce':'Oce. age',
        'APlu':'Plu. age',
        'ARen':'Ren. age',
        'maxiRiv':'Riverine water',
        'maxiOce':'Oceanic water',
        'maxiPlu':'Plume water',
        'maxiRen':'Renewed water',
        'maxaRiv':'Riv. age conc.',
        'maxaOce':'Oce. age conc.',
        'maxaPlu':'Plu. age conc.',
        'maxared':'Ren. age conc.',
        'maxARiv':'Riv. age',
        'maxAOce':'Oce. age',
        'maxAPlu':'Plu. age',
        'maxARen':'Ren. age',
        }

# These strings will appear as units in plots
UNITS = {'salt': 'psu',
        'temp': 'degC',
        'cond': '',
        'elev': 'm',
        'pres': 'Pa',
        'fluores': '',
        'kine': 'm^2/s^2',
        'vdff': 'm^2/s',
        'tdff': 'm^2/s',
        'mixl': 'm',
        'hvel': 'm/s',
        'vert': 'm/s',
        'u': 'm/s',
        'v': 'm/s',
        'w': 'm/s',
        'alongvel': 'm/s',
        'discharge': 'm^3/s',
        'tidal_range': 'm',
        'cui': '',
        'strat':'psu',
        'sil':'km',
        'srad':'W m-2',
        'sed':'mg/l',
        'bed_stress':'N/m^2',
        'bed_depth':'m',
        'sed_1':'kg/m^3',
        'sed_2':'kg/m^3',
        'sed_3':'kg/m^3',
        'sed_4':'kg/m^3',
        'sed_5':'kg/m^3',
        'di_sed':'kg/m^3',
        'NO3':'mmol N m-3',
        'NH4':'mmol N m-3',
        'PHYM':'mmol N m-3',
        'PHYF':'mmol N m-3',
        'SZOO':'mmol N m-3',
        'BZOO':'mmol N m-3',
        'DETN':'mmol N m-3',
        'DETC':'mmol C m-3',
        'BACT':'mmol N m-3',
        'DON':'mmol N m-3',
        'DOC':'mmol C m-3',
        'CHLM':'mg Chl m-3',
        'CHLF':'mg Chl m-3',
        'bnth_1':'mmol N m-3',
        'bnth_2':'mmol N m-3',
        'bnth_3':'mmol N m-3',
        'bnth_4':'mmol O m-3',
        'Diag':' ',
        'nem':'mmol O m-3 day-1 ',
        'nem_DAVG':'mmol O m-3 day-1 ',
        'nem_DI':'mmol O m-2 day-1 ',
        'nemb':'mmol O m-2 day-1 ',
        'nemi':'mmol O m-2 day-1 ',
        'totalN':'mmol N m-3',
        'trcr_5':'kg/m^3',
        'phy': 'mmol N m-3',
        'zoo': 'mmol N m-3',
        'det': 'mmol N m-3',
        'oxy': 'mmol m-3',
        'plumevol': '10^9 m^3',
        'plumethick': 'm',
        'flntu':'ntu',
        'ecored':'m-1 sr-1',
        'ecoblue':'m-1 sr-1',
        'dens':'kg m-3',
        'turbidity': 'kg/m^3',
        'bath':'m',
        'depth':'m',
        'sil_1':'km',
        'sil_2':'km',
        'sil_3':'km',
        'sil_4':'km',
        'sil_5':'km',
        'plume_area_12':'m2',
        'plume_volume_12':'m3',
        'plume_center_12':'m',
        'plume_center_x_12':'m',
        'plume_center_y_12':'m',
        'plume_thickness_12':'m',
        'plume_area_28':'m2',
        'plume_volume_28':'m3',
        'plume_center_28':'m',
        'plume_center_x_28':'m',
        'plume_center_y_28':'m',
        'plume_thickness_28':'m',
        'plume_area_30':'m2',
        'plume_volume_30':'m3',
        'plume_center_30':'m',
        'plume_center_x_30':'m',
        'plume_center_y_30':'m',
        'plume_thickness_30':'m',
        'volume': 'm3',
        'flux': 'm3 s-1',
        'bed_flux':'m3',
        'iRiv':'frac',
        'iOce':'frac',
        'iPlu':'frac',
        'iRen':'frac',
        'aRiv':'frac s',
        'aOce':'frac s',
        'aPlu':'frac s',
        'aRen':'frac s',
        'ARiv':'h',
        'AOce':'h',
        'APlu':'h',
        'ARen':'h',
        'maxiRiv':'frac',
        'maxiOce':'frac',
        'maxiPlu':'frac',
        'maxiRen':'frac',
        'maxaRiv':'frac s',
        'maxaOce':'frac s',
        'maxaPlu':'frac s',
        'maxaRen':'frac s',
        'maxARiv':'h',
        'maxAOce':'h',
        'maxAPlu':'h',
        'maxARen':'h',
        }

COLORS = ['b', 'g', 'k', 'm', 'c', 'y',
          'Darkorange','Orchid','SpringGreen','DarkTurquoise','DarkKhaki',
          'Coral','CadetBlue', 'GoldenRod'] # HTML names
MARKERS = ['o','^','*','d','s','v','h','>','<','p','H','+','x']

#-------------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------
def getAnnotationStrings( sColl, var, varPrefix ) :
  # var for figure title
  varTitle = VARS.get(var,var).capitalize()
  # var for image filename
  varStr = var
  if varPrefix :
    varTitle = varPrefix.capitalize() +' '+ varTitle
    varStr = varPrefix[:4]+'_'+var
  # start/end date for filename
  sT = str(sColl.startTime.date())
  eT = str(sColl.endTime.date())
  return varTitle,varStr,sT,eT

def makeColorsForModels( sColl ):
  """Given stationCollection, returns a color dictionary of
      unique colors where the key is modelTag.
  """
  modTags = sorted(sColl.getModelTags(),reverse=True)
  mColors = dict(zip(modTags,COLORS))
  mColors[sColl.getObsTag()] = 'r' # add hard-coded color for obs
  return mColors

def getComparableStations( coll, dataType='timeseries', requireObs=True ) :
  if not requireObs :
    return coll.getAttributes( 'location' )
  else :
    # All stations that are in obs and any model
    obsStations = set()
    modStations = set()
    for key in coll.getKeys( dataType=dataType, tag=coll.getObsTag() ) :
      obsStations.add( key['location'] )
    for key in coll.getKeys( dataType=dataType, tag=coll.getModelTags() ) :
      modStations.add( key['location'] )
    stations = list( obsStations.intersection( modStations ) )
    return stations

def makeMarkersForTaylor(coll):
  """Given models, stations and depth, returns a marker and color dictionary
      of unique markers where the key is the tuple (modelTag, station, msldepth).
  """
  markerDic ={} 
  colorDic ={} 
  colors = iter(COLORS)
  marks = iter(MARKERS)

  depths = dict()
  for sta in getComparableStations(coll) :
    depths[sta] = set()
    for key in coll.getKeys( location=sta ) :
      msldepth = key['msldepth']
      depths[sta].add( msldepth )
    depths[sta] = list(depths[sta])

  # color by model, marker by (station,msldepth)
  for modTag in sorted(coll.getModelTags(),reverse=True) :
    try:
      c = colors.next()
    except StopIteration as e:
      raise StopIteration('Not enough defined colors')
    marks = iter(MARKERS)
    for sta in getComparableStations(coll) :
      for msldepth in depths[sta] :
        try:
          m = marks.next()
        except StopIteration:
          c = colors.next()
          marks = iter(MARKERS)
          m = marks.next()
        markerDic[(modTag,sta,msldepth)] = m
        colorDic[(modTag,sta,msldepth)] = c

  return markerDic, colorDic


class Plots(object):
  """High level class to create timeseries and Taylor Diagram plots.
  
  Attributes:
    path -- String of path to image directory.
    coll -- StationCollection object
  """
  def __init__(self, path, sCollection):
    """Creates and saves a set of plots.

    Args:
      path -- (string) path to image directory.
      sCollection -- (StationCollection) collection of all data
    """
    if not os.path.isdir(path):
      raise Exception('Supplied image directory does not exist: '+path)
    self.path = path
    self.coll = sCollection

  def getComparableStations( self, dataType='timeseries', requireObs=True ) :
    return getComparableStations(self.coll)

  def getAnnotationStrings( self, var, varPrefix ) :
    return getAnnotationStrings( self.coll, var, varPrefix )

  def saveFigure( self, filename, extensions, **kwargs ) :
    """Saves currently open figure in the given format.
    If extensions is a list, figure is saved in multiple formats."""
    kw = dict(kwargs)
    kw.setdefault('verbose', True)
    kw.setdefault('bbox_tight', True)
    plotBase.saveFigure(self.path, filename, extensions, **kw)

  def makeSaltIntrusion(self, filetype='png', fPrefix='sil', varPrefix=None, ylim=None, err_ylim=None, **kwargs) :
    """Creates timeseries of salt intrusion length and saves to image directory.

    Files saved to <self.path>/sil_<location>.png OR
    Files saved to <self.path>/<fPrefix>_<station>.png

    Args:
      kwargs  -- Dictionary of plot options.
      fPrefix -- String of prefix for image file
      xlim    -- limits for x axis, [s, e] where s and e are datetime objects
      ylim    -- Dictionary of y axis limits for each variable
      err_ylim -- Dictionary of error axis limits for each variable
    Outputs:
      Saves images to self.path directory
    """
    modelColors = self.makeColorsForModels()
    # create one plot for each SIL transect ('location')
    locList = self.coll.getAttributes('location',dataType='sil')
    for location in locList :
      keys = self.coll.getKeys(dataType='sil',location=location)
      var = 'sil'
      # ----- Prep plot annotation -----
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      kw = dict(kwargs)
      if not 'title' in kw:
        kw['title'] = ('%s %s\n %s - %s (PST)' %
                            (location, varTitle, sT, eT))
      fig = plt.figure(figsize=(7.5,3))
      dia = timeSeriesPlot.timeSeriesPlotDC(varTitle, unit=UNITS[var], ylim=ylim, **kw)

      for modKey in keys :
        m = self.coll.getSample( **modKey )
        isMaxSIL = m.getMetaData('variable')=='max_sil'
        linewidth = 1.0 if isMaxSIL else 0.4
        lineTag = 'max sil' if isMaxSIL else 'sil'
        dia.addSample(m, color=modelColors[modKey['tag']], label=modKey['tag']+' '+lineTag, linewidth=linewidth)
      dia.makePlot()

      # ----- Save file -----
      file = '_'.join([fPrefix,location,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeTimeSeries(self, filetype='png', fPrefix='ts', varPrefix=None, ylim=None, err_ylim=None, **kwargs) :
    """Creates timeseries of offerings and saves to image directory.

    Files saved to <self.path>/ts_<station>_<var>_<depth>.png OR
    Files saved to <self.path>/<fPrefix>_<station>_<var>_<depth>.png 

    Args:
      kwargs  -- Dictionary of plot options.
      fPrefix -- String of prefix for image file
      xlim    -- limits for x axis, [s, e] where s and e are datetime objects
      ylim    -- Dictionary of y axis limits for each variable
      err_ylim -- Dictionary of error axis limits for each variable
    Outputs:
      Saves images to self.path directory
    """
    modelColors = self.makeColorsForModels()

    # create one plot for each (observation, models) pair
    comKeys = self.coll.getComparableKeys(dataType='timeseries',requireObs=False)
    for entry, obsKey, modKeys in comKeys :
      station,var,msldepth = entry
      if obsKey :
        o = self.coll.getSample( **obsKey )
        if len(o.time) < 3 :
          print 'observation time series too short:',obsKey
          continue
      else :
        o = None
      
      # ----- Prep plot annotation -----
      depth = max(float(msldepth),0)/100
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      kw = dict(kwargs)
      if not 'title' in kw:
        kw['title'] = ('%s %s Depth: %.1f m\n %s - %s (PST)' %
                            (station.upper(), varTitle, depth, sT, eT))
      yl = ylim.get(var) if ylim else None
      erryl = err_ylim.get(var) if err_ylim else None

      # ----- Make image -----
      if obsKey :
        dia=timeSeriesPlot.timeSeriesComboPlotDC(varTitle, o, 'Obs', unit=UNITS[var], ylim=yl, err_ylim=erryl, **kw)
      else :
        fig = plt.figure(figsize=(7.5,3))
        dia=timeSeriesPlot.timeSeriesPlotDC(varTitle, unit=UNITS[var], ylim=yl, **kw)
      for modKey in modKeys :
        m = self.coll.getSample( **modKey )
        dia.addSample(m, color=modelColors[modKey['tag']], label=modKey['tag'])
      dia.makePlot()

      # ----- Save file -----
      file = '_'.join([fPrefix,station,varStr,msldepth,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeTimeSeriesStack(self, stationX, filetype='png', fPrefix='tsStack', varPrefix=None, plotError=False, **kwargs) :
    """Creates a time series plot with a subplot for each station.
    """
    modelColors = self.makeColorsForModels()

    for var in self.coll.getAttributes('variable', dataType='timeseries') :
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      titleStr = ' '.join([ varTitle, sT+' - '+eT])
      if 'title' not in kwargs :
        kwargs['title'] = titleStr
      dia = timeSeriesPlot.timeSeriesStackPlotDC(**kwargs)
      stations = list( set( self.coll.getAttributes( 'location', variable=var ) ).intersection( set( stationX.keys() ) ) )
      orderedSta = sorted( stations, key = lambda k: stationX[k] )
      for sta in orderedSta :
        comKeys = self.coll.getComparableKeys( dataType='timeseries',variable=var, location=sta, requireObs=False )
        for entry, obsKey, modKeys in comKeys :
          station,var,msldepth = entry

          #dia.addSubplot( sta, VARS[var], UNITS[var] ) # ylabel = var+unit
          #dia.addSubplot( sta, sta.capitalize(), '' ) # ylabel = sta
          dia.addSubplot( sta, sta.upper()+' '+msldepth, '' ) # ylabel = sta msldepth
          if obsKey :
            o = self.coll.getSample( **obsKey )
            if len(o.time) < 3 :
              print 'observation time series too short:',obsKey
              continue
            dia.addSample( sta, o, color='r', label='obs')
          for modKey in modKeys :
            m = self.coll.getSample( **modKey )
            dia.addSample( sta, m, color=modelColors[modKey['tag']], label=modKey['tag'])
            if obsKey and plotError :
              try :
                err = o.computeError( m )
                errColor = (0.5,0.5,0.5)
                #dia.addSample( sta, err, linestyle='dashed', color=modelColors[modTag], label='err '+modTag)
                dia.addSample( sta, err, color=errColor, label='err '+modKey['tag'])
              except Exception as e:
                print 'Error could not be computed, skipping:', modKey['tag'],sta,var,msldepth
                print e
      if dia.isEmpty() :
        print 'TimeSeriesStack: diagram is empty, skipping', var
        continue
      dia.makePlot()

      # ----- Save file -----
      file = '_'.join([fPrefix,varStr,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeStationExtremaPlot(self, filetype='png', fPrefix='minmax', varPrefix=None, stationCoords=None, requireObs=True, **kwargs) :
    """Creates and saves stationExtrema plot, showing the min/mean/max for each station, versus distance from river mouth. One figure is produced for each variable.
    
    Args:
      kwargs  -- Dictionary of plot options.
      fPrefix -- String of prefix for image file
    """
    
    modelColors = self.makeColorsForModels()

    # generate station coordinates TODO remove use stationX
    # for now, just take x coord
    if not stationCoords :
      stationCoords = dict()
      sf = StationFile()
      for station in self.getComparableStations(requireObs=False) :
        x, y = sf.getLocation( station )
        stationCoords[station] = x

    for var in self.coll.getAttributes('variable') :
      dia = stationExtremaPlot.stationExtremaPlotDC(VARS[var], stationCoords, unit=UNITS[var], **kwargs)
      comKeys = self.coll.getComparableKeys(dataType='timeseries', requireObs=False)
      for entry, obsKey, modKeys in comKeys :
        station,var,msldepth = entry
        keys = []
        if obsKey :
          keys.append( obsKey )
        keys.extend( modKeys )
        for key in keys :
          s = self.coll.getSample( **key )
          dia.plotSample( key['tag'], station, s, color=modelColors[key['tag']])
      if dia.isEmpty() :
        print 'StationExtremaPlot: diagram is empty, skipping', var
        continue
    
      dia.showLegend(prop=dict(size='small'), loc='best')
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      title = ('%s Dates: %s - %s (PST)' % (varTitle, sT, eT))
      dia.addTitle(title, size=12)

      # ----- Save file -----
      file = '_'.join([fPrefix,varStr,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeErrorHistograms(self, filetype='png', fPrefix='eh', varPrefix=None, **kwargs) :
    """Creates an error histogram plot for all stations/variables/depths"""
    modelColors = self.makeColorsForModels()

    # create one plot for each (observation, models) pair
    comKeys = self.coll.getComparableKeys( dataType='timeseries', requireObs=True )
    for entry, obsKey, modKeys in comKeys :
      station,var,msldepth = entry

      o = self.coll.getSample( **obsKey )
      if len(o.time) < 3 :
        print 'observation time series too short:',obsKey
        continue

      # ----- Make plot -----
      dia = errorHistogram.errorHistogramDC(o, unit=UNITS[var], label='Obs', **kwargs)
      for modKey in modKeys :
        m = self.coll.getSample( **modKey )
        dia.addSample(m, label=modKey['tag'], color=modelColors[modKey['tag']])
      dia.makePlot()
      dia.showLegend()

      # ----- Prep plot annotation -----
      depth = max(float(msldepth),0)/100
      dia.showLegend(prop=dict(size='small'))
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      title = ('%s %s Depth: %.1f m Dates: %s - %s (PST)' %
                  (station.upper(), varTitle, depth, sT, eT))
      dia.addTitle(title, size=12)

      # ----- Save file -----
      file = '_'.join([fPrefix,station,varStr,msldepth,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeErrorHistogramStack(self, stationX, filetype='png', fPrefix='ehStack', varPrefix=None, **kwargs) :
    modelColors = self.makeColorsForModels()

    for var in self.coll.getAttributes('variable',dataType='timeseries') :
      dia = errorHistogram.stackHistogramDC(unit=UNITS[var],**kwargs)
      stations = list( set( self.coll.getAttributes( 'location', variable=var ) ).intersection( set( stationX.keys() ) ) )
#      stations = self.coll.getAttributes( 'location',variable=var )
      orderedSta = sorted( stations, key = lambda k: stationX[k] )
      for sta in orderedSta :
        comKeys = self.coll.getComparableKeys( dataType='timeseries', variable=var, location=sta, requireObs=True )
        for entry, obsKey, modKeys in comKeys :
          station,var,msldepth = entry
          o = self.coll.getSample( **obsKey )
          if len(o.time) < 3 :
            print 'observation time series too short:',obsKey
            continue

          #dia.addSubplot( sta, VARS[var], UNITS[var] ) # ylabel = var+unit
          #dia.addSubplot( sta, sta.capitalize(), '' ) # ylabel = sta
          #dia.addSubplot( sta, sta.upper()+' '+msldepth, '' ) # ylabel = sta msldepth
          dia.addPlot(sta+msldepth,o, title=sta.upper()+' '+msldepth, ylabel='freq')
          for modKey in modKeys :
            m = self.coll.getSample( **modKey )
            dia.addSample(sta+msldepth,m,modKey['tag'],color=modelColors[modKey['tag']])
      if dia.isEmpty() :
        print 'ErrorHistogramStack: diagram is empty, skipping', var
        continue
      dia.makePlot()
      dia.showLegend()
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      titleStr = ' '.join([ varTitle, sT+' - '+eT])
      dia.addTitle( titleStr, size=14)

      # ----- Save file -----
      file = '_'.join([fPrefix,varStr,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeErrorSpectralPlots(self, filetype='png', fPrefix='esp', varPrefix=None, **kwargs) :
    """Creates an error histogram plot for all stations/variables/depths"""
    modelColors = self.makeColorsForModels()
    constsToPlot = ['M2','S2','N2','L2','K2','K1','O1','P1','M4','M6','M8','Mf','Mm']
    # create one plot for each (observation, models) pair
    comKeys = self.coll.getComparableKeys( dataType='timeseries', requireObs=True )
    for entry, obsKey, modKeys in comKeys :
      station,var,msldepth = entry
      o = self.coll.getSample( **obsKey )
      if len(o.time) < 3 :
        print 'observation time series too short:',obsKey
        continue

      # ----- Make plot -----
      dia = spectralPlot.spectralPlotDC(xunit='hours')
      for modKey in modKeys :
        m = self.coll.getSample( **modKey )
        try :
          err = o.computeError( m )
          dia.addSample(err,label=modKey['tag'], color=modelColors[modKey['tag']])
        except Exception as e:
          print 'Error could not be computed, skipping:', modKey['tag'],station,var,msldepth
          print e
      #dia.addConstituentIndicators( constsToPlot )
      dia.showLegend()

      # ----- Prep plot annotation -----
      depth = max(float(msldepth),0)/100
      dia.showLegend(prop=dict(size='small'))
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      title = ('Error spectrum: %s %s Depth: %.1f m Dates: %s - %s (PST)' %
                  (station.upper(), varTitle, depth, sT, eT))
      dia.addTitle(title, size=12)

      # ----- Save file -----
      file = '_'.join([fPrefix,station,varStr,msldepth,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeHarmonicAnalysisErrorPlots(self, stationX, filetype='png', fPrefix='haStack', varPrefix=None, **kwargs) :
    """Creates an harmonic analysis plot for all stations/variables/depths"""
    modelColors = self.makeColorsForModels()
    #constsToPlot = ['M2','S2','N2','L2','K2','K1','O1','P1','M4','M6','M8','Mf','Mm']
    constsToPlot = ['M2','S2','N2','L2','K2','K1','O1','P1','J1','NO1','OO1','M4','M6','M8','Mf','Mm','MSf']
    for var in self.coll.getAttributes('variable',dataType='timeseries') :
      dia = stackAmplitudePhasePlot(unit=UNITS[var],**kwargs)
#      stations = self.coll.getAttributes( 'location',variable=var)
      stations = list( set( self.coll.getAttributes( 'location', variable=var ) ).intersection( set( stationX.keys() ) ) )
      orderedSta = sorted( stations, key = lambda k: stationX[k] )
      for sta in orderedSta :
        comKeys = self.coll.getComparableKeys( dataType='timeseries', variable=var, location=sta, requireObs=True )
        for entry, obsKey, modKeys in comKeys :
          station,var,msldepth = entry

          for modKey in modKeys :
            tc = self.coll.getHarmonicAnalysis( modKey, obsKey )
            if tc :
              if not dia.hasPlot(sta+msldepth) :
                dia.addPlot(sta+msldepth,ylabel=sta.upper())
              dia.addSample(sta+msldepth,tc,label=modKey['tag'],color=modelColors[modKey['tag']])
      if dia.isEmpty() :
        print 'HarmonicAnalysisPlots: diagram is empty, skipping', var
        continue
      dia.makePlot(include=constsToPlot)
      #dia.makePlot(ampThreshold=0.01)
      dia.showLegend()
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      titleStr = ' '.join([ varTitle, sT+' - '+eT])
      dia.addTitle( titleStr, size=14)

      # ----- Save file -----
      file = '_'.join([fPrefix,varStr,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeHarmonicAnalysisPlots(self, stationX, filetype='png', fPrefix='tcStack', varPrefix=None, **kwargs) :
    """Creates an harmonic analysis plot for all stations/variables/depths"""
    modelColors = self.makeColorsForModels()
    constsToPlot = ['M2','M4','M6','S2','K1','O1','MSf']
    for var in self.coll.getAttributes('variable',dataType='timeseries') :
      dia = stackAmplitudePhasePlot(unit=UNITS[var],**kwargs)
      #stations = self.coll.getAttributes( 'location',variable=var)
      stations = list( set( self.coll.getAttributes( 'location', variable=var ) ).intersection( set( stationX.keys() ) ) )
      orderedSta = sorted( stations, key = lambda k: stationX[k] )
      orderedSta = sorted( stations, key = lambda k: stationX[k] )
      for sta in orderedSta :
        comKeys = self.coll.getComparableKeys( dataType='timeseries', variable=var, location=sta, requireObs=True )
        for entry, obsKey, modKeys in comKeys :
          station,var,msldepth = entry
          tc = self.coll.getHarmonicAnalysis( obsKey )
          if not tc :
            continue
          if not dia.hasPlot(sta+msldepth) :
            dia.addPlot(sta+msldepth,ylabel=sta.upper())
          dia.addSample(sta+msldepth,tc,label=self.coll.getObsTag(),color='r')
          for modKey in modKeys :
            tc = self.coll.getHarmonicAnalysis( modKey )
            if tc :
              dia.addSample(sta+msldepth,tc,label=modKey['tag'],color=modelColors[modKey['tag']])
      if dia.isEmpty() :
        print 'HarmonicAnalysisPlots: diagram is empty, skipping', var
        continue
      dia.makePlot(include=constsToPlot)
      #dia.makePlot(ampThreshold=0.01)
      dia.showLegend()
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      titleStr = ' '.join([ varTitle, sT+' - '+eT])
      dia.addTitle( titleStr, size=14)

      # ----- Save file -----
      file = '_'.join([fPrefix,varStr,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeErrorSpectralStack(self, stationX, filetype='png', fPrefix='espStack', varPrefix=None, **kwargs) :
    modelColors = self.makeColorsForModels()

    for var in self.coll.getAttributes('variable',dataType='timeseries') :
      dia = spectralPlot.stackSpectralPlotDC(xunit='hours',**kwargs)
      #stations = self.coll.getAttributes( 'location',variable=var)
      stations = list( set( self.coll.getAttributes( 'location', variable=var ) ).intersection( set( stationX.keys() ) ) )
      orderedSta = sorted( stations, key = lambda k: stationX[k] )
      orderedSta = sorted( stations, key = lambda k: stationX[k] )
      for sta in orderedSta :
        comKeys = self.coll.getComparableKeys( dataType='timeseries', variable=var, location=sta, requireObs=True )
        for entry, obsKey, modKeys in comKeys :
          station,var,msldepth = entry
          o = self.coll.getSample( **obsKey )
          if len(o.time) < 3 :
            print 'observation time series too short:',obsKey
            continue

          dia.addPlot(sta+msldepth, title=sta.upper()+' '+msldepth)
          for modKey in modKeys :
            m = self.coll.getSample( **modKey )
            try :
              err = o.computeError( m )
              dia.addSample(sta+msldepth,err,label=modKey['tag'],color=modelColors[modKey['tag']])
            except Exception as e:
              print 'Error could not be computed, skipping:', modKey['tag'],sta,var,msldepth
              print e
      if dia.isEmpty() :
        print 'ErrorSpectralStack: diagram is empty, skipping', var
        continue
      dia.makePlot()
      #dia.addConstituentIndicators( ['M2','S2','K1','O1','M4','M6','M8','Mf'])
      dia.showLegend()
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      titleStr = ' '.join([ varTitle, sT+' - '+eT])
      dia.addTitle( titleStr, size=14)

      # ----- Save file -----
      file = '_'.join([fPrefix,varStr,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeTaylorDiagrams(self, filetype='png', fPrefix='td', varPrefix=None, **kwargs):
    """Creates and saves a Taylor Diagram comparing model and obs data, one for
    each variable and station.
    Note: This is non-normalized Taylor diagram.
   
     Args:
      kwargs  -- Dictionary of plot options.
      fPrefix -- String of prefix for image file
    Outputs:
      Saves images to self.path directory 
    """
    modelColors = self.makeColorsForModels()

    # create one plot for each (observation, models) pair
    comKeys = self.coll.getComparableKeys( dataType='timeseries', requireObs=True )
    for entry, obsKey, modKeys in comKeys :
      station,var,msldepth = entry
      o = self.coll.getSample( **obsKey )
      if len(o.time) < 3 :
        print 'observation time series too short:',obsKey
        continue

      # ----- Make plot -----
      dia = taylorDiagram.statisticsDiagramDC(o, 'Obs', unit=UNITS[var], **kwargs)
      for modKey in modKeys :
        m = self.coll.getSample( **modKey )
        dia.plotSample(m, label=modKey['tag'], color=modelColors[modKey['tag']])

      # ----- Prep plot annotation -----
      depth = max(float(msldepth),0)/100
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      title = ('%s %s Depth: %.1f m Dates: %s - %s (PST)' %
                  (station.upper(), varTitle, depth, sT, eT))
      dia.addTitle(title, size=12)
      dia.showLegend(prop=dict(size='small'), loc='lower left')

      # ----- Save file -----
      file = '_'.join([fPrefix,station,varStr,msldepth,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeTaylorDiagramsVar(self, filetype='png', fPrefix='tdVar', varPrefix='', **kwargs):
    """Creates and saves a normalized Taylor Diagram comparing model and obs data
      for a single variable at all stations.
     Args:
      kwargs -- Dictionary of plot options.
      fPrefix -- String of prefix for image file
    Outputs:
      Saves images to self.path directory 
   Args:
    """
    markers, colors = self.makeMarkersForTaylor()
    
    for var in self.coll.getAttributes('variable',dataType='timeseries') :
      dia = taylorDiagram.normalizedStatisticsDiagramDC(figsize=(9,6))
      for sta in self.coll.getAttributes( 'location', variable=var ) :
        comKeys = self.coll.getComparableKeys( dataType='timeseries',variable=var, location=sta, requireObs=True )
        for entry, obsKey, modKeys in comKeys :
          msldepth = obsKey['msldepth']
          o = self.coll.getSample( **obsKey )
          if len(o.time) < 3 :
            print 'observation time series too short:',obsKey
            continue
          # ----- Make plot -----
          dia.addDataSet((sta, msldepth), o, 'Obs', **kwargs)
          for modKey in modKeys :
            modTag = modKey['tag']
            mLabel = '%s %s %s' % (modTag, sta.upper(), msldepth)
            col = colors[(modTag,sta,msldepth)]
            mar = markers[(modTag,sta,msldepth)]
            m = self.coll.getSample( **modKey )
            dia.plotSample((sta, msldepth), m, label=mLabel, color=col, marker=mar)

      if dia.isEmpty() :
        print 'TaylorDiagramsVar: diagram is empty, skipping', var
        continue
      # ----- Prep plot annotation -----
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      title = ('%s  Dates: %s - %s (PST)' % (varTitle, sT, eT))
      dia.addTitle(title, size=12)
      dia.showLegend()

      # ----- Save file -----
      file = '_'.join([fPrefix,varStr,sT,eT])
      self.saveFigure( file, filetype,  )
      plt.close()

  def makeVertProfileTimeSeriesPlot(self, filetype='png', fPrefix='vprof', varPrefix='', clim=None, **kwargs) :
    """Plots all extracted vertical profiles versus time, with observation time series (if any)"""
    kwargs.setdefault( 'plotheight',3.0 )
    # find all profiles (any model, no observations)
    allProf = self.coll.getKeys(dataType='profile')
    # get station var pairs
    # (sta,var) -> list of model keys
    staVarProf = dict()
    for modKey in allProf :
      key = (modKey['location'],modKey['variable'])
      staVarProf.setdefault( key, [] ).append( modKey )
    # sort models in same order for consistency
    for key in staVarProf :
      staVarProf[key] = sorted( staVarProf[key] )
    # for each profile
    for station,var in staVarProf :
      modKeys = staVarProf[(station,var)]
      # find all obs that match station and var
      obsKeys = self.coll.getKeys( dataType='timeseries', tag=self.coll.getObsTag(), location=station, variable=var )
      # create plot
      if clim :
        climVar=clim.get(var, None)
      else :
        climVar=None
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      logScale = True if var in ['kine','vdff','tdff','bed_stress'] else False
      climIsLog = True if var in ['kine','vdff','tdff','bed_stress'] else False
      dia = profilePlot.stackProfileTimeSeriesDC(clabel=varTitle,unit=UNITS[var], clim=climVar, logScale=logScale,climIsLog=climIsLog,**kwargs)
      for modKey in modKeys :
        # add model(s)
        m = self.coll.getSample( **modKey ).copy()
        #m.z = np.tile(m.z.max(axis=0),(m.z.shape[0],1))-m.z # z coords to depth
        for field in m.fieldNames :
          # loop over fields (e.g. u,v in hvel)
          try:
            mf = m.getFields(field)
            dia.addSample( modKey['tag']+field,mf, plotType='contourf',clabel=VARS[field] )
            # add all observations
            for obsKey in obsKeys :
              o = self.coll.getSample( **obsKey ).copy()
              if o.getMetaData('bracket') == 'F' :
                # omit floating platform time series for now
                continue
              #o.z = -o.z # z coords to depth
              dia.addOverlay( modKey['tag']+field, o, plotType='scatter' )
            dia.showColorBar()
            dia.addTitle( '%s %s Dates: %s - %s (PST)'%(station.upper(),modKey['tag'],sT,eT), tag=modKey['tag']+field )
          except:
            pass

      # ----- Save file -----
      file = '_'.join([fPrefix,station,varStr,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()

  def makeProfilerPlots(self, filetype='png', fPrefix='profiler', varPrefix=None, clim=None, **kwargs) :
    """Plots winched profiler data (if available) as trackPlot."""
    kwargs.setdefault( 's', 5 ) # s= dot size in points^2
    comKeys = self.coll.getComparableKeys(dataType='profiler',requireObs=False)
    for entry, obsKey, modKeys in comKeys :
      station,var,msldepth = entry
      keys = []
      if obsKey :
        keys.append( obsKey )
      keys.extend( modKeys )
      varTitle,varStr,sT,eT = self.getAnnotationStrings(var,varPrefix)
      if clim :
        climVar=clim.get(var, None)
      else :
        climVar=None
      kwargs.setdefault( 'ylabel', 'depth below surface' )
      dia = trackPlot.stackTrackPlotDC(clabel=varTitle,unit=UNITS[var],clim=climVar, **kwargs)
      for key in keys :
        sample = self.coll.getSample( **key )
        dia.addPlot( key['tag'] )
        dia.addSample( key['tag'], sample )
        dia.addTitle( '%s %s Dates: %s - %s (PST)'%(station.upper(),key['tag'],sT,eT),tag=key['tag'] )

      # ----- Save file -----
      file = '_'.join([fPrefix,station,varStr,sT,eT])
      self.saveFigure( file, filetype )
      plt.close()
    
  def makeConditionPlot(self, filetype='png', fPrefix='conditions', **kwargs) :
    """Plots wind, river discharge and tidal conditions over the given period"""

    try :
      discharge = self.coll.getDischargeData()
    except Exception as e:
      print e
      discharge = None

    try :
      tidalRange = self.coll.getTidalData()
    except Exception as e:
      print e
      tidalRange = None
    try :
      cui = self.coll.getCoastalUpwellingData()
    except Exception as e:
      print e
      cui = None
    if not discharge and not tidalRange and not cui :
      return

    if not 'title' in kwargs :
      varTitle,varStr,sT,eT = self.getAnnotationStrings('discharge','')
      title = ('Dates: %s - %s (PST)' % (sT, eT))
      kwargs['title'] = ('Dates: %s - %s (PST)' % (sT, eT))
    dia = timeSeriesPlot.timeSeriesStackPlotDC(**kwargs)
    # discharge
    if discharge :
      station = discharge.getMetaData('location')
      ylab = ('%s %s' % ('Discharge',station.upper()))
      dia.addSubplot( 'disch', ylab, 'm3/s' )
      dia.addSample( 'disch', discharge, color='r', label='discharge', ylim=[2000,15000] )
    # tidal range
    if tidalRange :
      station = tidalRange.getMetaData('location')
      ylab = ('%s %s' % ('Tidal range',station.upper()))
      dia.addSubplot( 'tidalr', ylab, 'm' )
      dia.addSample( 'tidalr', tidalRange, color='g', label='tidal range', ylim=[1,4], title=title )
    # cui
    if cui :
      station = cui.getMetaData('location')
      ylab = ('%s %s' % ('Upwelling index',station.upper()))
      dia.addSubplot( 'cui', ylab, '' )
      dia.addSample( 'cui', cui, color='b', label='upwelling index', ylim=[-1300,430] )
    dia.makePlot()

    # ----- Save file -----
    file = '_'.join([fPrefix,sT,eT])
    self.saveFigure( file, filetype )
    plt.close()

  def makeColorsForModels(self):
    """Given models, returns a color dictionary of
       unique colors where the key is modelTag.
    """
    return makeColorsForModels( self.coll )

  def makeMarkersForTaylor(self):
    """Given models, stations and depth, returns a marker and color dictionary
       of unique markers where the key is the tuple (modelTag, station, msldepth).
    """
    return makeMarkersForTaylor(self.coll)

#-------------------------------------------------------------------------------
# Tests 
#-------------------------------------------------------------------------------
def test_Plots():
  """Test class creation and plotting"""
  import time
  from data.stationSet import StationSet
  from data.stationCollection import StationCollection, createDirectory
  startTime = time.time()

  I_PATH='tmp_images'
  M_FILE='/home/workspace/users/lopezj/data/test/process/stationSet_model_2012-04-23-2012-05-20.npz'
  O_FILE='/home/workspace/users/lopezj/data/test/process/stationSet_obs_2012-04-23-2012-05-20.npz'

  createDirectory(I_PATH)
  print 'Loading data'
  oSet = StationSet.loadFromDisk(O_FILE)
  mSet = StationSet.loadFromDisk(M_FILE)
  sc = StationCollection( oSet.startTime, oSet.endTime )
  sc.addStationSet( 'obs', oSet )
  sc.addStationSet( 'model', mSet )
  loadTime = time.time()

  print 'Creating plots'
  tsTime = time.time()
  plot = Plots(I_PATH, sc)
  plot.makeTimeSeries() 
  etsTime = time.time()

  tdTime = time.time()
  plot.makeTaylorDiagrams()
  etdTime = time.time()

  tdvTime = time.time()
  plot.makeTaylorDiagramsVar()
  endTime = time.time()

  print '\n Times \n'
  print 'Elapsed time (Loading data) %s' % (loadTime - startTime)
  print 'Elapsed time (Time series) %s' % (etsTime - tsTime)
  print 'Elapsed time (Taylor diagram) %s' % (etdTime - tdvTime)
  print 'Elapsed time (Taylor diagram (var)) %s' % (endTime - tdvTime)
  print 'Elapsed time (Total) %s' % (endTime - startTime)

#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
if __name__ == '__main__':
  test_Plots()
