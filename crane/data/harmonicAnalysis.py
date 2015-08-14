#!/usr/bin/python

"""
Implementation of harmonic analysis.
Uses tappy for computing the constituents
http://sourceforge.net/apps/mediawiki/tappy/indetap.php?title=Main_Page

Tuomas Karna 2012-09-14
"""

import sys
import numpy as np
import os

import data.timeArray as timeArray
import data.dataContainer as dataContainer
import time
import datetime

def computeTidalConstituents( signal ) :
  """Computes tidal constituents for given signal using TAPPY.
  signal is a dataContainer object containing single time series.
  """
  import data.tappy as tappy
  rayleigh = 1.0
  tap = tappy.tappy(
        outputts=False,
        outputxml='',
        quiet=False,
        debug=False,
        ephemeris=False,
        rayleigh=rayleigh,
        print_vau_table=False,
        missing_data='ignore',
        linear_trend=False,
        remove_extreme=False,
        zero_ts=None,
        filter=None,
        pad_filters=None,
        include_inferred=True,
        )
  
  dates = np.array([ timeArray.epochToDatetime(t) for t in signal.time.asEpoch().array ])
  array = np.squeeze(signal.data)
  
  #print dates
  #print array
  tap.setElevation(dates, array)
  
  if tap.remove_extreme:
      tap.remove_extreme_values()

  package = tap.astronomic(tap.dates)
  (tap.zeta, tap.nu, tap.nup, tap.nupp, tap.kap_p, tap.ii, tap.R, tap.Q, tap.T, tap.jd, tap.s, tap.h, tap.N, tap.p, tap.p1) = package

  if rayleigh:
      ray = float(rayleigh)
  else:
      ray = 1.0
  (tap.speed_dict, tap.key_list) = tap.which_constituents(len(tap.dates),
                                                    package,
                                                    rayleigh_comp = ray)
  if tap.zero_ts:
      # FIX - have to run the constituents package here in order to have
      # filters available , and then run AGAIN later on.
      tap.constituents()
      print len(tap.dates), len(tap.elevation)
      tap.dates_filled, tap.elevation_filled = tap.missing('fill', tap.dates, tap.elevation)
      print len(tap.dates_filled), len(tap.elevation_filled)
      tap.dates, filtered = tap.filters(zero_ts,
                                    tap.dates_filled,
                                    tap.elevation_filled)
      print len(tap.dates), len(filtered)
      tap.elevation = tap.elevation_filled - filtered
      package = tap.astronomic(tap.dates)
      (tap.zeta, tap.nu, tap.nup, tap.nupp, tap.kap_p, tap.ii, tap.R, tap.Q, tap.T, tap.jd, tap.s, tap.h, tap.N, tap.p, tap.p1) = package
      (tap.speed_dict, tap.key_list) = tap.which_constituents(len(tap.dates),
                                                        package,
                                                        rayleigh_comp = ray)
  tap.constituents()
  
  return tap.r, tap.phase

def printConstituents(r, ph) :
  constituents = list(r.keys())
  constituents.sort()
  for k in constituents :
    print '{0:5s} {1:7.4f} {2:6.2f}'.format(k, r[k], ph[k])

class tidalConstituents(object) :
  """A container for harmonic analysis results."""
  TIMEFORMAT = '%Y-%m-%dT%H:%M:%S.%f'
  def __init__( self, description, ampDict, phaseDict, fieldName, startTime, endTime ) :
    """Create new object using the given amplitudes and phases.
    ampDict     -- dictionary of amplitudes (with unit)
    phaseDict   -- dictionary of phases (in degrees) for each constituent
    description -- description string
    fieldName   -- varible name (string)
    startTime   -- start of the analysis period (datetime object)
    endTime     -- end of the analysis period (datetime object)
    """
    if ampDict.keys() != phaseDict.keys() :
      raise Exception( 'ampDict and phaseDict are not compatible' )
    self.constituents = ampDict.keys()
    self.constituents.sort()
    self.amplitude = ampDict
    self.phase = phaseDict
    self.description = description
    self.fieldName = fieldName
    self.startTime = startTime
    self.endTime = endTime
    
  def getConstituent( self, constName ) :
    """Returns amplitude and phase for given constituent.
    Returns (None,None) if constituent is not found."""
    if constName in self.constituents :
      return self.amplitude[constName], self.phase[constName]
    return None, None
    
  def saveToDisk( self, filename, path='' ) :
    """Saves data to disk in numpy array format"""
    if path and not os.path.isdir(path) :
      os.makedirs(path)
      #raise Exception( 'given path does not exist: '+path )
    filename = os.path.join(path,filename)
    print 'writing to '+filename
    kwargs = {}
    kwargs['description'] = self.description
    kwargs['fieldName'] = self.fieldName
    kwargs['startTime'] = self.startTime
    kwargs['endTime'] = self.endTime
    kwargs['constituents'] = self.amplitude.keys()
    kwargs['amplitude'] = np.array(self.amplitude.values())
    kwargs['phase'] = np.array(self.phase.values())    
    np.savez( filename, **kwargs )
  
  @classmethod
  def loadFromDisk( cls, filename ) :
    """Creates a new object from a saved npz file."""
    print 'loading '+filename
    # add npz extension if missing
    namestem,extension = os.path.splitext(filename)
    if not extension :
      filename += '.npz'
    bundle = np.load( filename )
    description = str(bundle['description'])
    fieldName = str(bundle['fieldName'])
    startTime = bundle['startTime'].tolist()
    endTime = bundle['endTime'].tolist()
    amp = list(bundle['amplitude'])
    pha = list(bundle['phase'])
    constituents = bundle['constituents'].tolist()
    ampDict = dict(zip(constituents,amp))
    phaseDict = dict(zip(constituents,pha))
    return cls(description, ampDict, phaseDict, fieldName, startTime, endTime)

  def saveToASCII( self, filename, path='',delimiter=',' ) :
    """Saves data to disk in ascii format"""
    if path and not os.path.isdir(path) :
      os.makedirs(path)
      #raise Exception( 'given path does not exist: '+path )
    filename = os.path.join(path,filename)
    print 'writing to '+filename
    f = open(filename, 'w')
    f.write( self.description+delimiter )
    f.write( self.fieldName+delimiter )
    f.write( self.startTime.strftime(self.TIMEFORMAT)+delimiter )
    f.write( self.endTime.strftime(self.TIMEFORMAT)+'\n' )
    for c in self.amplitude.keys() :
      f.write( c+delimiter )
      f.write( repr(self.amplitude[c])+delimiter )
      f.write( repr(self.phase[c])+'\n' )

  @classmethod
  def loadFromASCII( cls, filename, delimiter=',' ) :
    """Creates a new object from a saved ASCII file."""
    print 'loading '+filename
    f = open( filename, 'r' )
    line = f.readline()
    parts = line.split( delimiter )
    description = parts[0].strip()
    fieldName = parts[1].strip()
    startTime = datetime.datetime.strptime(parts[2].strip(), cls.TIMEFORMAT)
    endTime = datetime.datetime.strptime(parts[3].strip(), cls.TIMEFORMAT)
    tags = []
    amps = []
    phas = []
    for line in f :
      parts = line.split( delimiter )
      tags.append( parts[0] )
      amps.append( float(parts[1]) )
      phas.append( float(parts[2]) )
    ampDict = dict(zip(tags,amps))
    phaseDict = dict(zip(tags,phas))
    return cls(description, ampDict, phaseDict, fieldName, startTime, endTime)
    
  @classmethod
  def computeFromData( cls, dataContainer ) :
    """Computes harmonic analysis for given dataContainer and returns a tidalConstituents object."""
    amp, pha = computeTidalConstituents( dataContainer )
    t = dataContainer.time
    return cls( dataContainer.description, amp, pha, dataContainer.fieldNames[0], t.getDatetime(0), t.getDatetime(-1) )
    
  def printConstituents( self ) :
    print self.description,':', self.fieldName
    print self.startTime,'-', self.endTime
    for c in self.constituents :
      print '{0:5s} {1:7.4f} {2:6.2f}'.format(c, self.amplitude[c], self.phase[c])
      

    
if __name__=='__main__' :
  import datetime

  # generate example data
  startTime = datetime.datetime(2010,1,1,0,0,0)
  endTime   = datetime.datetime(2010,1,13,22,0,0)
  dt = 900.0
  ta = timeArray.generateSimulationTimeArray(startTime,endTime,dt).asEpoch()
  t = ta.array

  T = 44714
  ref = np.sin(2*np.pi*t/T) + 0.95*np.sin(0.95*t/T)
  
  #d0 = dataContainer.dataContainer.fromTimeSeries( 'Observation', ta, ref, ['elev'] )
  dataDir = '/home/tuomas/workspace/cmop/projects/cathlamet_circulation/db22/process/'
  d1 = dataContainer.dataContainer.loadFromDisk(dataDir+'db22_cbnc3_elev_650_2002.npz')
  d0 = d1.timeWindow(datetime.datetime(2002,3,10), datetime.datetime(2002,4,12))
  
  # low level function
  r, phase = computeTidalConstituents( d0 )
  printConstituents(r,phase)
  
  # class object
  tc = tidalConstituents.computeFromData( d0 )
  tc.printConstituents()
  tc.saveToDisk('ha_example')
  
  td = tidalConstituents.loadFromDisk('ha_example')
  td.printConstituents()
  
  if os.path.isfile('ha_example.npz') :
    print 'removing ', 'ha_example.npz'
    os.remove('ha_example.npz')
  
