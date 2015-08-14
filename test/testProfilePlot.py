"""
Tests profile plots.

All test cases should be independent, so that they can be executed
without any other component in the sofware bundle. This ensures that
modules remain independent with clean interfaces.

Tuomas Karna 2013-11-05
"""
import unittest
import numpy as np
from imgTestCaseBase import *
# NOTE this is the only allowed import from the processing lib
from plotting.profilePlot import verticalProfilePlot,verticalProfilePlotDC
from plotting.profilePlot import profileTimeSeries,profileTimeSeriesDC
from plotting.profilePlot import stackProfileTimeSeries
from plotting.profilePlot import stackProfileTimeSeriesDC

class profileTestBase(imgTestCaseBase) :
  """A base class with example profile data sets"""
  def setUp(self):
    imgTestCaseBase.setUp(self)
    """Create all data structures necessary to run the tests."""

    # simple vertical profile data set with no time dependency
    self.staticProfile1 = self.dataSet()
    self.staticProfile1.z = np.linspace(-17.5,0.2,240)
    self.staticProfile1.data = np.sin(self.staticProfile1.z)
    self.staticProfile2 = self.dataSet()
    self.staticProfile2.z = np.linspace(-16.2,0.0,18)
    self.staticProfile2.data = np.sin(0.1*self.staticProfile2.z)
    
    # time dependent faux dataContainer
    self.dc1 = self.dataSet()
    nTime = 245
    nZ = 25
    self.dc1.time = np.linspace(0,80e3,nTime)+951942660.0 # some epoch time
    self.dc1.z = np.zeros((nZ,nTime))
    self.dc1.data = np.zeros((nZ,2,nTime))
    eta = np.sin(self.dc1.time/2/3600)
    h = np.ones_like(eta)*14.5
    for i in range(nTime) :
      self.dc1.z[:,i] = np.linspace(eta[i],h[i],nZ)
    T = np.tile(self.dc1.time,(nZ,1))
    self.dc1.data[:,0,:] = 32*(np.sin(T/2/3600)*np.sin(self.dc1.z/10)+1)/2
    self.dc1.data[:,1,:] = 34*(np.sin(T/2/3600)*self.dc1.z/14.5)
    def detectGaps( *args, **kwargs ):
      return [], [[0,nTime]], self.dc1.time
    self.dc1.detectGaps = detectGaps
    self.dc1.zDependsOnTime=True
    
    # time series at certain depth
    self.ts1 = self.dataSet()
    nTime = 245
    self.ts1.time = np.linspace(0,80e3,nTime)+951942660.0 # some epoch time
    self.ts1.z = np.array([8.2])
    self.ts1.data = np.zeros((1,2,nTime))
    self.ts1.data[:,0,:] = 32*np.sin(self.ts1.time/2/3600)
    self.ts1.data[:,1,:] = 34*np.sin(self.ts1.time/2/3600)
    def detectGaps( *args, **kwargs ):
      return [], [[0,nTime-1]], self.ts1.time
    self.ts1.detectGaps = detectGaps
    self.ts1.zDependsOnTime=False

    self.ts2 = self.dataSet()
    self.ts2.time = np.linspace(0,80e3,nTime)+951942660.0 # some epoch time
    self.ts2.z = np.array([2.2])
    self.ts2.data = np.zeros((1,1,nTime))
    self.ts2.data[:,0,:] = 15*np.sin(self.ts2.time/1/3600)
    def detectGaps( *args, **kwargs ):
      return [], [[0,nTime-1]], self.ts2.time
    self.ts2.detectGaps = detectGaps
    self.ts2.zDependsOnTime=False

#@unittest.skip('skipped')
class testVerticalProfilePlot(profileTestBase):
  """Tests (var,z) profile plots."""

  def testArraysBasicCase(self):
    """Tests basic functionality with array inputs"""
    ax = plt.figure().add_subplot(111)
    dia = verticalProfilePlot(xlabel='Salinity',xunit='psu')
    dia.setAxes(ax)
    dia.addSample(self.staticProfile1.z,self.staticProfile1.data,'sample1',
                  lw=2.0)
    dia.addSample(self.staticProfile2.z,self.staticProfile2.data,'sample2',
                  linestyle='dashed')
    dia.addTitle('arrayBasicCase')
    dia.showLegend(loc=1)
    self.assertCurrentImageByName('arrayBasicCase')
    
  def testDataContainerBasicCase(self):
    """Tests vertical profile plot with dataContainer inputs"""
    ax = plt.figure().add_subplot(111)
    dia = verticalProfilePlotDC(xlabel='Salinity',xunit='psu')
    dia.setAxes(ax)
    dia.addSample(self.dc1,'sample1',iTime=13,iField=0,lw=2.0)
    dia.addSample(self.dc1,'sample2',iTime=13,iField=1,linestyle='dashed')
    dia.addTitle('dcBasicCase')
    dia.showLegend(loc=1)
    self.assertCurrentImageByName('dcBasicCase')

class testProfileTimeSeriesPlot(profileTestBase) :
  """Tests (time,z) time series profile plots."""

  def testArraysBasicCase(self) :
    """Tests basic functionality with array inputs"""
    ax = plt.figure(figsize=(15,7)).add_subplot(111)
    dia = profileTimeSeries(clabel='Salinity',unit='psu',plotType='contourf',
                            clim=[0,34],N=7, ylim=[-16,2.0])
    dia.setAxes( ax )
    dia.addSample( self.dc1.time, -self.dc1.z, self.dc1.data[:,0,:] )
    dia.showColorBar( )
    dia.addSample( self.dc1.time, -self.dc1.z, self.dc1.data[:,1,:],
                  plotType='contour', colors='k', linewidths=1.5, N=5 )
    dia.addTitle('arrayBasicCase')
    self.assertCurrentImageByName('arrayBasicCase')

  def testArraysWithOverlayTimeSeries(self) :
    """Tests with overlayed time series"""
    ax = plt.figure(figsize=(15,7)).add_subplot(111)
    dia = profileTimeSeries(clabel='Salinity',unit='psu',plotType='contourf',
                            clim=[0,34],N=7, ylim=[-16,2.0])
    dia.setAxes( ax )
    dia.addSample( self.dc1.time, -self.dc1.z, self.dc1.data[:,0,:] )
    dia.showColorBar( )
    dia.addOverlay( self.ts1.time, -self.ts1.z, self.ts1.data[:,0,:] )
    dia.addOverlay( self.ts2.time, -self.ts2.z, self.ts2.data[:,0,:] )
    dia.addTitle('arrayWithOverlay')
    self.assertCurrentImageByName('arrayWithOverlay')

  def testDataContainerBasicCase(self) :
    """Tests basic functionality with dataContainer inputs"""
    ax = plt.figure(figsize=(15,7)).add_subplot(111)
    dia = profileTimeSeriesDC(clabel='Salinity',unit='psu',plotType='contourf',
                              clim=[0,34],N=7,invert_yaxis=True)
    dia.setAxes( ax )
    dia.addSample( self.dc1 )
    dia.showColorBar( )
    dia.addSample( self.dc1, iField=1, plotType='contour', colors='k',
                  linewidths=1.5, N=5 )
    dia.addTitle('dcBasicCase')
    self.assertCurrentImageByName('dcBasicCase')

  def testDataContainerWithOverlay(self) :
    """Tests with overlayed time series"""
    ax = plt.figure(figsize=(15,7)).add_subplot(111)
    dia = profileTimeSeriesDC(clabel='Salinity',unit='psu',plotType='contourf',
                              clim=[0,34],N=7,invert_yaxis=True)
    dia.setAxes( ax )
    dia.addSample( self.dc1 )
    dia.showColorBar( )
    dia.addOverlay( self.ts1 )
    dia.addOverlay( self.ts2 )
    dia.addTitle('dcWithOverlay')
    self.assertCurrentImageByName('dcWithOverlay')

class testStackProfileTimeSeriesPlot(profileTestBase) :
  """Tests stacked (time,z) time series profile plots."""

  def testArraysBasicCase(self) :
    """Tests basic functionality with array inputs"""
    dia = stackProfileTimeSeries(clabel='Salinity',unit='psu',plotType='contourf',
                                 clim=[0,34],N=7, ylim=[-16,2.0])
    dia.addSample( 'one', self.dc1.time, -self.dc1.z, self.dc1.data[:,0,:] )
    dia.addSample( 'two', self.dc1.time, -self.dc1.z, self.dc1.data[:,1,:] )
    dia.showColorBar( )
    dia.addSample( 'one', self.dc1.time, -self.dc1.z, self.dc1.data[:,1,:],
                  plotType='contour', colors='k', linewidths=1.5, N=5 )
    dia.addOverlay( 'two', self.ts1.time, -self.ts1.z, self.ts1.data[:,0,:] )
    dia.addOverlay( 'two', self.ts2.time, -self.ts2.z, self.ts2.data[:,0,:] )
    dia.addTitle('arrayBasicCase')
    self.assertCurrentImageByName('arrayBasicCase')

  def testDataContainerBasicCase(self) :
    """Tests basic functionality with dataContainer inputs"""
    dia = stackProfileTimeSeriesDC(clabel='Salinity',unit='psu',plotType='contourf',
                                 clim=[0,34],N=7,invert_yaxis=True)
    dia.addSample( 'one', self.dc1, iField=0 )
    dia.addSample( 'two', self.dc1, iField=1 )
    dia.showColorBar( )
    dia.addSample( 'one', self.dc1, iField=1,
                  plotType='contour', colors='k', linewidths=1.5, N=5 )
    dia.addOverlay( 'two', self.ts1 )
    dia.addOverlay( 'two', self.ts2 )
    dia.addTitle('dcBasicCase')
    self.assertCurrentImageByName('dcBasicCase')

  def testDataContainerLogScale(self) :
    """Tests log plots with dataContainer inputs"""
    dia = stackProfileTimeSeriesDC(clabel='Salinity',unit='psu',plotType='contourf',
                                 logScale=True,N=7,invert_yaxis=True)
    dia.addSample( 'one', self.dc1, iField=0 )
    dia.addSample( 'two', self.dc1, iField=1 )
    dia.showColorBar( )
    dia.addSample( 'one', self.dc1, iField=1,
                  plotType='contour', colors='k', linewidths=1.5, N=5 )
    dia.addOverlay( 'two', self.ts1 )
    dia.addOverlay( 'two', self.ts2 )
    dia.addTitle('dcLogScale')
    self.assertCurrentImageByName('dcLogScale')

if __name__ == '__main__':
  """Run all tests"""
  unittest.main()
