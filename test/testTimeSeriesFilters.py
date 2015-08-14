"""
Tests timeSeriesFilters

All test cases should be independent, so that they can be executed
without any other component in the sofware bundle. This ensures that
modules remain independent with clean interfaces.

Tuomas Karna 2013-11-07 
"""
import numpy as np
import unittest
from data.timeSeriesFilters import *

class motherTestCase(unittest.TestCase) :
  """Generic test case class with some common helper functions"""
  class dataSet(object):
    """Tiny data container object with no members"""
    pass

  def setUp(self) :
    def allclose(a,b,msg=None) :
      if not np.allclose(a,b,rtol=1e-8,atol=1e-12) :
        msg = 'Arrays do not match %s != %s'%(a,b)
        raise self.failureException(msg)
    self.addTypeEqualityFunc(np.ndarray,allclose)
    
class testTimeSeriesFilters(motherTestCase) :
  """Tests functions in timeSeriesFilters"""

  def setUp(self) :
    """Build test arrays"""
    # register function for comparing numpy arrays
    motherTestCase.setUp(self)
    # arrays
    self.time1 = np.linspace(0,1,11)
    self.vals1 = self.time1*2.0
    # faux dataContainer
    self.dc1 = self.dataSet()
    self.dc1.time = self.time1
    self.dc1.data = self.vals1[None,None,:]
    def detectGaps( *args, **kwargs ):
      return [], np.array([[0,len(self.dc1.time)]]), self.dc1.time
    self.dc1.detectGaps = detectGaps
    self.dc1.x = 0.0
    self.dc1.y = np.linspace(0,1,11)[None,:]
    self.dc1.z = 0.0
    self.dc1.xDependsOnTime=False
    self.dc1.yDependsOnTime=True
    self.dc1.zDependsOnTime=False
    def getMetaData() :
      return {}
    self.dc1.getMetaData = getMetaData
    self.dc1.fieldNames = ['field']
    self.dc1.coordSys = ''

  def testArrayMean(self) :
    t,x = computeRunningMean(self.time1,self.vals1,T=0.2-1e-12)
    self.assertEqual(t,np.linspace(0.1,0.9,9))
    self.assertEqual(x,np.linspace(0.2,1.8,9))

  def testArrayMax(self) :
    t,x = computeRunningMax(self.time1,self.vals1,T=0.2-1e-12)
    self.assertEqual(t,np.linspace(0.1,0.9,9))
    self.assertEqual(x,np.linspace(0.4,2.0,9))

  def testArrayMin(self) :
    t,x = computeRunningMin(self.time1,self.vals1,T=0.2-1e-12)
    self.assertEqual(t,np.linspace(0.1,0.9,9))
    self.assertEqual(x,np.linspace(0.0,1.6,9))

  def testArrayRange(self) :
    t,x = computeRunningRange(self.time1,self.vals1,T=0.2-1e-12)
    self.assertEqual(t,np.linspace(0.1,0.9,9))
    self.assertEqual(x,np.ones((9,))*0.4)

  def testDCMean(self) :
    dc = runningMean(self.dc1,T=0.2-1e-12)
    self.assertEqual(dc.time.array,np.linspace(0.1,0.9,9))
    self.assertEqual(dc.data[0,0,:],np.linspace(0.2,1.8,9))
    self.assertEqual(dc.y[0,:],np.linspace(0.1,0.9,9))

  def testDCMax(self) :
    dc = runningMax(self.dc1,T=0.2-1e-12)
    self.assertEqual(dc.time.array,np.linspace(0.1,0.9,9))
    self.assertEqual(dc.data[0,0,:],np.linspace(0.4,2.0,9))
    self.assertEqual(dc.y[0,:],np.linspace(0.1,0.9,9))

  def testDCMin(self) :
    dc = runningMin(self.dc1,T=0.2-1e-12)
    self.assertEqual(dc.time.array,np.linspace(0.1,0.9,9))
    self.assertEqual(dc.data[0,0,:],np.linspace(0.0,1.6,9))
    self.assertEqual(dc.y[0,:],np.linspace(0.1,0.9,9))
  
  def testDCRange(self) :
    dc = runningRange(self.dc1,T=0.2-1e-12)
    self.assertEqual(dc.time.array,np.linspace(0.1,0.9,9))
    self.assertEqual(dc.data[0,0,:],np.ones((9,))*0.4)
    self.assertEqual(dc.y[0,:],np.linspace(0.1,0.9,9))

if __name__ == '__main__':
  """Run all tests"""
  unittest.main()
