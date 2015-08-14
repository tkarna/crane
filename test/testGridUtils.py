"""
Tests gridUtils

All test cases should be independent, so that they can be executed
without any other component in the sofware bundle. This ensures that
modules remain independent with clean interfaces.

Tuomas Karna 2013-11-05
"""
import unittest
import numpy as np
from data.gridUtils import interpolateInVertical

class testVerticalInterp(unittest.TestCase) :
  """Tests vertical interpolation routine"""
  
  def setUp(self) :
    """Creates data arrays for testing"""
    # register function for comparing numpy arrays
    def allclose(a,b,msg=None) :
      if not np.allclose(a,b,rtol=1e-8,atol=1e-12) :
        msg = 'Arrays do not match %s != %s'%(a,b)
        raise self.failureException(msg)
    self.addTypeEqualityFunc(np.ndarray,allclose)
    # one vertical
    self.Z1 = np.linspace(-10,2,4)[:,None]
    self.V1 = np.linspace(1,2,4)[:,None]
    # three verticals
    self.Z2 = np.vstack( (np.linspace(-10,2.5,5),
                          np.linspace(-1,0,5),
                          np.linspace(-5,5,5)) ).T
    self.V2 = np.vstack( (np.linspace(1,2,5),
                            np.linspace(1,2,5),
                            np.linspace(1,2,5)) ).T
    # NaN-padded z
    self.Z3 = self.Z2.copy()
    self.Z3[0,0] = np.nan
    self.V3 = self.V2.copy()
    # NaN-padded v also
    self.Z4 = self.Z3.copy()
    self.V4 = self.V3.copy()
    self.V4[-1,1] = np.nan
    
  def testZ(self) :
    v,z = interpolateInVertical(self.Z1,self.V1,z=2.0)
    self.assertEqual(v,np.array([2.00]))
    self.assertEqual(z,np.array([2.0]))
    v,z = interpolateInVertical(self.Z1,self.V1,z=-1.0)
    self.assertEqual(v,np.array([1.75]))
    self.assertEqual(z,np.array([-1.0]))
    v,z = interpolateInVertical(self.Z1,self.V1,z=-11.0)
    self.assertEqual(v,np.array([1.00]))
    self.assertEqual(z,np.array([-10.0]))
  
  def testRelToSurf(self) :
    v,z = interpolateInVertical(self.Z1,self.V1,z=0.0,zRelToSurf=True)
    self.assertEqual(v,np.array([2.00]))
    self.assertEqual(z,np.array([2.0]))
    v,z = interpolateInVertical(self.Z1,self.V1,z=3.0,zRelToSurf=True)
    self.assertEqual(v,np.array([1.75]))
    self.assertEqual(z,np.array([-1.0]))
    v,z = interpolateInVertical(self.Z1,self.V1,z=13.0,zRelToSurf=True)
    self.assertEqual(v,np.array([1.00]))
    self.assertEqual(z,np.array([-10.0]))

  def testK(self) :
    v,z = interpolateInVertical(self.Z1,self.V1,k=1)
    self.assertEqual(v,np.array([1.00]))
    self.assertEqual(z,np.array([-10.0]))
    v,z = interpolateInVertical(self.Z1,self.V1,k=-1)
    self.assertEqual(v,np.array([2.00]))
    self.assertEqual(z,np.array([2.0]))
    v,z = interpolateInVertical(self.Z1,self.V1,k=2)
    self.assertEqual(v,np.array([4.0/3]))
    self.assertEqual(z,np.array([-6.0]))
    v,z = interpolateInVertical(self.Z1,self.V1,k=-2)
    self.assertEqual(v,np.array([5.0/3]))
    self.assertEqual(z,np.array([-2.0]))

  def testMultipleArrays(self):
    v,z = interpolateInVertical(self.Z2,self.V2,k=1)
    self.assertEqual(v,np.array([1.,1.,1.]))
    self.assertEqual(z,np.array([-10.0,-1.0,-5.0]))
    v,z = interpolateInVertical(np.flipud(self.Z2),np.flipud(self.V2),k=1)
    self.assertEqual(v,np.array([1.,1.,1.]))
    self.assertEqual(z,np.array([-10.0,-1.0,-5.0]))
    v,z = interpolateInVertical(self.Z2,self.V2,z=0)
    self.assertEqual(v,np.array([1.8,2.,1.5]))
    self.assertEqual(z,np.array([0.0,0.0,0.0]))
    v,z = interpolateInVertical(self.Z2,self.V2,z=-11)
    self.assertEqual(v,np.array([1.,1.,1.]))
    self.assertEqual(z,np.array([-10.0,-1.0,-5.0]))
    
  def testNaNPaddedZ(self):
    v,z = interpolateInVertical(self.Z3,self.V3,z=0)
    self.assertEqual(v,np.array([1.8,2.,1.5]))
    self.assertEqual(z,np.array([0.0,0.0,0.0]))

  def testNaNPaddedZAndV(self):
    v,z = interpolateInVertical(self.Z4,self.V4,z=0)
    self.assertEqual(v,np.array([1.8,1.75,1.5]))
    self.assertEqual(z,np.array([0.0,-0.25,0.0]))

if __name__ == '__main__':
  """Run all tests"""
  unittest.main()
