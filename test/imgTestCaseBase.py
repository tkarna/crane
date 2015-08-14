"""
A base class for testing image generation.
Uses matplotlib compare_images routine to check similarity of png images

Tuomas Karna 2013-11-05
"""
import os
import glob
import shutil
import unittest
import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images

class ImageComparisonFailure(Exception):
  """Dummy exception class for image comparison failures."""
  pass

class imgTestCaseBase(unittest.TestCase):
  """Base class for all test that need to generate, store and compare images"""
  def setUp(self):
    """Makes temp directory where all images will be stored."""
    self.tmpDir = 'tmp_'+self.__class__.__name__
    self.assertDir = os.path.join('baseline_images',self.__class__.__name__)
    if not os.path.isdir(self.tmpDir):
      os.mkdir(self.tmpDir)
    if not os.path.isdir(self.assertDir):
      os.makedirs(self.assertDir)
    
  def removeDirTree( self, d ):
    """Removes a directory tree d"""
    if os.path.isdir( d ):
      for f in glob.glob(os.path.join(d,'*')):
        os.remove(f)
      os.rmdir(d)

  def tearDown(self):
    """Removes temp directory and all its content"""
    if os.path.isdir(self.tmpDir):
      self.removeDirTree(self.tmpDir)

  def getNewFigFilename(self, figname, ext='png'):
    """Generates a path to a png file in the temp directory"""
    return os.path.join( self.tmpDir, figname+'.'+ext )

  def getExpectedFigFilename(self, figname, ext='png'):
    """Generates a path to a png file in the temp directory"""
    return os.path.join( self.assertDir, self.__class__.__name__+'_'+figname+'.'+ext )

  def assertCurrentImageByName(self, figname, tol=0.0001):
    """Stores current figure in temporary location and compares to the
    expected image file.
    
    If the expected image file does not exist (e.g. in the first pass), the
    current figure is stored as the correct image and IOError is raised to
    warn the user. The test will pass next time."""
    actual = self.getNewFigFilename(figname)
    plt.savefig(actual)
    expected = self.getExpectedFigFilename(figname)
    if not os.path.isfile(expected):
      shutil.copy(actual,expected)
      msg = 'Sample image to compare against does not exist. Copying currently\n'
      msg += 'generated image to location:\n'
      msg += '\''+expected+'\''
      msg += '\nThis image will be used as the correct sample from now on.'
      msg += '\nPlease check the image and delete it if not acceptable.'
      raise IOError(msg)
    self.assertImage(actual,expected,tol)
    plt.close('all')

  def assertImage(self,expected,actual,tol=0.0001):
    """Tests wether images are equal (up to a tolerance)
    
    Parameters
    ----------
    expected : str
        Filename of the correct image
    actual : str
        Filename of the actual generated image
    tol : float, optional
        Tolerance
    
    Raises
    ------
    ImageComparisonFailure
        If images differ more than the given tolerance
    """
    if not os.path.isfile(expected):
      raise IOError('File '+expected+' does not exist')
    if not os.path.isfile(actual):
      raise IOError('File '+expected+' does not exist')
    err = compare_images(expected,actual,tol)
    if err :
      raise ImageComparisonFailure( err )

  class dataSet(object):
    """Tiny data container object with no members"""
    pass
