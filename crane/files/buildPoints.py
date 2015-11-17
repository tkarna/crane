"""
Handles reading and writing of build point files. (.bp) used for many
pre- and post-processing applications.
"""
#-------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------
import os
import numpy as np

#-------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Classes and functions
#-------------------------------------------------------------------------


class BuildPoint(object):
    """Handles reading and writing of build point files."""

    def __init__(self):
        """Just returns a blank object"""
        self.path = ''
        self.description = ''
        self.coordSys = ''
        self.vDatum = ''
        self.points = []

    def regenerateIndices(self):
        """Rebuilds point running index"""
        self.points[:, 0] = np.arange(self.points.shape[0]) + 1

    def sortByX(self):
        """Sort points by x coordinate"""
        ix = np.argsort(self.points[:, 1])
        self.points = self.points[ix, :]
        self.regenerateIndices()

    def reverseOrder(self):
        self.points = np.flipud(self.points)
        self.regenerateIndices()

    def getX(self):
        return self.points[:, 1]

    def getY(self):
        return self.points[:, 2]

    def setXYZ(self, x, y, z):
        nPoints = len(x)
        self.points = np.zeros((nPoints, 4))
        self.points[:, 1] = x
        self.points[:, 2] = y
        self.points[:, 3] = z
        self.points[:, 0] = np.arange(nPoints) + 1

    def readFileFromDisk(self, file):
        """Reads a build point file.

        Args:
          file - String of path to a .bp file
        """
        if not os.path.isfile(file):
            raise Exception('File does not exist: ' + file)
        self.file = file

        f = open(file, 'r')
        self.description = f.readline().rstrip('\n')
        nPoints = int(f.readline().rstrip('\n'))
        self.points = np.empty([nPoints, 4])
        for p in xrange(nPoints):
            tmp = f.readline().split()
            self.points[p, 0] = int(tmp[0])
            self.points[p, 1] = float(tmp[1])
            self.points[p, 2] = float(tmp[2])
            self.points[p, 3] = float(tmp[3])
        f.close()
        self.idCoordSys()

    def writeFileToDisk(self, file, overwrite=False):
        """Writes object to *.bp file.

        Args:
          file -- String of the path to build point file
        """
        if os.path.isfile(file) and not overwrite:
            raise Exception('File already exists: ' + file)

        f = open(file, 'w')
        f.write(self.description + '\n')
        tmp = str(self.points.shape[0])
        f.write(tmp + '\n')
        for p in xrange(self.points.shape[0]):
            pStr = '%d %.6f %.6f %f\n' % (p + 1, self.points[p, 1],
                                          self.points[p, 2], self.points[p, 3])
            f.write(pStr)
        f.close()

    def idCoordSys(self):
        """Ids coordinate system as 'spcs' or 'latlon' based on x,y values."""
        if isinstance(self.points, list):
            raise Exception('Must load build points before coordinate '
                            'system can be identified.')

        if abs(self.points[0, 0]) > 180 and abs(self.points[0, 1]) > 180:
            self.coordSys = 'spcs'
        elif abs(self.points[0, 0]) <= 180 and abs(self.points[0, 1]) <= 180:
            self.coordSys = 'latlon'

if __name__ == '__main__':
    # Test class
    print 'Instantiating class\n'
    bP = BuildPoint()

    print 'Read build point file\n'
    BP_FILE = '/home/workspace/users/lopezj/data/test/fg.bp'
    TEST_WRITE = '/home/workspace/users/lopezj/data/test/test_mouth1.bp'
    bP.readFileFromDisk(BP_FILE)
    print 'Description: ' + bP.description
    print 'File: ' + bP.file
    print 'Coordinate system: ' + bP.coordSys
    print 'Vertical datum: ' + bP.vDatum
    print 'Number of points: ' + str(bP.points.shape[0])

    print 'Write new file: ' + TEST_WRITE
    bP.writeFileToDisk(TEST_WRITE)

    print '\nInstantiating new class\n'
    bP2 = BuildPoint()

    print 'Read new build point file\n'
    bP2.readFileFromDisk(TEST_WRITE)
    print 'Description: ' + bP2.description
    print 'File: ' + bP2.file
    print 'Coordinate system: ' + bP2.coordSys
    print 'Vertical datum: ' + bP2.vDatum
    print 'Number of points: ' + str(bP2.points.shape[0])

    print '\nTest for equivalence\n'
    good = True
    if bP2.description != bP.description:
        good = False
        print 'Fail description'
    if bP2.coordSys != bP.coordSys:
        good = False
        print 'Fail coordinate system'
    if bP2.vDatum != bP.vDatum:
        good = False
        print 'Fail vertical datum'
    if bP2.points.shape[0] != bP.points.shape[0]:
        good = False
        print 'Fail number of build points'
    for p in xrange(bP2.points.shape[0]):
        if bP2.points[p, 0] != bP.points[p, 0]:
            good = False
            print 'Fail x -- bP: %f bP2: %f' % (bP.points[p, 0], bP2.points[p, 0])
        if bP2.points[p, 1] != bP.points[p, 1]:
            good = False
            print 'Fail y -- bP: %f bP2: %f' % (bP.points[p, 1], bP2.points[p, 1])
        if bP2.points[p, 2] != bP.points[p, 2]:
            good = False
            print 'Fail z -- bP: %f bP2: %f' % (bP.points[p, 2], bP2.points[p, 2])
        if good:
            print 'Equivalence test passed'
        else:
            print 'Equivalence test failed'

    # os.remove(TEST_WRITE)
