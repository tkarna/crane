"""
Tests csvStationFile

All test cases should be independent, so that they can be executed
without any other component in the sofware bundle. This ensures that
modules remain independent with clean interfaces.

Tuomas Karna 2013-11-05
"""
import os
import glob
import unittest
from files.csvStationFile import csvStationFile, csvStationFileWithDepth


class testCsvStationFile(unittest.TestCase):
    """Tests tupleCollection object"""

    def setUp(self):
        """Creates a temporary directory where example files are stored"""
        self.tmpDir = 'tmp_' + self.__class__.__name__
        if not os.path.isdir(self.tmpDir):
            os.mkdir(self.tmpDir)
        self.toydata = [('sta1', 3.4, 4.4),
                        ('sta2', 1.4, 1.4),
                        ('sta3', 5.5, 2.0),
                        ('sta4', 5.5, 2.0),
                        ]
        self.goodFileStr = """sta1, 3.40000, 4.40000
sta2, 1.40000, 1.40000
sta3, 5.50000, 2.00000"""
        fid = open(os.path.join(self.tmpDir, 'good.cvs'), 'w')
        fid.write(self.goodFileStr)
        fid.close()
        self.badFileStr = """sta1, 3.40000, 4.40000
sta2, 1.40000, 1.40000
sta3, 5.50000, 2.00000
sta4, 5.50000
"""
        fid = open(os.path.join(self.tmpDir, 'bad.cvs'), 'w')
        fid.write(self.badFileStr)
        fid.close()
        self.badFile2Str = """sta1, 3.40000, 4.40000
sta2, 1.40000, 1.40000
sta3, 5.50000, 2.00000
sta4, 5.50000, ssss
"""
        fid = open(os.path.join(self.tmpDir, 'bad2.cvs'), 'w')
        fid.write(self.badFile2Str)
        fid.close()
        self.depthFileStr = """sta1, 3.40000, 4.40000, -12.90000, z, salt
sta2, 1.40000, 1.40000, 2.20000, depth, temp
sta3, 5.50000, 2.00000, 0.00000, z, elev"""
        fid = open(os.path.join(self.tmpDir, 'depth.cvs'), 'w')
        fid.write(self.depthFileStr)
        fid.close()

    def removeDirTree(self, d):
        """Removes a directory tree d"""
        if os.path.isdir(d):
            for f in glob.glob(os.path.join(d, '*')):
                os.remove(f)
            os.rmdir(d)

    def tearDown(self):
        """Removes temp directory and all its content"""
        if os.path.isdir(self.tmpDir):
            self.removeDirTree(self.tmpDir)

    def testCreation(self):
        sf = csvStationFile()
        self.assertTrue(isinstance(sf, csvStationFile))

    def testAdd(self):
        sf = csvStationFile()
        for loc, x, y in self.toydata:
            sf.addSample(location=loc, x=x, y=y)
        data = sf.getTuples()
        self.assertEqual(data, self.toydata)

    def testWrite(self):
        sf = csvStationFile()
        for loc, x, y in self.toydata[:-1]:
            sf.addSample(location=loc, x=x, y=y)
        sf.writeToDisk(os.path.join(self.tmpDir, 'toy.cvs'))

    def testRead(self):
        sf = csvStationFile()
        sf.readFromFile(os.path.join(self.tmpDir, 'good.cvs'))
        data = sf.getTuples()
        self.assertEqual(data, self.toydata[:-1])

    def testReadFail(self):
        sf = csvStationFile()
        with self.assertRaises(Exception) as e:
            sf.readFromFile(os.path.join(self.tmpDir, 'bad.cvs'))
        self.assertEqual(
            e.exception.args[0],
            'Malformed line [\'sta4\', \'5.50000\'] in file tmp_testCsvStationFile/bad.cvs')

    def testReadFail2(self):
        sf = csvStationFile()
        with self.assertRaises(Exception) as e:
            sf.readFromFile(os.path.join(self.tmpDir, 'bad2.cvs'))
        self.assertEqual(
            e.exception.args[0],
            'Malformed line [\'sta4\', \'5.50000\', \'ssss\'] in file tmp_testCsvStationFile/bad2.cvs')

    def testReadDepthFile(self):
        sf = csvStationFile()
        sf.readFromFile(os.path.join(self.tmpDir, 'depth.cvs'))
        data = sf.getTuples()
        self.assertEqual(data, self.toydata[:-1])

    def testSubSet(self):
        sf = csvStationFile()
        sf.readFromFile(os.path.join(self.tmpDir, 'good.cvs'))
        sf2 = sf.getSubset(x=5.5)
        data = sf2.getTuples()
        self.assertEqual(data, self.toydata[2:-1])

    def testUpdate(self):
        sf = csvStationFile()
        sf.readFromFile(os.path.join(self.tmpDir, 'good.cvs'))
        sf2 = csvStationFile()
        sf2.addSample(('foo', 1.2, 1.2))
        sf.update(sf2)
        data = sf.getTuples()
        expected = self.toydata[:-1]
        expected.append(('foo', 1.2, 1.2))
        self.assertEqual(data, expected)

    def testGetXY(self):
        sf = csvStationFile()
        sf.readFromFile(os.path.join(self.tmpDir, 'good.cvs'))
        x = sf.getX('sta1')
        self.assertEqual(x, 3.4)
        y = sf.getY('sta1')
        self.assertEqual(y, 4.4)


class testCsvStationFileWithDepth(unittest.TestCase):
    """Tests tupleCollection object"""

    def setUp(self):
        """Creates a temporary directory where example files are stored"""
        self.tmpDir = 'tmp_' + self.__class__.__name__
        if not os.path.isdir(self.tmpDir):
            os.mkdir(self.tmpDir)
        self.toydata = [('sta1', 3.4, 4.4, -12.9, 'z', 'salt'),
                        ('sta2', 1.4, 1.4, 2.2, 'depth', 'temp'),
                        ('sta3', 5.5, 2.0, 0, 'z', 'elev'),
                        ('sta4', 5.5, 2.0, 2.9, 'z', None),
                        ]
        self.goodFileStr = """sta1, 3.40000, 4.40000, -12.90000, z, salt
sta2, 1.40000, 1.40000, 2.20000, depth, temp
sta3, 5.50000, 2.00000, 0.00000, z, elev"""
        fid = open(os.path.join(self.tmpDir, 'good.cvs'), 'w')
        fid.write(self.goodFileStr)
        fid.close()
        self.badFileStr = """sta1, 3.40000, 4.40000, -12.90000, z, salt
sta2, 1.40000, 1.40000, 2.20000, depth, temp
sta3, 5.50000, 2.00000, 0.00000, z, elev
sta4, 5.50000, 2.00000, 2.90000, z
"""
        fid = open(os.path.join(self.tmpDir, 'bad.cvs'), 'w')
        fid.write(self.badFileStr)
        fid.close()
        self.badFile2Str = """sta1, 3.40000, 4.40000, -12.90000, z, salt
sta2, 1.40000, 1.40000, 2.20000, depth, temp
sta3, 5.50000, 2.00000, 0.00000, z, elev
sta4, 5.50000, ssss, 2.90000, z, elev
"""
        fid = open(os.path.join(self.tmpDir, 'bad2.cvs'), 'w')
        fid.write(self.badFile2Str)
        fid.close()

    def removeDirTree(self, d):
        """Removes a directory tree d"""
        if os.path.isdir(d):
            for f in glob.glob(os.path.join(d, '*')):
                os.remove(f)
            os.rmdir(d)

    def tearDown(self):
        """Removes temp directory and all its content"""
        if os.path.isdir(self.tmpDir):
            self.removeDirTree(self.tmpDir)

    def testCreation(self):
        sf = csvStationFileWithDepth()
        self.assertTrue(isinstance(sf, csvStationFileWithDepth))

    def testAdd(self):
        sf = csvStationFileWithDepth()
        for loc, x, y, z, zType, var in self.toydata:
            sf.addSample(
                location=loc,
                x=x,
                y=y,
                z=z,
                zType=zType,
                variable=var)
        data = sf.getTuples()
        self.assertEqual(data, self.toydata)

    def testWrite(self):
        sf = csvStationFileWithDepth()
        for loc, x, y, z, zType, var in self.toydata[:-1]:
            sf.addSample(
                location=loc,
                x=x,
                y=y,
                z=z,
                zType=zType,
                variable=var)
        sf.writeToDisk(os.path.join(self.tmpDir, 'toy.cvs'))

    def testRead(self):
        sf = csvStationFileWithDepth()
        sf.readFromFile(os.path.join(self.tmpDir, 'good.cvs'))
        data = sf.getTuples()
        self.assertEqual(data, self.toydata[:-1])

    def testReadFail(self):
        sf = csvStationFileWithDepth()
        with self.assertRaises(Exception) as e:
            sf.readFromFile(os.path.join(self.tmpDir, 'bad.cvs'))
        self.assertEqual(
            e.exception.args[0],
            'Malformed line [\'sta4\', \'5.50000\', \'2.00000\', \'2.90000\', \'z\'] in file tmp_testCsvStationFileWithDepth/bad.cvs')

    def testReadFail2(self):
        sf = csvStationFileWithDepth()
        with self.assertRaises(Exception) as e:
            sf.readFromFile(os.path.join(self.tmpDir, 'bad2.cvs'))
        self.assertEqual(
            e.exception.args[0],
            'Malformed line [\'sta4\', \'5.50000\', \'ssss\', \'2.90000\', \'z\', \'elev\'] in file tmp_testCsvStationFileWithDepth/bad2.cvs')

    def testSubSet(self):
        sf = csvStationFileWithDepth()
        sf.readFromFile(os.path.join(self.tmpDir, 'good.cvs'))
        sf2 = sf.getSubset(x=5.5)
        data = sf2.getTuples()
        self.assertEqual(data, self.toydata[2:-1])

    def testUpdate(self):
        sf = csvStationFileWithDepth()
        sf.readFromFile(os.path.join(self.tmpDir, 'good.cvs'))
        sf2 = csvStationFileWithDepth()
        sf2.addSample(('foo', 1.2, 1.2, 1.2, 'z', 'salt'))
        sf.update(sf2)
        data = sf.getTuples()
        expected = self.toydata[:-1]
        expected.append(('foo', 1.2, 1.2, 1.2, 'z', 'salt'))
        self.assertEqual(data, expected)

    def testGetXY(self):
        sf = csvStationFileWithDepth()
        sf.readFromFile(os.path.join(self.tmpDir, 'good.cvs'))
        x = sf.getX('sta1')
        self.assertEqual(x, 3.4)
        y = sf.getY('sta1')
        self.assertEqual(y, 4.4)
        z = sf.getZ('sta1')
        self.assertEqual(z, -12.9)

if __name__ == '__main__':
    """Run all tests"""
    unittest.main()
