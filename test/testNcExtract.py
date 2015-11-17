import pstats
import cProfile
import numpy as np
import data.ncExtract as NCE
from data.timeArray import *

import unittest

# path='/home/tuomas/workspace/cmop/projects/turb_tests/cre_open_channel/real_bath/fluxTest/combined'
path = '/home/workspace/ccalmr53/karnat/projects/turb_tests/fluxtest/fluxTest/combined/'

# all good
#staX = np.array([10.0,4600.0])
#staY = np.array([70.0,60.0])
#staZ = np.array([-5.2,-3.5])
#staT = np.array([300.0,92600.0])+1333700100.0
#staNames = ['foo','bar']

# one bad
#staX = np.array([10.0,4600.0,5000.0])
#staY = np.array([70.0,60.0,200.0])
#staZ = np.array([-5.2,-3.5,-2.5])
#staT = np.array([300.0,92600.0,100000.0])+1333700100.0
#staNames = ['foo','bar','bad']

# one bad, two near each other
staX = np.array([33.4, 4600.0, 4600.0, 5000.0])
staY = np.array([50.0, 60.0, 60.5, 200.0])
staZ = np.array([-5.2, -3.5, -3.5, -2.5])
staT = np.array([300.0, 92600.0, 92601.0, 100000.0]) + 1333700100.0
staNames = ['foo', 'bar', 'bar2', 'bad']
nGoodSta = 3

st = datetime.datetime(2012, 4, 5)
et = datetime.datetime(2012, 4, 9)

verbose = False


class TestTimeSeriesExtraction(unittest.TestCase):

    def testElev(self):
        se = NCE.selfeExtract(path, var='elev', verbose=verbose)
        dcs = se.extractTimeSeries(
            st, et, 'elev', staX, staY, staNames, staZ=staZ)
        self.assertEqual(len(dcs), nGoodSta)

    def testSalt(self):
        se = NCE.selfeExtract(path, var='salt', verbose=verbose)
        dcs = se.extractTimeSeries(
            st, et, 'salt', staX, staY, staNames, staZ=staZ)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvel(self):
        se = NCE.selfeExtract(path, var='hvel', verbose=verbose)
        dcs = se.extractTimeSeries(
            st, et, 'hvel', staX, staY, staNames, staZ=staZ)
        self.assertEqual(len(dcs), nGoodSta)

    def testSalt70(self):
        se = NCE.selfeExtract(
            path,
            var='salt',
            verbose=verbose,
            fileTypeStr='70')
        dcs = se.extractTimeSeries(
            st, et, 'salt', staX, staY, staNames, staZ=staZ)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvel67(self):
        se = NCE.selfeExtract(
            path,
            var='hvel',
            verbose=verbose,
            fileTypeStr='67')
        dcs = se.extractTimeSeries(
            st, et, 'hvel', staX, staY, staNames, staZ=staZ)
        self.assertEqual(len(dcs), nGoodSta)


class TestProfileExtraction(unittest.TestCase):

    def testElevFail(self):
        se = NCE.selfeExtract(path, var='elev', verbose=verbose)
        with self.assertRaises(Exception) as cm:
            dcs = se.extractVerticalProfile(
                st, et, 'elev', staX, staY, staNames)
        targetError = 'This is a 2D variable, cannot extract vertical profile: elev'
        self.assertEqual(cm.exception.message, targetError,
                         'Got wrong error message: ' + cm.exception.message)

    def testSalt(self):
        se = NCE.selfeExtract(path, var='salt', verbose=verbose)
        dcs = se.extractVerticalProfile(st, et, 'salt', staX, staY, staNames)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvel(self):
        se = NCE.selfeExtract(path, var='hvel', verbose=verbose)
        dcs = se.extractVerticalProfile(st, et, 'hvel', staX, staY, staNames)
        self.assertEqual(len(dcs), nGoodSta)

    def testSalt70(self):
        se = NCE.selfeExtract(
            path,
            var='salt',
            verbose=verbose,
            fileTypeStr='70')
        dcs = se.extractVerticalProfile(st, et, 'salt', staX, staY, staNames)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvel67(self):
        se = NCE.selfeExtract(
            path,
            var='hvel',
            verbose=verbose,
            fileTypeStr='67')
        dcs = se.extractVerticalProfile(st, et, 'hvel', staX, staY, staNames)
        self.assertEqual(len(dcs), nGoodSta)


class TestTransectExtraction(unittest.TestCase):

    def testElev(self):
        se = NCE.selfeExtract(path, var='elev', verbose=verbose)
        dc = se.extractTransect(st, et, 'elev', staX, staY, 'trans')
        self.assertEqual(dc.getMetaData('variable'), 'elev')

    def testSalt(self):
        se = NCE.selfeExtract(path, var='salt', verbose=verbose)
        dc = se.extractTransect(st, et, 'salt', staX, staY, 'trans')
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvel(self):
        se = NCE.selfeExtract(path, var='hvel', verbose=verbose)
        dc = se.extractTransect(st, et, 'hvel', staX, staY, 'trans')
        self.assertEqual(dc.getMetaData('variable'), 'hvel')

    def testSalt70(self):
        se = NCE.selfeExtract(
            path,
            var='salt',
            verbose=verbose,
            fileTypeStr='70')
        dc = se.extractTransect(st, et, 'salt', staX, staY, 'trans')
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvel70(self):
        se = NCE.selfeExtract(
            path,
            var='hvel',
            verbose=verbose,
            fileTypeStr='67')
        dc = se.extractTransect(st, et, 'hvel', staX, staY, 'trans')
        self.assertEqual(dc.getMetaData('variable'), 'hvel')


class TestTrackExtraction(unittest.TestCase):

    def testElev(self):
        se = NCE.selfeExtract(path, var='elev', verbose=verbose)
        dc = se.extractTrack('elev', staX, staY, staZ, staT, 'foo')
        self.assertEqual(dc.getMetaData('variable'), 'elev')

    def testSalt(self):
        se = NCE.selfeExtract(path, var='salt', verbose=verbose)
        dc = se.extractTrack('salt', staX, staY, staZ, staT, 'foo')
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvel(self):
        se = NCE.selfeExtract(path, var='hvel', verbose=verbose)
        dc = se.extractTrack('hvel', staX, staY, staZ, staT, 'foo')
        self.assertEqual(dc.getMetaData('variable'), 'hvel')

    def testSalt70(self):
        se = NCE.selfeExtract(
            path,
            var='salt',
            verbose=verbose,
            fileTypeStr='70')
        dc = se.extractTrack('salt', staX, staY, staZ, staT, 'foo')
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvel67(self):
        se = NCE.selfeExtract(
            path,
            var='hvel',
            verbose=verbose,
            fileTypeStr='67')
        dc = se.extractTrack('hvel', staX, staY, staZ, staT, 'foo')
        self.assertEqual(dc.getMetaData('variable'), 'hvel')


class TestSlabExtraction(unittest.TestCase):

    def testElevK(self):
        se = NCE.selfeExtract(path, var='elev', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'elev', k=1)
        self.assertEqual(dc.getMetaData('variable'), 'elev')

    def testSaltK(self):
        se = NCE.selfeExtract(path, var='salt', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'salt', k=1)
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvelK(self):
        se = NCE.selfeExtract(path, var='hvel', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'hvel', k=1)
        self.assertEqual(dc.getMetaData('variable'), 'hvel')

    def testSaltK70(self):
        se = NCE.selfeExtract(
            path,
            var='salt',
            verbose=verbose,
            fileTypeStr='70')
        dc = se.extractSlab(st, et, 'slab', 'salt', k=1)
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvelK67(self):
        se = NCE.selfeExtract(
            path,
            var='hvel',
            verbose=verbose,
            fileTypeStr='67')
        dc = se.extractSlab(st, et, 'slab', 'hvel', k=1)
        self.assertEqual(dc.getMetaData('variable'), 'hvel')

    def testElevKneg(self):
        se = NCE.selfeExtract(path, var='elev', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'elev', k=-1)
        self.assertEqual(dc.getMetaData('variable'), 'elev')

    def testSaltKneg(self):
        se = NCE.selfeExtract(path, var='salt', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'salt', k=-1)
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvelKneg(self):
        se = NCE.selfeExtract(path, var='hvel', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'hvel', z=-1)
        self.assertEqual(dc.getMetaData('variable'), 'hvel')

    def testSaltKneg70(self):
        se = NCE.selfeExtract(
            path,
            var='salt',
            verbose=verbose,
            fileTypeStr='70')
        dc = se.extractSlab(st, et, 'slab', 'salt', k=-1)
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvelKneg67(self):
        se = NCE.selfeExtract(
            path,
            var='hvel',
            verbose=verbose,
            fileTypeStr='67')
        dc = se.extractSlab(st, et, 'slab', 'hvel', z=-1)
        self.assertEqual(dc.getMetaData('variable'), 'hvel')

    def testElevZ(self):
        se = NCE.selfeExtract(path, var='elev', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'elev', z=-2.5)
        self.assertEqual(dc.getMetaData('variable'), 'elev')

    def testSaltZ(self):
        se = NCE.selfeExtract(path, var='salt', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'salt', z=-2.5)
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvelZ(self):
        se = NCE.selfeExtract(path, var='hvel', verbose=verbose)
        dc = se.extractSlab(st, et, 'slab', 'hvel', z=-2.5)
        self.assertEqual(dc.getMetaData('variable'), 'hvel')

    def testSaltZ70(self):
        se = NCE.selfeExtract(
            path,
            var='salt',
            verbose=verbose,
            fileTypeStr='70')
        dc = se.extractSlab(st, et, 'slab', 'salt', z=-2.5)
        self.assertEqual(dc.getMetaData('variable'), 'salt')

    def testHvelZ67(self):
        se = NCE.selfeExtract(
            path,
            var='hvel',
            verbose=verbose,
            fileTypeStr='67')
        dc = se.extractSlab(st, et, 'slab', 'hvel', z=-2.5)
        self.assertEqual(dc.getMetaData('variable'), 'hvel')

#-------------------------------------------------------------------------
# High-level interface
#-------------------------------------------------------------------------


class TestHLTimeSeriesExtractionXYZ(unittest.TestCase):

    def testElevZ(self):
        dcs = NCE.extractForXYZ(
            path,
            'elev',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testElevZRelSurf(self):
        dcs = NCE.extractForXYZ(
            path,
            'elev',
            st,
            et,
            staX,
            staY,
            z=-staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=True,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testSaltZ(self):
        dcs = NCE.extractForXYZ(
            path,
            'salt',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testSaltZRelSurf(self):
        dcs = NCE.extractForXYZ(
            path,
            'salt',
            st,
            et,
            staX,
            staY,
            z=-staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=True,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testSalt70Z(self):
        dcs = NCE.extractForXYZ(
            path,
            'salt.70',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testSalt70ZRelSurf(self):
        dcs = NCE.extractForXYZ(
            path,
            'salt.70',
            st,
            et,
            staX,
            staY,
            z=-staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=True,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvelZ(self):
        dcs = NCE.extractForXYZ(
            path,
            'hvel',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvelZRelSurf(self):
        dcs = NCE.extractForXYZ(
            path,
            'hvel',
            st,
            et,
            staX,
            staY,
            z=-staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=True,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvelZ(self):
        dcs = NCE.extractForXYZ(
            path,
            'hvel.67',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvelZRelSurf(self):
        dcs = NCE.extractForXYZ(
            path,
            'hvel.67',
            st,
            et,
            staX,
            staY,
            z=-staZ,
            stationNames=staNames,
            profile=False,
            zRelToSurf=True,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)


class TestHLProfileExtractionXYZ(unittest.TestCase):

    def testElevFail(self):
        dcs = NCE.extractForXYZ(
            path,
            'elev',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=True,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), 0)

    def testSalt(self):
        dcs = NCE.extractForXYZ(
            path,
            'salt',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=True,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testSalt70(self):
        dcs = NCE.extractForXYZ(
            path,
            'salt.70',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=True,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvel(self):
        dcs = NCE.extractForXYZ(
            path,
            'hvel',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=True,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)

    def testHvel67(self):
        dcs = NCE.extractForXYZ(
            path,
            'hvel.67',
            st,
            et,
            staX,
            staY,
            z=staZ,
            stationNames=staNames,
            profile=True,
            zRelToSurf=False,
            verbose=verbose)
        self.assertEqual(len(dcs), nGoodSta)


class TestHLTransectExtractionForCoords(unittest.TestCase):

    def testMultiple(self):
        dcs = NCE.extractTransectForCoords(staX, staY, path, ['salt', 'hvel'],
                                           st, et, 'foo', verbose=verbose)
        self.assertEqual(len(dcs), 2)

    def testMultipleAltDiscretization(self):
        dcs = NCE.extractTransectForCoords(
            staX, staY, path, ['salt.70', 'hvel.67'],
            st, et, 'foo', verbose=verbose)
        self.assertEqual(len(dcs), 2)


class TestHLSlabExtraction(unittest.TestCase):

    def testMultipleK(self):
        dcs = NCE.extractSlabForLevel(
            path, ['elev', 'salt', 'hvel'],
            st, et, 'foo', k=1, verbose=True)
        self.assertEqual(len(dcs), 3)

    def testSaltK(self):
        dcs = NCE.extractSlabForLevel(path, ['salt'], st, et, 'foo',
                                      k=1, verbose=True)
        self.assertEqual(len(dcs), 1)

    def testSaltKAlt(self):
        dcs = NCE.extractSlabForLevel(path, ['salt.70'], st, et, 'foo',
                                      k=1, verbose=True)
        self.assertEqual(len(dcs), 1)

    def testHvelK(self):
        dcs = NCE.extractSlabForLevel(path, ['hvel'], st, et, 'foo',
                                      k=1, verbose=True)
        self.assertEqual(len(dcs), 1)

    def testHvelKAlt(self):
        dcs = NCE.extractSlabForLevel(path, ['hvel.67'], st, et, 'foo',
                                      k=1, verbose=True)
        self.assertEqual(len(dcs), 1)

    def testSaltKneg5(self):
        dcs = NCE.extractSlabForLevel(path, ['salt'], st, et, 'foo',
                                      k=-5, verbose=True)
        self.assertEqual(len(dcs), 1)

    def testSaltKneg5Alt(self):
        dcs = NCE.extractSlabForLevel(path, ['salt.70'], st, et, 'foo',
                                      k=-5, verbose=True)
        self.assertEqual(len(dcs), 1)

    def testHvelKneg5(self):
        dcs = NCE.extractSlabForLevel(path, ['hvel'], st, et, 'foo',
                                      k=-5, verbose=True)
        self.assertEqual(len(dcs), 1)

    def testHvelKneg5Alt(self):
        dcs = NCE.extractSlabForLevel(path, ['hvel.67'], st, et, 'foo',
                                      k=-5, verbose=True)
        self.assertEqual(len(dcs), 1)

    def testMultipleKAlt(self):
        dcs = NCE.extractSlabForLevel(
            path, ['salt.70', 'hvel.67'],
            st, et, 'foo', k=1, verbose=True)
        self.assertEqual(len(dcs), 2)

    def testMultipleKneg(self):
        dcs = NCE.extractSlabForLevel(
            path, ['elev', 'salt', 'hvel'],
            st, et, 'foo', k=-1, verbose=True)
        self.assertEqual(len(dcs), 3)

    def testMultipleZ(self):
        dcs = NCE.extractSlabForLevel(
            path, ['elev', 'salt', 'hvel'],
            st, et, 'foo', z=-5.5, verbose=True)
        self.assertEqual(len(dcs), 3)

    def testMultipleZRelSurf(self):
        dcs = NCE.extractSlabForLevel(
            path, ['elev', 'salt', 'hvel'],
            st, et, 'foo', z=5.5, zRelToSurf=True, verbose=True)
        self.assertEqual(len(dcs), 3)

if __name__ == '__main__':
    # run all tests
    # unittest.main()

    loader = unittest.defaultTestLoader
    tsSuite = loader.loadTestsFromTestCase(TestTimeSeriesExtraction)
    profileSuite = loader.loadTestsFromTestCase(TestProfileExtraction)
    transectSuite = loader.loadTestsFromTestCase(TestTransectExtraction)
    trackSuite = loader.loadTestsFromTestCase(TestTrackExtraction)
    slabSuite = loader.loadTestsFromTestCase(TestSlabExtraction)
    tsHLSuite = loader.loadTestsFromTestCase(TestHLTimeSeriesExtractionXYZ)
    profileHLSuite = loader.loadTestsFromTestCase(TestHLProfileExtractionXYZ)
    transectHLSuite = loader.loadTestsFromTestCase(
        TestHLTransectExtractionForCoords)
    slabHLSuite = loader.loadTestsFromTestCase(TestHLSlabExtraction)
    # unittest.TextTestRunner().run(profileSuite)
    # unittest.TextTestRunner().run(transectSuite)
    # unittest.TextTestRunner().run(trackSuite)
    unittest.TextTestRunner().run(slabSuite)
    # unittest.TextTestRunner().run(tsHLSuite)
    # unittest.TextTestRunner().run(profileHLSuite)
    # unittest.TextTestRunner().run(transectHLSuite)
    # unittest.TextTestRunner().run(slabHLSuite)
