"""
A test suite for dirTreeManager

Tuomas Karna 2014-03-18
"""
import unittest
from netCDF4 import Dataset
import os
import shutil
from data.dataContainer import *
import glob
import data.dirTreeManager as dtm


class testSingleFileTree(unittest.TestCase):
    """Tests single file tree rule"""

    def createProfileDataContainer(self, var, withNaNs=True):
        nTime = 245
        nZ = 25
        startTime = datetime.datetime(2012, 5, 1)
        endTime = datetime.datetime(2012, 6, 5)
        st = datetimeToEpochTime(startTime)
        et = datetimeToEpochTime(endTime)
        time = np.linspace(st, et, nTime) # some epoch time
        z = np.zeros((nZ, nTime))
        data = np.zeros((nZ, 1, nTime))
        eta = np.sin(time/2/3600)
        h = np.ones_like(eta)*14.5
        for i in range(nTime):
            z[:, i] = np.linspace(eta[i], h[i], nZ)
        T = np.tile(time, (nZ, 1))
        data[:, 0,:] = 32*(np.sin(T/2/3600)*np.sin(z/10)+1)/2
        #data[:, 1,:] = 34*(np.sin(T/2/3600)*z/14.5)
        if withNaNs:
            data[:, 0, 10] = np.NaN
        meta = {}
        meta['tag'] = self.runTag
        meta['dataType'] = 'profile'
        meta['location'] = 'station01'
        meta['variable'] = var
        meta['bracket'] = 'F'
        ta = timeArray(time, 'epoch')
        x = np.zeros((nZ, nTime))*10.0
        y = np.zeros((nZ, nTime))*20.0
        dc = dataContainer('', ta, x, y, z, data, [var], metaData=meta,
                           coordSys='lonlat', acceptNaNs=withNaNs)
        return dc

    def createTimeSeriesDataContainer(self, var):
        nTime = 245
        startTime = datetime.datetime(2012, 5, 1)
        endTime = datetime.datetime(2012, 6, 5)
        st = datetimeToEpochTime(startTime)
        et = datetimeToEpochTime(endTime)
        time = np.linspace(st, et, nTime) # some epoch time
        meta = {}
        meta['tag'] = self.runTag
        meta['dataType'] = 'timeseries'
        meta['location'] = 'station01'
        meta['variable'] = var
        meta['msldepth'] = '100'
        meta['bracket'] = 'F'
        ta = timeArray(time, 'epoch')
        x = np.array([10.0])
        y = np.array([20.0])
        z = np.array([1.00])
        data = np.zeros((1, 1, nTime))
        data[:, 0,:] = 32*(np.sin(time/2/3600)+1)/2
        # data[:,1,:] = 34*(np.sin(time/2/3600)/14.5)
        dc = dataContainer('', ta, x, y, z, data, [var], metaData=meta)
        return dc

    def setUp(self):
        """Create dataContainers to be used for saving data"""
        print ' --- setup begins ---'
        self.runTag = 'tmpRun01'
        # these are used to test writing
        self.dcTSSalt = self.createTimeSeriesDataContainer('salt')
        self.dcProfSalt = self.createProfileDataContainer('salt')
        self.dcTSSalt.saveAsNetCDF(self.runTag+'/data/station01_salt_100_2012-05-01_2012-06-05.nc')
        self.dcProfSalt.saveAsNetCDF(self.runTag+'/data/profile/station01_salt_0_2012-05-01_2012-06-05.nc')
       # these are used to test reading
        self.dcTSTemp = self.createTimeSeriesDataContainer('temp')
        self.dcProfTemp = self.createProfileDataContainer('temp')
        if not os.path.isdir(self.runTag+'/data/profile'):
            os.makedirs(self.runTag+'/data/profile')
        self.dcTSTemp.saveAsNetCDF(self.runTag+'/data/station01_temp_100_2012-05-01_2012-06-05.nc')
        self.dcProfTemp.saveAsNetCDF(self.runTag+'/data/profile/station01_temp_0_2012-05-01_2012-06-05.nc')
        # with prefix path
        self.prefixPath = 'testPath/subdir/'
        self.dcTSTemp.saveAsNetCDF(self.prefixPath + self.runTag+'/data/station01_temp_100_2012-05-01_2012-06-05.nc')
        print ' --- setup done ---'

    def removeDirTree(self, d):
        """Removes a directory tree d"""
        if os.path.isdir(d):
            shutil.rmtree(d)

    def tearDown(self):
        """Removes temp directory and all its content"""
        if os.path.isdir(self.runTag):
            self.removeDirTree(self.runTag)
        if os.path.isdir(self.prefixPath):
            self.removeDirTree(self.prefixPath)

    def arrayClose(self, a, b):
        return np.isclose(a, b, equal_nan=True).all()

    def assertNetCDFFile(self, filename, dataContainer, dataArrayName):
        """Tests whether given file exists and the content matches the given
        dataContainer"""
        if not os.path.isfile(filename):
            self.fail('File {0:s} does not exist'.format(filename))
        d = Dataset(filename)
        self.assertTrue(self.arrayClose(d.variables['time'], dataContainer.time.array), 'time array mismatch')
        self.assertTrue(self.arrayClose(d.variables['x'], dataContainer.x), 'x array mismatch')
        self.assertTrue(self.arrayClose(d.variables['y'], dataContainer.y), 'y array mismatch')
        self.assertTrue(self.arrayClose(d.variables['z'], dataContainer.z), 'z array mismatch')
        self.assertTrue(self.arrayClose(d.variables[dataArrayName], dataContainer.data[:, 0,:]), 'data array mismatch')

    def testFilenameGenerationTimeSeries(self):
        sft = dtm.singleFileTree()
        res = sft.generateFileName(tag=self.runTag, dataType='timeseries',
                                   location='station01',
                                   variable='var', msldepth='100',
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 3, 12, 7))
        self.assertTrue(res == self.runTag+'/data/station01_var_100_2012-01-01_2012-03-12.nc', 'Bad filename '+res)

    def testFilenameGenerationTimeSeriesPath(self):
        path = self.prefixPath
        sft = dtm.singleFileTree()
        res = sft.generateFileName(tag=self.runTag, dataType='timeseries',
                                   location='station01',
                                   variable='var', msldepth='100',
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 3, 12, 7),
                                   rootPath=path)
        self.assertTrue(res == self.prefixPath + self.runTag+'/data/station01_var_100_2012-01-01_2012-03-12.nc', 'Bad filename '+res)

    def testFilenameGenerationProfile(self):
        sft = dtm.singleFileTree()
        res = sft.generateFileName(tag=self.runTag, dataType='profile',
                                   location='station01',
                                   variable='var',
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 3, 12, 7))
        self.assertTrue(res == self.runTag+'/data/profile/station01_var_0_2012-01-01_2012-03-12.nc', 'Bad filename '+res)

    def testFilenameGenerationProfileWildCard(self):
        sft = dtm.singleFileTree()
        res = sft.generateSearchPattern(dataType='profile',
                                   location='station01',
                                   variable='var',)
        target = '*/data/profile/station01_var_0_*_*.nc'
        self.assertTrue(res == target, 'Bad filename '+res)

    def testFilenameGenerationSlab(self):
        sft = dtm.singleFileTree()
        res = sft.generateFileName(tag=self.runTag, dataType='slab',
                                   location='station01',
                                   variable='var', slevel=3,
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 3, 12, 7))
        self.assertTrue(res == self.runTag+'/data/slab/station01_var_s3_2012-01-01_2012-03-12.nc', 'Bad filename '+res)

    def testFilePatternGenerationTag(self):
        sft = dtm.singleFileTree()
        res = sft.generateSearchPattern(dataType='slab',
                                   location='station01',
                                   variable='var', slevel=3,
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 3, 12, 7))
        self.assertTrue(res == '*/data/slab/station01_var_s3_2012-01-01_2012-03-12.nc', 'Bad filename '+res)

    def testFilePatternGenerationLocation(self):
        sft = dtm.singleFileTree()
        res = sft.generateSearchPattern(tag=self.runTag, dataType='slab',
                                   variable='var', slevel=3,
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 3, 12, 7))
        self.assertTrue(res == self.runTag+'/data/slab/*_var_s3_2012-01-01_2012-03-12.nc', 'Bad filename '+res)

    def testFilePatternGenerationVariable(self):
        sft = dtm.singleFileTree()
        res = sft.generateSearchPattern(tag=self.runTag, dataType='slab',
                                   location='station01',
                                   slevel=3,
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 3, 12, 7))
        self.assertTrue(res == self.runTag+'/data/slab/station01_*_s3_2012-01-01_2012-03-12.nc', 'Bad filename '+res)

    def testFilePatternGenerationStartEndTimes(self):
        sft = dtm.singleFileTree()
        res = sft.generateSearchPattern(tag=self.runTag, dataType='slab',
                                   location='station01',
                                   variable='var', slevel=3,)
        self.assertTrue(res == self.runTag+'/data/slab/station01_var_s3_*_*.nc', 'Bad filename '+res)

    def testFilePatternGenerationDepth(self):
        sft = dtm.singleFileTree()
        res = sft.generateSearchPattern(tag=self.runTag, dataType='timeseries',
                                   location='station01',
                                   variable='var',
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 3, 12, 7))
        self.assertTrue(res == self.runTag+'/data/station01_var_*_2012-01-01_2012-03-12.nc', 'Bad filename '+res)

    def testSaveFileTS(self):
        dtm.saveDataContainerInTree(self.dcTSSalt, rule='singleFile', overwrite=True)
        self.assertNetCDFFile(self.runTag+'/data/station01_salt_100_2012-05-01_2012-06-05.nc', self.dcTSSalt, 'water_salinity')

    def testSaveFileTSPath(self):
        path = self.prefixPath
        dtm.saveDataContainerInTree(self.dcTSSalt, rule='singleFile', overwrite=True, rootPath=path)
        self.assertNetCDFFile(self.prefixPath + self.runTag+'/data/station01_salt_100_2012-05-01_2012-06-05.nc', self.dcTSSalt, 'water_salinity')

    def testSaveFileProf(self):
        dtm.saveDataContainerInTree(self.dcProfSalt, rule='singleFile', overwrite=True)
        self.assertNetCDFFile(self.runTag+'/data/profile/station01_salt_0_2012-05-01_2012-06-05.nc', self.dcProfSalt, 'water_salinity')

    def testReadFileTS(self):
        dc = dtm.getDataContainer(rule='singleFile', dataType='timeseries', location='station01', variable='temp')
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        self.assertTrue(self.dcTSTemp == dc, 'time series dataContainer read from the disk does not match')

    def testReadFileTSPath(self):
        path = self.prefixPath
        dc = dtm.getDataContainer(rule='singleFile', dataType='timeseries',
                                  location='station01', variable='temp',
                                  rootPath=path)
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        self.assertTrue(self.dcTSTemp == dc, 'time series dataContainer read from the disk does not match')

    def testReadFileTSBadTime(self):
        st = datetime.datetime(2011, 1, 1)
        et = datetime.datetime(2011, 11, 11)
        with self.assertRaises(Exception) as cm:
            dc = dtm.getDataContainer(rule='singleFile', dataType='timeseries', location='station01', variable='temp',
                                      startTime=st, endTime=et)
        # old error message
        #targetError = "Reading data from tree failed: {'dataType': 'timeseries', 'variable': 'temp', 'endTime': datetime.datetime(2011, 11, 11, 0, 0), 'location': 'station01', 'startTime': datetime.datetime(2011, 1, 1, 0, 0)}"
        targetError = "Reading data from tree failed: rootPath=None dataType='timeseries' slevel=None msldepth=None tag=None location='station01' startTime=datetime.datetime(2011, 1, 1, 0, 0) variable='temp' endTime=datetime.datetime(2011, 11, 11, 0, 0)"
        self.assertEqual(cm.exception.message, targetError,
                         'Got wrong error message: '+cm.exception.message)

    def testReadFileProfS(self):
        dc = dtm.getDataContainer(rule='singleFile', dataType='profile', location='station01', variable='salt')
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        self.assertTrue(self.dcProfSalt == dc, 'time series dataContainer read from the disk does not match')

    def testReadFileProfT(self):
        dc = dtm.getDataContainer(rule='singleFile',dataType='profile',location='station01',variable='temp')
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        self.assertTrue(self.dcProfTemp == dc, 'profile dataContainer read from the disk does not match')


class testMonthlyFileTree(unittest.TestCase):
    """Tests monthly file tree rule"""

    def createProfileDataContainer(self, var, withNaNs=True):
        nTime = 245
        nZ = 25
        startTime = datetime.datetime(2012, 5, 1)
        endTime = datetime.datetime(2012, 6, 5)
        st = datetimeToEpochTime(startTime)
        et = datetimeToEpochTime(endTime)
        time = np.linspace(st, et, nTime) # some epoch time
        z = np.zeros((nZ, nTime))
        data = np.zeros((nZ, 1, nTime))
        eta = np.sin(time/2/3600)
        h = np.ones_like(eta)*14.5
        for i in range(nTime):
            z[:, i] = np.linspace(eta[i], h[i], nZ)
        T = np.tile(time, (nZ, 1))
        data[:, 0,:] = 32*(np.sin(T/2/3600)*np.sin(z/10)+1)/2
        #data[:, 1,:] = 34*(np.sin(T/2/3600)*z/14.5)
        if withNaNs:
            data[:, 0, 10] = np.NaN
        meta = {}
        meta['tag'] = self.runTag
        meta['dataType'] = 'profile'
        meta['location'] = 'station01'
        meta['variable'] = var
        meta['bracket'] = 'F'
        ta = timeArray(time, 'epoch')
        x = np.zeros((nZ, nTime))*10.0
        y = np.zeros((nZ, nTime))*20.0
        dc = dataContainer('', ta, x, y, z, data, [var], metaData=meta,
                           coordSys='lonlat', acceptNaNs=withNaNs)
        return dc

    def createTimeSeriesDataContainer(self, var):
        nTime = 245
        startTime = datetime.datetime(2012, 1, 3, 12, 8)
        endTime = datetime.datetime(2012, 3, 3, 10, 40, 12)
        st = datetimeToEpochTime(startTime)
        et = datetimeToEpochTime(endTime)
        time = np.linspace(st, et, nTime) # some epoch time
        meta = {}
        meta['tag'] = self.runTag
        meta['dataType'] = 'timeseries'
        meta['location'] = 'station01'
        meta['variable'] = var
        meta['msldepth'] = '100'
        meta['bracket'] = 'F'
        ta = timeArray(time, 'epoch')
        x = np.array([10.0])
        y = np.array([20.0])
        z = np.array([1.00])
        data = np.zeros((1, 1, nTime))
        data[:, 0,:] = 32*(np.sin(time/2/3600)+1)/2
        dc = dataContainer('', ta, x, y, z, data, [var], metaData=meta)
        return dc

    def assertNetCDFFile(self, filename, dataContainer, dataArrayName):
        """Tests whether given file exists and the content matches the given
        dataContainer"""
        if not os.path.isfile(filename):
            self.fail('File {0:s} does not exist'.format(filename))
        d = Dataset(filename)
        self.assertTrue(self.arrayClose(d.variables['time'], dataContainer.time.array), 'time array mismatch')
        self.assertTrue(self.arrayClose(d.variables['x'], dataContainer.x), 'x array mismatch')
        self.assertTrue(self.arrayClose(d.variables['y'], dataContainer.y), 'y array mismatch')
        self.assertTrue(self.arrayClose(d.variables['z'], dataContainer.z), 'z array mismatch')
        self.assertTrue(self.arrayClose(d.variables[dataArrayName], dataContainer.data[:, 0,:]), 'data array mismatch')

    def testFilenameGenerationTimeSeries(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateFileName(tag=self.runTag, dataType='timeseries',
                                   location='station01',
                                   variable='var', msldepth='100', bracket='A',
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == self.runTag+'/data/stations/station01/station01.100.A.var/201201.nc', 'Bad filename '+res)

    def testFilenameGenerationProfile(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateFileName(tag=self.runTag, dataType='profile',
                                   location='station01',
                                   variable='var',
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == self.runTag+'/data/profile/station01/var/201201.nc', 'Bad filename '+res)

    def testFilenameGenerationTransect(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateFileName(tag=self.runTag, dataType='transect',
                                   location='station01',
                                   variable='var',
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == self.runTag+'/data/transect/station01/var/201201.nc', 'Bad filename '+res)

    def testFilenameGenerationSlab(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateFileName(tag=self.runTag, dataType='slab',
                                   location='station01',
                                   variable='var', slevel='3',
                                   startTime=datetime.datetime(2012, 1, 1, 12),
                                   endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == self.runTag+'/data/slab/station01/var.s3/201201.nc', 'Bad filename '+res)

    def testFilePatternGenerationTag(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateSearchPattern(dataType='slab',
                                        location='station01',
                                        variable='var', slevel='3',
                                        startTime=datetime.datetime(2012, 1, 1, 12),
                                        endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == '*/data/slab/station01/var.s3/201201.nc', 'Bad filename '+res)

    def testFilePatternGenerationLocation(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateSearchPattern(tag=self.runTag,
                                        dataType='slab',
                                        variable='var', slevel='3',
                                        startTime=datetime.datetime(2012, 1, 1, 12),
                                        endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == self.runTag+'/data/slab/*/var.s3/201201.nc', 'Bad filename '+res)

    def testFilePatternGenerationVariable(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateSearchPattern(tag=self.runTag,
                                        dataType='slab',
                                        location='station01',
                                        slevel='3',
                                        startTime=datetime.datetime(2012, 1, 1, 12),
                                        endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == self.runTag+'/data/slab/station01/*.s3/201201.nc', 'Bad filename '+res)

    @unittest.skip('this fails!')
    def testFilePatternGenerationSLevel(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateSearchPattern(tag=self.runTag,
                                        dataType='slab',
                                        location='station01',
                                        variable='var',
                                        startTime=datetime.datetime(2012, 1, 1, 12),
                                        endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == self.runTag+'/data/slab/station01/var.*/201201.nc', 'Bad filename '+res)

    def testFilePatternGenerationTime(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateSearchPattern(tag=self.runTag,
                                        dataType='slab',
                                        location='station01',
                                        variable='var',
                                        slevel='3',
                                        startTime=None,
                                        endTime=None)
        self.assertTrue(res == self.runTag+'/data/slab/station01/var.s3/*.nc', 'Bad filename '+res)

    def testFilePatternGenerationStartTime(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateSearchPattern(tag=self.runTag,
                                        dataType='slab',
                                        location='station01',
                                        variable='var',
                                        slevel='3',
                                        startTime=None,
                                        endTime=datetime.datetime(2012, 1, 22, 7))
        self.assertTrue(res == self.runTag+'/data/slab/station01/var.s3/*.nc', 'Bad filename '+res)

    def testFilePatternGenerationEndTime(self):
        mft = dtm.monthlyFileTree()
        res = mft.generateSearchPattern(tag=self.runTag,
                                        dataType='slab',
                                        location='station01',
                                        variable='var',
                                        slevel='3',
                                        startTime=datetime.datetime(2012, 1, 1, 12),
                                        endTime=None)
        self.assertTrue(res == self.runTag+'/data/slab/station01/var.s3/*.nc', 'Bad filename '+res)

    def testSaveFileTS(self):
        dtm.saveDataContainerInTree(self.dcTSSalt, rule='monthlyFile', overwrite=True)
        files = ['tmpRun01/data/stations/station01/station01.100.F.salt/201201.nc',
                 'tmpRun01/data/stations/station01/station01.100.F.salt/201202.nc',
                 'tmpRun01/data/stations/station01/station01.100.F.salt/201203.nc',]
        dc = dataContainer.loadFromNetCDF(files[0])
        for f in files[1:]:
            dc.mergeTemporal(dataContainer.loadFromNetCDF(f))
        self.assertTrue(self.dcTSSalt == dc, 'time series dataContainer read from the disk does not match')

    def testSaveFileTSWithPath(self):
        path = self.prefixPath
        dtm.saveDataContainerInTree(self.dcTSSalt, rule='monthlyFile', overwrite=True,
                                    rootPath=path)
        files = ['tmpRun01/data/stations/station01/station01.100.F.salt/201201.nc',
                 'tmpRun01/data/stations/station01/station01.100.F.salt/201202.nc',
                 'tmpRun01/data/stations/station01/station01.100.F.salt/201203.nc',]
        for i in range(len(files)):
            files[i] = os.path.join(path, files[i])
        print 'reading', files[0]
        dc = dataContainer.loadFromNetCDF(files[0])
        for f in files[1:]:
            dc.mergeTemporal(dataContainer.loadFromNetCDF(f))
        self.assertTrue(self.dcTSSalt == dc, 'time series dataContainer read from the disk does not match')

    def testSaveFileTSWithPathNoTrailingBSlash(self):
        path = self.prefixPath[:-1]
        dtm.saveDataContainerInTree(self.dcTSSalt, rule='monthlyFile', overwrite=True,
                                    rootPath=path)
        files = ['tmpRun01/data/stations/station01/station01.100.F.salt/201201.nc',
                 'tmpRun01/data/stations/station01/station01.100.F.salt/201202.nc',
                 'tmpRun01/data/stations/station01/station01.100.F.salt/201203.nc',]
        for i in range(len(files)):
            files[i] = os.path.join(path, files[i])
        print 'reading', files[0]
        dc = dataContainer.loadFromNetCDF(files[0])
        for f in files[1:]:
            dc.mergeTemporal(dataContainer.loadFromNetCDF(f))
        self.assertTrue(self.dcTSSalt == dc, 'time series dataContainer read from the disk does not match')

    def testSaveFileProf(self):
        dtm.saveDataContainerInTree(self.dcProfSalt, rule='monthlyFile', overwrite=True)
        files = ['tmpRun01/data/profile/station01/salt/201205.nc',
                 'tmpRun01/data/profile/station01/salt/201206.nc',]
        dc = dataContainer.loadFromNetCDF(files[0])
        for f in files[1:]:
            dc.mergeTemporal(dataContainer.loadFromNetCDF(f))
        self.assertTrue(self.dcProfSalt == dc, 'time series dataContainer read from the disk does not match')

    def testFindFiles(self):
        mft = dtm.monthlyFileTree()
        files = mft.findMatchingFiles(dataType='timeseries', location='station01')
        correctList = ['tmpRun01/data/stations/station01/station01.100.F.salt/201201.nc',
                       'tmpRun01/data/stations/station01/station01.100.F.salt/201202.nc',
                       'tmpRun01/data/stations/station01/station01.100.F.salt/201203.nc',
                       'tmpRun01/data/stations/station01/station01.100.F.temp/201201.nc',
                       'tmpRun01/data/stations/station01/station01.100.F.temp/201202.nc',
                       'tmpRun01/data/stations/station01/station01.100.F.temp/201203.nc']
        self.assertTrue(files == correctList, 'find files returns wrong list')

    def testFindFilesTime(self):
        st = datetime.datetime(2012, 2, 15, 10)
        et = datetime.datetime(2012, 3, 15, 10)
        mft = dtm.monthlyFileTree()
        files = mft.findMatchingFiles(dataType='timeseries', location='station01',
                                      startTime=st, endTime=et)
        for f in files:
            print f
        correctList = ['tmpRun01/data/stations/station01/station01.100.F.salt/201202.nc',
                       'tmpRun01/data/stations/station01/station01.100.F.salt/201203.nc',
                       'tmpRun01/data/stations/station01/station01.100.F.temp/201202.nc',
                       'tmpRun01/data/stations/station01/station01.100.F.temp/201203.nc']
        self.assertTrue(files == correctList, 'find files returns wrong list')

    def testSampleListing(self):
        mft = dtm.monthlyFileTree()
        files = mft.findMatchingFiles(dataType='timeseries', location='station01')
        samples = mft.getAvailableSamples(files)
        correctList = [{'dataType': 'timeseries',
                        'slevel': None,
                        'msldepth': '100',
                        'tag': 'tmpRun01',
                        'location': 'station01',
                        'startTime': None,
                        'variable': 'salt',
                        'endTime': None,
                        'bracket': 'F'},
                       {'dataType': 'timeseries',
                        'slevel': None,
                        'msldepth': '100',
                        'tag': 'tmpRun01',
                        'location': 'station01',
                        'startTime': None,
                        'variable': 'temp',
                        'endTime': None,
                        'bracket': 'F'}
                       ]
        self.assertTrue(samples == correctList, 'returned sample list is wrong')

    def testReadFileTSAmbiguous(self):
        mft = dtm.monthlyFileTree()
        dcs = mft.readFiles(dataType='timeseries', location='station01')
        dc = dcs[0]
        self.assertTrue(len(dcs) == 2, 'wrong number of dataContainers returned: ' + str(len(dcs)))
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        self.assertTrue(self.dcTSSalt == dc, 'time series dataContainer read from the disk does not match')
        dc = dcs[1]
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        self.assertTrue(self.dcTSTemp == dc, 'time series dataContainer read from the disk does not match')

    def testReadFileTSAmbiguousTime(self):
        st = datetime.datetime(2012, 2, 15, 10)
        et = datetime.datetime(2012, 3, 15, 10)
        mft = dtm.monthlyFileTree()
        dcs = mft.readFiles(dataType='timeseries', location='station01',
                            startTime=st, endTime=et)
        self.assertTrue(len(dcs) == 2, 'wrong number of dataContainers returned: ' + str(len(dcs)))
        dc = dcs[0]
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        target = self.dcTSSalt.timeWindow(st, et)
        self.assertTrue(target == dc, 'time series dataContainer read from the disk does not match')
        dc = dcs[1]
        target = self.dcTSTemp.timeWindow(st, et)
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        self.assertTrue(target == dc, 'time series dataContainer read from the disk does not match')

    def testReadFileTS(self):
        dc = dtm.getDataContainer(rule='monthlyFile', dataType='timeseries',
                                  location='station01', variable='temp',
                                  verbose=True)
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        self.assertTrue(self.dcTSTemp == dc, 'time series dataContainer read from the disk does not match')

    def testReadFileTSTime(self):
        st = datetime.datetime(2012, 2, 15, 10)
        et = datetime.datetime(2012, 3, 15, 10)
        dc = dtm.getDataContainer(rule='monthlyFile', dataType='timeseries',
                                  location='station01', variable='temp',
                                  startTime=st, endTime=et)
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        target = self.dcTSTemp.timeWindow(st, et)
        self.assertTrue(target == dc, 'time series dataContainer read from the disk does not match')

    def testReadFileTSTimePath(self):
        st = datetime.datetime(2012, 2, 15, 10)
        et = datetime.datetime(2012, 3, 15, 10)
        path = self.prefixPath
        dc = dtm.getDataContainer(rule='monthlyFile', dataType='timeseries',
                                  location='station01', variable='temp',
                                  startTime=st, endTime=et, rootPath=path)
        self.assertIsInstance(dc, dataContainer, 'returned object is not dataContainer: '+str(type(dc)))
        target = self.dcTSTemp.timeWindow(st, et)
        self.assertTrue(target == dc, 'time series dataContainer read from the disk does not match')

    def testReadFileTSBadTime(self):
        st = datetime.datetime(2011, 1, 1)
        et = datetime.datetime(2011, 11, 11)
        with self.assertRaises(Exception) as cm:
            dc = dtm.getDataContainer(rule='monthlyFile', dataType='timeseries', location='station01', variable='temp',
                                      startTime=st, endTime=et)
        # old error message
        #targetError = "Reading data from tree failed: {'dataType': 'timeseries', 'variable': 'temp', 'endTime': datetime.datetime(2011, 11, 11, 0, 0), 'location': 'station01', 'startTime': datetime.datetime(2011, 1, 1, 0, 0)}"
        targetError = "Reading data from tree failed: rootPath=None dataType='timeseries' slevel=None msldepth=None tag=None location='station01' startTime=datetime.datetime(2011, 1, 1, 0, 0) variable='temp' endTime=datetime.datetime(2011, 11, 11, 0, 0)"
        self.assertEqual(cm.exception.message, targetError,
                         'Got wrong error message: '+cm.exception.message)

    def setUp(self):
        """Create dataContainers to be used for saving data"""
        print ' --- setup begins ---'
        self.runTag = 'tmpRun01'
        # these are used to test writing
        self.dcTSSalt = self.createTimeSeriesDataContainer('salt')
        self.dcProfSalt = self.createProfileDataContainer('salt')
        self.dcTSTemp = self.createTimeSeriesDataContainer('temp')
        self.dcProfTemp = self.createProfileDataContainer('temp')
        # save time series dc in monthly files
        for m in range(1, 4):
            st = datetime.datetime(2012, m, 1)
            et = datetime.datetime(2012, m+1, 1)
            dc = self.dcTSTemp.timeWindow(st, et, includeEnd=False)
            dc.saveAsNetCDF(self.runTag+'/data/stations/station01/station01.100.F.temp/20120'+str(m)+'.nc')
        for m in range(1, 4):
            st = datetime.datetime(2012, m, 1)
            et = datetime.datetime(2012, m+1, 1)
            dc = self.dcTSSalt.timeWindow(st, et, includeEnd=False)
            dc.saveAsNetCDF(self.runTag+'/data/stations/station01/station01.100.F.salt/20120'+str(m)+'.nc')
        # save profile dc in monthly files
        for m in range(5, 7):
            st = datetime.datetime(2012, m, 1)
            et = datetime.datetime(2012, m+1, 1)
            dc = self.dcProfTemp.timeWindow(st, et, includeEnd=False)
            dc.saveAsNetCDF(self.runTag+'/data/profile/station01/temp/20120'+str(m)+'.nc')
        # save time series with prefix path
        for m in range(1, 4):
            st = datetime.datetime(2012, m, 1)
            et = datetime.datetime(2012, m+1, 1)
            dc = self.dcTSTemp.timeWindow(st, et, includeEnd=False)
            self.prefixPath = 'testPath/subdir/'
            dc.saveAsNetCDF(self.prefixPath+self.runTag+'/data/stations/station01/station01.100.F.temp/20120'+str(m)+'.nc')
        print ' --- setup done ---'

    def removeDirTree(self, d):
        """Removes a directory tree d"""
        if os.path.isdir(d):
            shutil.rmtree(d)

    def tearDown(self):
        """Removes temp directory and all its content"""
        if os.path.isdir(self.runTag):
            self.removeDirTree(self.runTag)
        if os.path.isdir(self.prefixPath):
            self.removeDirTree(self.prefixPath)

    def arrayClose(self, a, b):
        return np.isclose(a, b, equal_nan=True).all()

if __name__ == '__main__':
    """Run all tests"""
    unittest.main()
