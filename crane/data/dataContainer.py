#!/usr/bin/python
"""
Generic implementation for a data container.
Can represent time series, 2D (transsect), vertical profile data.
Contains data, necessary information to interpret it, and some basic
data manipulation methods.

Tuomas Karna 2012-08-15
"""

import os
import numpy as np

from crane.data import timeArray
from crane.data import netcdfIO


class dataContainer(object):
    """
    Generic data container.
    Contains timeArray, spatial coordinates and data.
    """

    def __init__(self, description, time, x, y, z, data, fieldNames,
                 coordSys='whatever', metaData={}, acceptNaNs=False,
                 checkDataXDim=True, dtype=None):
        """Creates new dataContainer object

        time        -- timeArray object
        x, y, z     -- (scalar) or (nPoints,) or (nPoints,nTime) coordinate arrays
                       2-dimensional array if coodinates depend on time, scalar if
                       only one static point. All arrays must have the same number
                       of rows (nPoints). It is possible to set some arrays time
                       dependent (e.g. z) while other are static.
        data        -- (nPoints,nFields,nTime) array of p fields
        fieldNames  -- list of p strings
        coordSys    -- string representing the used coordinate system
        acceptNaNs  -- if False, raises an error if NaNs/Infs found in data array
        """
        if dtype is not None:
            self.dtype = dtype
        else:
            self.dtype = data.dtype
        self.acceptNaNs = acceptNaNs
        if len(time) == 0:
            raise Exception('time array is empty')
        if np.isnan(time.array).any() or np.isinf(time.array).any():
            raise Exception('time array contains nans/infs')
        if not acceptNaNs and (np.isnan(data).any() or np.isinf(data).any()):
            raise Exception('data array contains nans/infs')
        if len(data.shape) != 3:
            raise Exception('data must be of shape (nPoints, nFields, nTime)')
        if np.any(np.array(data.shape) == 0):
            print data.shape
            raise Exception('data is empty (has zero dimension)')

        if not hasattr(x, 'shape') or not x.shape:
            x = np.array([x])
        if not hasattr(y, 'shape') or not y.shape:
            y = np.array([y])
        if not hasattr(z, 'shape') or not z.shape:
            z = np.array([z])
        nX = x.shape[0]
        nY = y.shape[0]
        nZ = z.shape[0]
        if nY != nX or nZ != nX:
            print x.shape, y.shape, z.shape
            raise Exception('x, y, z must have same number of rows')
        self.xDependsOnTime = len(x.shape) == 2
        self.yDependsOnTime = len(y.shape) == 2
        self.zDependsOnTime = len(z.shape) == 2
        self.coordsDependOnTime = (self.xDependsOnTime or
                                   self.yDependsOnTime or
                                   self.zDependsOnTime)

        nTime = len(time)
        nPoints = nX
        nFields = data.shape[1]

        if self.xDependsOnTime and x.shape[1] != nTime:
            print x.shape, len(time)
            raise Exception('If x has 2nd dimension, it must match time')
        if self.yDependsOnTime and y.shape[1] != nTime:
            print y.shape, len(time)
            raise Exception('If y has 2nd dimension, it must match time')
        if self.zDependsOnTime and z.shape[1] != nTime:
            print z.shape, len(time)
            raise Exception('If z has 2nd dimension, it must match time')

        if checkDataXDim and data.shape[0] != nPoints:
            print data.shape, x.shape
            raise Exception('number of rows in data must match x,y,z')
        if data.shape[2] != nTime:
            print data.shape, time.array.shape
            raise Exception('3rd dimension of data array must match time')
        if not isinstance(fieldNames, list):
            raise Exception('fieldNames must be a list of strings')
        if len(fieldNames) != nFields:
            print len(fieldNames), data.shape
            raise Exception('length of fieldNames must match cols in data')

        self.time = time
        self.x = x if dtype is None else x.astype(self.dtype)
        self.y = y if dtype is None else y.astype(self.dtype)
        self.z = z if dtype is None else z.astype(self.dtype)
        self.data = data if dtype is None else data.astype(self.dtype)
        self.description = str(description)
        self.coordSys = str(coordSys)
        self.fieldNames = fieldNames
        self.metaData = dict(metaData)

    def setMetaData(self, name, value=None):
        """Assigns metaData for given name."""
        if isinstance(name, dict):
            for n in name:
                self.setMetaData(n, name[n])
        else:
            self.metaData[name] = value

    def getMetaData(self, name=None, suppressError=False):
        """Returns metaData corresponding to the given name. Raises error if metadata is not found. If suppressError=True, returns None instead of raising an error."""
        if name is None:
            return dict(self.metaData)
        if suppressError:
            return self.metaData.get(name, None)
        return self.metaData[name]

    def hasMetaData(self, name):
        """Returns True if metaData for given name is present."""
        return name in self.metaData.keys()

    @classmethod
    def fromTimeSeries(cls, description, time, data, fieldNames,
                       x=np.nan, y=np.nan, z=np.nan,
                       timeFormat='', coordSys='', metaData=None):
        """Alternative constructor.
        Creates new dataContainer object from time series.

        x, y, z    -- coordinates, scalars of arrays of shape (1,)
        data       -- time series data, shape (nTime, )
                      or (nTime,nFields) or (nFields, nTime)
        time       -- singleton time array or timeArray object
        timeFormat -- time format string if time is an array
        """
        if isinstance(time, np.ndarray):
            if timeFormat == '':
                raise Exception(
                    'If time is an array, timeFormat must be specified')
            time = timeArray.timeArray(time, timeFormat)

        nTime = len(time)
        nPoints = 1
        nFields = 1
        # TODO a better way for doing this?
        if len(data.shape) > 1:
            # assume that time is along rows
            if data.shape[0] == nTime:
                nFields = data.shape[1]
            elif data.shape[1] == nTime:
                nFields = data.shape[0]
            else:
                print time.array.shape, data.shape
                raise Exception(
                    'data array size does not match time series specs')
        elif data.shape[0] != nTime:
            print len(time), data.shape
            raise Exception('data array length does not match time series')

        data = data.reshape((nPoints, nFields, nTime))  # correct shape
        return cls(
            description,
            time,
            x,
            y,
            z,
            data,
            fieldNames,
            coordSys,
            metaData)

    def copy(self):
        """Deep copy, all numpy arrays are copied instead of referenced."""
        return dataContainer(
            self.description,
            self.time.copy(),
            self.x.copy(),
            self.y.copy(),
            self.z.copy(),
            self.data.copy(),
            list(
                self.fieldNames),
            self.coordSys,
            metaData=dict(
                self.metaData),
            acceptNaNs=self.acceptNaNs,
            dtype=self.dtype)

    def __eq__(self, other):
        """True if all data is equal to other."""
        if not isinstance(other, type(self)):
            return False
        if self.time != other.time:
            return False
        if not np.allclose(self.x, other.x):
            return False
        if not np.allclose(self.y, other.y):
            return False
        if not np.allclose(self.z, other.z):
            return False
        try:
            # raises error if arrays differ, treats NaNs as numbers
            np.isclose(self.data, other.data, equal_nan=True).all()
            #np.testing.assert_array_equal( self.data, other.data )
        except AssertionError as e:
            return False
        if self.description != other.description:
            return False
        if self.fieldNames != other.fieldNames:
            return False
        return True

    def __ne__(self, other):
        """True if some data is equal to other."""
        return not self.__eq__(other)

    def __str__(self):
        """
        Dumps a string summarizing the data content
        """
        outstr = ' * ' + self.__class__.__name__ + '\n'
        outstr += 'coordSys: ' + self.coordSys + '\n'
        for n in self.metaData:
            outstr += n + ': ' + str(self.metaData[n]) + '\n'
        outstr += str(self.time) + '\n'
        outstr += '{0:d} points'.format(self.x.shape[0]) + '\n'
        outstr += arrayRangeStr('x', self.x) + '\n'
        outstr += arrayRangeStr('y', self.y) + '\n'
        outstr += arrayRangeStr('z', self.z) + '\n'
        outstr += '{0:d} fields'.format(self.data.shape[1])
        for i in range(self.data.shape[1]):
            outstr += '\n'
            outstr += arrayRangeStr(self.fieldNames[i], self.data[:, i, :])
        return outstr

    def xyzMatch(self, other):
        """True if x,y,z arrays of the two containers are the same"""
        if (np.isnan(self.x).all() and np.isnan(other.x).all() and
                np.isnan(self.y).all() and np.isnan(other.y).all() and
                np.isnan(self.z).all() and np.isnan(other.z).all()):
            return (self.x.shape == other.x.shape and
                    self.y.shape == other.y.shape and
                    self.z.shape == other.z.shape)
        else:
            return (
                np.array_equal(
                    np.nan_to_num(
                        self.x),
                    np.nan_to_num(
                        other.x)) and np.array_equal(
                    np.nan_to_num(
                        self.y),
                    np.nan_to_num(
                        other.y)) and np.array_equal(
                    np.nan_to_num(
                        self.z),
                    np.nan_to_num(
                        other.z)))

    def getFieldArray(self, fieldName):
        """
        Returns the requested field in an array.

        Parameters
        ----------
        fieldName : str or int
            fieldName to extract. If int, the index of the field to extract.
        """
        if isinstance(fieldName, str):
            if fieldName in self.fieldNames:
                return self.data[:, self.fieldNames.index(fieldName), :]
            else:
                raise Exception('given field not found', fieldName)
        else:
            # assume index
            return self.data[:, fieldName, :]

    def extractFields(self, *fields, **kwargs):
        """
        Returns a dataContainer containing only the requested field.

        Parameters
        ----------
        fields : str or int
            fieldName to extract. If int, the index of the field to extract.
        copy   : bool
            if True copies the data array instead of using a view
        """
        copy = kwargs.get(copy, False)
        indices = []
        names = []
        for f in fields:
            if isinstance(f, str) or isinstance(f, unicode):
                # deduce index
                i = self.fieldNames.index(f)
            else:
                i = f
            indices.append(i)
            names.append(self.fieldNames[i])
        if copy:
            data = self.data[:, indices, :].copy()
        else:
            data = self.data[:, indices, :]

        return dataContainer(
            self.description,
            self.time,
            self.x,
            self.y,
            self.z,
            data,
            names,
            self.coordSys,
            self.metaData,
            acceptNaNs=True)

    def mergeTemporal(self, other, testSanity=True, acceptDuplicates=False):
        """
        Appends other data in the end of this container's time dimension.
        Only time and data fields are changed.
        """
        # TODO force increasing time (sort) ?
        if testSanity and self.fieldNames != other.fieldNames:
            print self.fieldNames, other.fieldNames
            raise Exception('Cannot merge, fieldNames do not agree')

        if self.xDependsOnTime:
            if testSanity and not other.xDependsOnTime:
                print self.xDependsOnTime, other.xDependsOnTime
                raise Exception('Cannot merge, x time dimension incompatible')
            self.x = np.concatenate((self.x, other.x), axis=1)
        elif testSanity and not np.array_equal(self.x, other.x):
            print self.x.shape, other.x.shape
            raise Exception('Cannot merge, x arrays do not match ')

        if self.yDependsOnTime:
            if testSanity and not other.yDependsOnTime:
                print self.yDependsOnTime, other.yDependsOnTime
                raise Exception('Cannot merge, y time dimension incompatible')
            self.y = np.concatenate((self.y, other.y), axis=1)
        elif testSanity and not np.array_equal(self.y, other.y):
            print self.y.shape, other.y.shape
            raise Exception('Cannot merge, y arrays do not match ')

        if self.zDependsOnTime:
            if testSanity and not other.zDependsOnTime:
                print self.zDependsOnTime, other.zDependsOnTime
                raise Exception('Cannot merge, z time dimension incompatible')
            self.z = np.concatenate((self.z, other.z), axis=1)
        elif testSanity and not np.array_equal(self.z, other.z):
            print self.z.shape, other.z.shape
            raise Exception('Cannot merge, z arrays do not match ')

        self.time.merge(other.time, acceptDuplicates=acceptDuplicates)
        self.data = np.concatenate((self.data, other.data), axis=2)

    def mergeFields(self, other):
        """
        Appends fields from other to this container
        """
        if not self.xyzMatch(other):
            raise Exception('xyz coordinates do not agree')
        if self.time != other.time:
            raise Exception('time arrays do not agree')

        self.fieldNames += other.fieldNames
        self.data = np.hstack((self.data, other.data))

    def changeTimeFormat(self, timeFormat, startDate=None):
        """
        Changes time stamps to given format. Data is not changed.
        """
        self.time = self.time.toFormat(timeFormat, startDate)

    def interpolateInTime(self, newTime, acceptNaNs=False):
        """
        Interpolates data to newTime. newTime is a timeArray object.
        Returns a new dataContainer object.
        Note: x,y,z are referenced rather than copied!
        """

        # test if extrapolation is needed
        if (newTime.asEpoch().max() > self.time.asEpoch().max() or
                newTime.asEpoch().min() < self.time.asEpoch().min()):
            print 'self.time:', self.time
            print 'newtime:', newTime
            print (newTime.asEpoch().min() - self.time.asEpoch().min(),
                   self.time.asEpoch().max() - newTime.asEpoch().max())
            raise Exception('Data does not cover reguested time period.')
        # reshape data to compatible time series form (nPoints*nFields,nTime)
        nPoints = self.x.shape[0]
        nFields = self.data.shape[1]
        data2 = self.data.reshape((nPoints * nFields, len(self.time)))
        # interpolate
        from scipy.interpolate import interp1d
        newdata = interp1d(self.time.toFormat(newTime), data2)(newTime)
        #from scipy.interpolate import UnivariateSpline
        # k is spline degree, s = smoothing condition for optimizing knots
        # s = 0 -> do not optimize, use data locations as knots
        #newdata =  UnivariateSpline( self.time.toFormat(newTime), np.copy(data2), k, s=0 )(newTime)
        newdata = newdata.reshape((nPoints, nFields, len(newTime)))

        newx = self.x
        newy = self.y
        newz = self.z
        # interpolate coordinates if time-dependent
        if self.xDependsOnTime:
            newx = interp1d(self.time.toFormat(newTime), self.x)(newTime)
        if self.yDependsOnTime:
            newy = interp1d(self.time.toFormat(newTime), self.y)(newTime)
        if self.zDependsOnTime:
            newz = interp1d(self.time.toFormat(newTime), self.z)(newTime)

        return dataContainer(
            self.description, newTime, newx, newy, newz, newdata, list(
                self.fieldNames), self.coordSys, dict(
                self.metaData), acceptNaNs=acceptNaNs)

    # TODO call timeArray.detectGaps
    def detectGaps(self, dt=None, gapFactor=5):
        """Detects gaps in the time series.

        Args:
        dt        -- (float) data sampling period. If None, taken as a mean step between data points.
        gapFactor -- (float) A factor to determine a minimum gap: gapFactor*dt

        Returns:
        gaps      -- (array) Indices that mark the beginning of each gap. (ngaps,)
        ranges    -- (array) Start and end indices of each contiguous data block. (ngaps+1,2)
        t         -- (array) time array in epoch format. (ntime,)
        """
        t = self.time.asEpoch().array
        steps = np.diff(t)
        if not dt:
            dt = np.mean(steps)
        gaps = np.nonzero(steps >= gapFactor * dt)[0]
        # convert to data ranges
        ranges = np.zeros((len(gaps) + 1, 2), dtype=int)
        ranges[1:, 0] = gaps + 1  # start indices
        ranges[:-1, 1] = gaps  # end indices
        ranges[-1, 1] = len(self.time) - 1
        return gaps, ranges, t

    def subsample(
            self,
            timeStamps=None,
            skipFactor=None,
            targetDt=None,
            currentDt=None,
            gapFactor=5):
        """Subsamples the data with the given skipFactor or targetDt.
        If timeStamps is not given, appropriate time indices will be estimated.
        Returns a new dataContainer.

        Args:
        timeStamps -- (ndarray) array of time indices to include
        skipFactor -- (int) subsampling factor, accept every skipFactor data point.
        targetDt -- (float) alternatively estimate skipFactor based on targetDt (sec).
        currentDt -- (float) data sampling period. If None, taken as a mean step between data points. Passed to detectGaps.
        gapFactor -- (float) A factor to determine a minimum gap: gapFactor*dt. Passed to detectGaps.

        Returns:
        newDC      -- (array) subsampled version of this dataContainer
        """
        if timeStamps is None:
            if skipFactor is None and targetDt is None:
                raise Exception(
                    'Either skipFactor or targetDt is required',
                    skipFactor,
                    targetDt)
            # detect gaps
            gaps, ranges, t = self.detectGaps(currentDt, gapFactor)
            if currentDt:
                dt = currentDt
            else:
                # estimate dt based on each data range
                dt_list = []
                for r in range(ranges.shape[0]):
                    dt = np.mean(np.diff(self.time[ranges[r, 0]:ranges[r, 1]]))
                    if np.isfinite(dt):
                        dt_list.append(dt)
                dt = np.array(dt_list).mean()
            if dt == 0 or ~np.isfinite(dt):
                print ranges[0, 0], ranges[0, 1]
                print self.time[ranges[0, 0]:ranges[0, 1]]
                raise Exception('weird dt:' + str(dt))
            if skipFactor is None:
                skipFactor = int(round(targetDt / dt))
            if skipFactor <= 1:
                print dt, targetDt, skipFactor
                print 'Warning: Requested subsampling dt smaller or equal to the original, skipping'
                return self
            print 'subsampling: dt', dt, targetDt, skipFactor
            # determine subsampling time stamps
            timeStamps = []
            for r in range(ranges.shape[0]):
                timeStamps.append(
                    np.arange(
                        ranges[
                            r, 0], ranges[
                            r, 1] + 1, skipFactor, dtype=int))
            timeStamps = np.concatenate(tuple(timeStamps), axis=0)
        # create new dataContainer
        newt = timeArray.timeArray(
            self.time[timeStamps],
            'epoch', acceptDuplicates=True)
        newData = self.data[:, :, timeStamps]
        x = self.x[:, timeStamps] if self.xDependsOnTime else self.x
        y = self.y[:, timeStamps] if self.yDependsOnTime else self.y
        z = self.z[:, timeStamps] if self.zDependsOnTime else self.z
        return dataContainer(
            self.description,
            newt,
            x,
            y,
            z,
            newData,
            self.fieldNames,
            self.coordSys,
            self.metaData,
            acceptNaNs=True,
            checkDataXDim=False)

    def alignTimes(self, other):
        """Aligns time stamps with the other dataContainer. In both dataContainers
        gaps are detected first and then mutually overlapping data ranges are
        computed. Data in other is restricted on the overlapping ranges. Data in
        self is additionally interpolated to the time steps of other. Returns
        aligned (self,other) dataContainers.
        Raises an exception if overlapping time intervals are not found.

        Args:
        other : dataContainer
            other container to interpolate to

        Returns:
        alignedSelf : dataContainer
            aligned version of self, interpolated to other's time steps
        alignedOther : dataContainer
            aligned version of other
        """
        # TODO get sampling rate from db and save in dataContainer
        selfDt = otherDt = None
        if (self.getMetaData('location', suppressError=True) == 'saturn03' and
                self.getMetaData('tag', suppressError=True) == 'obs'):
            selfDt = 360
        if (other.getMetaData('location', suppressError=True) == 'saturn03' and
                other.getMetaData('tag', suppressError=True) == 'obs'):
            otherDt = 360
        oTimeStamps = self.time.getAlignedTimeIndices(
            other.time, selfDt, otherDt)
        newTA = other.time.copy()
        newTA.array = newTA.array[oTimeStamps]
        if len(newTA) < 3:
            raise Exception('Aligned time series too short.')
        newSelf = self.interpolateInTime(newTA)
        x = other.x[:, oTimeStamps] if other.xDependsOnTime else other.x
        y = other.y[:, oTimeStamps] if other.yDependsOnTime else other.y
        z = other.z[:, oTimeStamps] if other.zDependsOnTime else other.z
        # TODO what to do with metadata ?
        newOther = dataContainer(
            other.description, newTA, x, y, z, other.data
            [:, :, oTimeStamps],
            list(other.fieldNames),
            other.coordSys, other.metaData)
        return newSelf, newOther

    def computeError(self, other):
        """Computes error versus other data container.
        Data in this container is interpolated on the time steps of the other dataContainer.
        Gaps are automatically detected and handled.
        Error is defined as other-self. Typically call err = obs.computeError( model ).
        Raises an exception if overlapping time intervals are not found.
        """
        if self.time != other.time:
            r, o = self.alignTimes(other)
        else:
            r = self  # reference
            o = other  # signal to compare
        errArray = o.data - r.data
        # TODO what kind of metadata?
        meta = o.metaData
        meta['errorSignal'] = True
        return dataContainer('error_' + o.description, o.time, o.x, o.y, o.z,
                             errArray, list(o.fieldNames), o.coordSys, meta)

    def timeWindow(self, startDate, endDate, includeEnd=False):
        """Returns a dataContainer, restricted to the given interval.
        startDate is included in the interval, endDate is not.
        startDate and endDate are datetime objects. Data in the returned array is refenreced, not copied.
        To obtain a copy use timeWindow(startDate, endDate).copy()
        """
        goodIx = self.time.getRangeIndices(startDate, endDate, includeEnd)
        if len(goodIx) == 0:
            raise Exception('Given time window out of range')
        t = timeArray.timeArray(
            self.time.array[goodIx],
            self.time.timeFormat,
            acceptDuplicates=True)
        d = self.data[:, :, goodIx]
        x = self.x[:, goodIx] if self.xDependsOnTime else self.x
        y = self.y[:, goodIx] if self.yDependsOnTime else self.y
        z = self.z[:, goodIx] if self.zDependsOnTime else self.z
        return dataContainer(
            self.description,
            t,
            x,
            y,
            z,
            d,
            self.fieldNames,
            self.coordSys,
            acceptNaNs=True,
            metaData=dict(
                self.metaData),
            checkDataXDim=False)

    def saveAsNetCDF(self, filename, dtype=None, overwrite=True,
                     compress=False, digits=None):
        """Saves data in netCDF format.
        dtype sets data type to any numpy format (except float16)
        overwrite=True forces writing over existing file, otherwise error is raised
        compress=True sets lossless zlib compression for netcdf4 dataset
        digits rounds data to given decimals before storing
        """
        if dtype is None:
            dtype = self.dtype

        nc = netcdfIO.netcdfIO(filename)
        nc.saveDataContainer(self, dtype, overwrite, compress, digits)

    @classmethod
    def loadFromNetCDF(
            cls,
            filename,
            startTime=None,
            endTime=None,
            includeEnd=False,
            verbose=True):
        """Creates a new dataContainer from netCDF file.
        """
        nc = netcdfIO.netcdfIO(filename)
        dc = nc.readToDataContainer(
            startTime,
            endTime,
            includeEnd=includeEnd,
            verbose=verbose)
        return dc

    def isTrack(self):
        """Returns true if the contents represent a time dependent track in space (x,y,z,t)."""
        if self.xDependsOnTime or self.yDependsOnTime or self.zDependsOnTime:
            result = True
            if self.xDependsOnTime:
                result *= self.x.shape[0] == 1 and self.x.shape[1] > 1
            if self.yDependsOnTime:
                result *= self.y.shape[0] == 1 and self.y.shape[1] > 1
            if self.zDependsOnTime:
                result *= self.z.shape[0] == 1 and self.z.shape[1] > 1
            return result
        return False

# ------ helper functions ------


def arrayRangeStr(fieldName, array):
    """creates a string that summarizes data ranges"""
    outstr = ''
    uarray = np.unique(array.flatten())
    N = len(uarray)
    if N == 1:
        tmpStr = '{name:s}: {shape:s} value: {v:g}'
        outstr = tmpStr.format(
            name=fieldName,
            shape=array.shape,
            v=float(
                uarray[0]))
    else:
        uarray = uarray[np.logical_not(np.isnan(uarray))]
        tmpStr = '{name:s}: {shape:s} range: {m:g} ... {M:g}'
        outstr = tmpStr.format(
            name=fieldName, shape=array.shape, m=float(
                np.min(uarray)), M=float(
                np.max(uarray)))
    return outstr
