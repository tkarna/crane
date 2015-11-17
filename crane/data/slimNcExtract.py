import os
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from glob import glob
import traceback

#from data.extractStation import fieldNameList
from crane.data import timeArray
from crane.data.dataContainer import *
from crane.data import gridUtils
from crane.data import meshContainer
from crane.data.selfeGridUtils import *
from crane.files import buildPoints
from crane.physicalVariableDefs import fieldNameList
import datetime
from crane.files import csvStationFile

# use consistent field names throughout the skill assessment package

VARS2D = ['elev', 'dahv']


def getNCVariableName(varStr):
    nc_names = {'elev': 'eta',
                'salt': 'S',
                'temp': 'T',
                }
    return nc_names.get(varStr, varStr)


class slimExtractBase(object):
    """Base class for all netcdf extract objects."""
    # TODO morf into generic base class, leave all model dependent stuff
    # undefined

    def __init__(self, path, var=None, verbose=False):
        """Intializes reader object."""
        fieldNameToFilename = {
            'temp': 'T',
            'elev': 'eta',
            'salt': 'S',
            'kine': '',
            'vdff': '',
            'tdff': '',
            'mixl': '',
            'hvel': '',
            'vert': '',
            'dens': '',
            'trcr_1': '',
            'trcr_2': '',
            'trcr_3': '',
            'trcr_4': '',
            'turbidity': ''}  # TODO belongs to derived class
        self.path = path
        self.component = 0
        self.fileTypeStr = fieldNameToFilename[var]
        self.headerIsRead = False
        self.verbose = verbose

    def generateFileName(self, iStack=None, fileTypeStr=None):
        """Returns full path to the netcdf file for iStack.
        If iStack==None, returns a pattern with '*' as a wildcard."""
        # TODO raise NotImplementedError('This method must be defined in the
        # derived class')
        if fileTypeStr is None:
            fileTypeStr = self.fileTypeStr
        stackStr = '*' if iStack is None else '{0:05d}'.format(iStack)
        fname = '{typeStr:s}_{stack:s}_COMP_{comp:d}.nc'.format(
            typeStr=fileTypeStr, stack=stackStr,
            comp=self.component)
        return os.path.join(self.path, fname)

    def getNCFile(self, iStack=None, fileTypeStr=None):
        """Opens netcdf file corresponding to the given stack number.
        If no stack number is given opens first matching file."""
        if fileTypeStr is None:
            fileTypeStr = self.fileTypeStr
        f = self.generateFileName(iStack, fileTypeStr)
        if iStack is None:
            # try to find a file that matches file name pattern
            pattern = f
            files = sorted(glob(pattern))
            if len(files) == 0:
                raise Exception('no files found in ' + pattern)
            f = files[0]
        if not os.path.isfile(f):
            raise IOError('File not found: ' + f)
        else:
            if self.verbose:
                print 'Opening file', f
            return NetCDFFile(f, 'r')

    def readHeader(self, ncfile=None):
        """
        Reads header of the netcdf file and prepares data structures.

        If ncfile is given, will read its header. Otherwise will search for first
        matching netcdf file in the path.
        """
        # TODO raise NotImplementedError('This method must be defined in the
        # derived class')
        ncfileGiven = ncfile is not None
        if self.verbose:
            print 'Reading header'
        if not ncfileGiven:
            ncfile = self.getNCFile()

        # read
        faceNOffset = 0  # ncfile.variables['face_nodes'].start_index
        self.faceNodes = ncfile.variables[
            'face_nodes'][:].astype(int) - faceNOffset
        self.nodeX = ncfile.variables['node_x'][:]
        self.nodeY = ncfile.variables['node_y'][:]

        self.nNodes = len(ncfile.dimensions['node'])
        self.nFaces = len(ncfile.dimensions['face'])
        self.nElemNodes = len(ncfile.dimensions['nFaceNodes'])  # ==3 always
        self.nTime = len(ncfile.dimensions['time'])
        if 'layers' in ncfile.dimensions:
            self.nVert = len(ncfile.dimensions['layers'])
        else:
            self.nVert = 0
        if self.verbose:
            print 'nodes', self.nNodes
            print 'elems', self.nFaces
            print 'verts', self.nVert
            print 'elem nodes', self.nElemNodes
            print 'time stamps', self.nTime
        timeStr = ' '.join(ncfile.variables['time'].base_date.split()[2:4])
        self.simulationStartTime = datetime.datetime.strptime(
            timeStr, '%Y-%m-%d %H:%M:%S')

        if 'node_lon' in ncfile.variables:
            self.node_lon = ncfile.variables['node_lon'][:]
        if 'node_lat' in ncfile.variables:
            self.node_lat = ncfile.variables['node_lat'][:]

        if 'edge_nodes' in ncfile.variables:
            self.edgeNodes = ncfile.variables['edge_nodes'][:]
        if 'edge_x' in ncfile.variables:
            self.edge_x = ncfile.variables['edge_x'][:]
            self.edge_y = ncfile.variables['edge_y'][:]
        if 'edge_lon' in ncfile.variables:
            self.edge_lon = ncfile.variables['edge_lon'][:]
            self.edge_lat = ncfile.variables['edge_lat'][:]

        # construct mesh search object
        self.meshSearch2d = gridUtils.meshSearch2d(
            self.nodeX, self.nodeY, self.faceNodes)
        self.headerIsRead = True
        if not ncfileGiven:
            ncfile.close()

    def getTime(self, ncfile):
        """Returns time stamps from given netCDF file in epoch format."""
        # TODO raise NotImplementedError('This method must be defined in the
        # derived class')
        nTime = len(ncfile.dimensions['time'])
        startTime = ' '.join(ncfile.variables['time'].base_date.split()[2:4])
        startTime = datetime.datetime.strptime(startTime, '%Y-%m-%d %H:%M:%S')
        time = timeArray.simulationToEpochTime(
            ncfile.variables['time'][:], startTime)
        return time

    def getZCoordinates(self, iStack):
        """Returns vertical coordinates for the given stack."""
        # TODO raise NotImplementedError('This method must be defined in the
        # derived class')
        ncfile = self.getNCFile(iStack, 'z')
        Z = ncfile.variables[getNCVariableName('z')][:]
        ncfile.close()
        return Z

    def getVerticalProfile(
            self,
            iStack,
            varStr,
            x,
            y,
            stationNames=None,
            horzInterp=None):
        """
        Extracts vertical profiles for the given locations from the given ncfile.

        Parameters
        ----------
        iStack : int
               Stack number of the netCDF file to process
        varStr : string
               Variable to extract
        x,y    : array_like (nPoints,)
               Coordinates of the points where to extract
        stationNames : list of strings, optional
               Names of the points for debugging

        Returns
        -------
        time  : array_like (nTime,)
              Time stamps of the extracted data in epoch format
        vals  : array_like (nPoints, nVert, nTime)
              Values of the extracted profiles. For 2d variables nVert=1.
              Missing data is filled with NaNs.
        zcoords : array_like (nPoints, nVert, nTime)
              Z-coordinates for the vertical profiles
        """

        if horzInterp is None:
            horzInterp = gridUtils.horizontalInterpolator(
                self.meshSearch2d, x, y, stationNames)

        ncfile = self.getNCFile(iStack)
        time = self.getTime(ncfile)

        V = ncfile.variables[getNCVariableName(varStr)]
        is3d = len(V.shape) == 3
        if not is3d:  # 2D variable
            if self.verbose:
                print '2D var', V.shape
            nodalValues = V[:][:, None, :]  # expand to (nNodes,nVert,nTime)
        else:
            if self.verbose:
                print '3D var', V.shape
            # NOTE reading the whole array may be inefficient?
            nodalValues = V[:]
        ncfile.close()

        vals = horzInterp.evaluateArray(nodalValues)
        if is3d:
            Z = self.getZCoordinates(iStack)
            zcoords = horzInterp.evaluateArray(Z)
        else:
            zcoords = np.zeros_like(vals)

        vals = np.ma.masked_invalid(vals)
        zcoords = np.ma.masked_invalid(zcoords)
        return time, vals, zcoords

    def getVerticalProfileForStacks(
            self, stacks, varStr, x, y, stationNames=None):
        """Extracts vertical profile for the given netcdf file stacks

        Parameters
        ----------
        stacks : list of int
               Stack numbers to process
        varStr : string
               Variable to extract
        x,y    : array_like (nPoints,)
               Coordinates of the points where to extract
        stationNames : list of strings, optional
               Names of the points for debugging

        Returns
        -------
        time  : array_like (nTime,)
              Time stamps of the extracted data in epoch format
        vals  : array_like (nPoints, nVert, nTime)
              Values of the extracted profiles. For 2d variables nVert=1.
              Missing data is filled with NaNs.
        zcoords : array_like (nPoints, nVert, nTime)
              Z-coordinates for the vertical profiles
        """
        # create interpolator object for recycling
        horzInterp = gridUtils.horizontalInterpolator(
            self.meshSearch2d, x, y, stationNames)

        time = []
        vals = []
        zcoords = []
        for stack in stacks:
            # extract for individual stacks
            try:
                ti, vi, zi = self.getVerticalProfile(
                    stack, varStr, x, y, horzInterp=horzInterp)
                time.append(ti)
                vals.append(vi)
                zcoords.append(zi)
            except Exception as e:
                print e
        # concatenate time axis
        time = np.concatenate(tuple(time), axis=0)  # (nTime,)
        vals = np.concatenate(tuple(vals), axis=2)  # (nProfiles,nVert,nTime)
        zcoords = np.concatenate(
            tuple(zcoords),
            axis=2)  # (nProfiles,nVert,nTime)
        time = np.ma.masked_invalid(time)
        vals = np.ma.masked_invalid(vals)
        zcoords = np.ma.masked_invalid(zcoords)
        return time, vals, zcoords

    def getTimeSeriesFromProfiles(
            self,
            vals,
            zcoords,
            z=None,
            k=None,
            zRelToSurf=None):
        """Interpolates vertical profiles in vertical at given depth.

        Parameters
        ----------
        vals  : array_like (nPoints, nVert, nTime)
              Values of the extracted profiles. For 2d variables nVert=1.
              Missing data is filled with NaNs.
        zcoords : array_like (nPoints, nVert, nTime)
              Z-coordinates for the vertical profiles
        z : float, array_like (nProfiles,), optional
          z coordinate where each vertical profile is evaluated. z coordinates
          increase upwards.
        k : int, array_like (nProfiles,), optional
          Instead of interpolating, take k-th nodal value from bottom.
          k=1 stands for bottom, k=-1 stands for surface
        zRelToSurf : bool, array_like (nProfiles,), optional
          If True z coordinate is taken depth below free surface instead of
          static z coordinate

        Returns
        -------
        vals : array_like (nProfiles,nTime,)
            Interpolated values
        z_actual : array_like (nProfiles,nTime,)
            The z coordinate at which the interpolation actually took place
        """
        vertInterp = verticalInterpolator(z, k, zRelToSurf)
        v, z_actual = vertInterp.evaluateArray(zcoords, vals)
        v = np.ma.masked_invalid(v)
        return v, z_actual

    def getTimeSeries(
            self,
            iStack,
            varStr,
            x,
            y,
            stationNames=None,
            z=None,
            k=None,
            zRelToSurf=None):
        """Extracts time series from the iStack netcdf file"""
        time, vprof, zcoords = self.getVerticalProfile(
            iStack, varStr, x, y, stationNames)
        if varStr in VARS2D:
            vals = vprof[:, 0, :]
            vals = np.ma.masked_invalid(vals)
            z_actual = np.zeros_like(vals)
        else:
            vals, z_actual = self.getTimeSeriesFromProfiles(
                vprof, zcoords, z, k, zRelToSurf)
        return time, vals, z_actual

    def getTimeSeriesForStacks(
            self,
            stacks,
            varStr,
            x,
            y,
            stationNames=None,
            z=None,
            k=None,
            zRelToSurf=None):
        """Extracts time series for the given stacks.

        Parameters
        ----------
        stacks : list of int
               Stack numbers to process
        varStr : string
               Variable to extract
        x,y    : array_like (nPoints,)
               Coordinates of the points where to extract
        stationNames : list of strings, optional
               Names of the points for debugging
        z : float, array_like (nProfiles,), optional
          z coordinate where each vertical profile is evaluated. z coordinates
          increase upwards.
        k : int, array_like (nProfiles,), optional
          Instead of interpolating, take k-th nodal value from bottom.
          k=1 stands for bottom, k=-1 stands for surface
        zRelToSurf : bool, array_like (nProfiles,), optional
          If True z coordinate is taken depth below free surface instead of
          static z coordinate

        Returns
        -------
        time  : array_like (nTime,)
              Time stamps of the extracted data in epoch format
        vals : array_like (nProfiles,nTime,)
            Interpolated values
        z_actual : array_like (nProfiles,nTime,)
            The z coordinate at which the interpolation actually took place
        """
        time, vprof, zcoords = self.getVerticalProfileForStacks(
            stacks, varStr, x, y, stationNames)
        if varStr in VARS2D:
            vals = vprof[:, 0, :]
            vals = np.ma.masked_invalid(vals)
            z_actual = np.zeros_like(vals)
        else:
            vals, z_actual = self.getTimeSeriesFromProfiles(
                vprof, zcoords, z, k, zRelToSurf)
        return time, vals, z_actual

    def getSlab(self, iStack, varStr, z=None, k=None, zRelToSurf=None):
        """
        Extracts a horizontal slice from the given ncfile.

        Parameters
        ----------
        iStack : int
               Stack number of the netCDF file to process
        varStr : string
               Variable to extract
        z      : float, array_like (nProfiles,), optional
               z coordinate where each vertical profile is evaluated. z coordinates
               increase upwards.
        k      : int, array_like (nProfiles,), optional
               Instead of interpolating, take k-th nodal value from bottom.
               k=-1 stands for bottom, k=0 stands for surface
        zRelToSurf : bool, array_like (nProfiles,), optional
               If True z coordinate is taken depth below free surface instead of
               static z coordinate

        Returns
        -------
        time  : array_like (nTime,)
              Time stamps of the extracted data in epoch format
        vals  : array_like (nPoints, nTime)
              Values of the extracted horizontal slice.
        """

        ncfile = self.getNCFile(iStack)
        time = self.getTime(ncfile)

        V = ncfile.variables[getNCVariableName(varStr)]
        is3d = len(V.shape) == 3
        if not is3d:  # 2D variable
            if self.verbose:
                print '2D var', V.shape
            vals = V[:]  # take
            z_actual = np.zeros_like(vals)
        else:
            if self.verbose:
                print '3D var', V.shape
            # NOTE reading the whole array may be inefficient?
            nodalValues = V[:]
            zcoords = self.getZCoordinates(iStack)
            nNodes, nZ, nTime = nodalValues.shape
            # nProfiles,nZ,nTime
            # TODO add vertical interpolation
            if k is not None:
                # bottom: k=1 kk=0, surface: k=0 kk=0
                #kk = k-1 if k>=0 else nZ+k
                kk = k
                vals = nodalValues[:, kk, :]
                z_actual = zcoords[:, kk, :]
            else:  # interpolate in vertical
                zArray = np.ones((nodalValues.shape[0],)) * z
                zRelArray = np.ones(
                    (nodalValues.shape[0],),
                    dtype=int) * zRelToSurf
                vertInterp = verticalInterpolator(zArray, None, zRelArray)
                vals, z_actual = vertInterp.evaluateArray(zcoords, nodalValues)
            vals = np.ma.masked_invalid(vals)
        ncfile.close()

        return time, vals, z_actual

    def getSlabForStacks(
            self,
            stacks,
            varStr,
            z=None,
            k=None,
            zRelToSurf=None):
        """Returns slab for the given stacks"""
        time = []
        vals = []
        zcoords = []
        for stack in stacks:
            # extract for individual stacks
            try:
                ti, vi, zi = self.getSlab(stack, varStr, z, k, zRelToSurf)
                time.append(ti)
                vals.append(vi)
                zcoords.append(zi)
            except Exception as e:
                print e
        if len(time) == 0:
            raise Exception('Extraction Failed: no time steps were retrieved')
        # concatenate time axis
        time = np.concatenate(tuple(time), axis=0)  # (nTime,)
        vals = np.concatenate(tuple(vals), axis=1)  # (nProfiles,nTime)
        zcoords = np.concatenate(tuple(zcoords), axis=1)  # (nProfiles,nTime)
        time = np.ma.masked_invalid(time)
        vals = np.ma.masked_invalid(vals)
        zcoords = np.ma.masked_invalid(zcoords)
        return time, vals, zcoords

    def getStacks(
            self,
            startTime,
            endTime,
            ncfile=None,
            firstPointIncluded=True):
        """Returns a list of file stack numbers that covers the given
        time period [startTime,endTime].

        Simulation start time is read from the netcdf header.
        """
        # deduce correct stacks
        ncfileGiven = ncfile is not None
        if self.verbose:
            print 'Reading header'
        if not ncfileGiven:
            ncfile = self.getNCFile()

        if not self.headerIsRead:
            self.readHeader(ncfile)
        time = self.getTime(ncfile)
        nSpool = 1  # number of exports in each file
        exportDt = 15 * 60  # time interval between exports
        spoolDt = nSpool * exportDt

        startDelta = (startTime - self.simulationStartTime).total_seconds()
        endDelta = (endTime - self.simulationStartTime).total_seconds()
        if not firstPointIncluded and nSpool > 1:
            startDelta -= exportDt
        if firstPointIncluded and nSpool > 1:
            endDelta += exportDt
        startStack = int(np.floor(startDelta / spoolDt)) + 1
        endStack = int(np.ceil(endDelta / spoolDt))
        if not ncfileGiven:
            ncfile.close()
        return range(startStack, endStack + 1)


class slimExtract(slimExtractBase):
    """
    This class contains only high-level extraction routines and returns the
    data in dataContainer with metadata.
    """

    def __init__(self, path, var=None, verbose=False):
        slimExtractBase.__init__(self, path, var, verbose)

    def extractTimeSeries(
            self,
            startTime,
            endTime,
            var,
            staX,
            staY,
            stationNames,
            staZ=None,
            k=None,
            zRelToSurf=None):
        """Extracts time series for the given time range."""
        stacks = self.getStacks(startTime, endTime)
        time, vals, actualZ = self.getTimeSeriesForStacks(
            stacks, var, staX, staY, stationNames, staZ, k, zRelToSurf)

        # build dataContainer for each station
        dcs = []
        for iSta in range(len(staX)):
            data = vals[iSta, :]
            # remove nans
            goodIx = np.logical_not(data.mask)
            if not goodIx.any():
                # all bad data
                print 'all bad data', stationNames[iSta]
                continue
            # NOTE data,time must be ndarray not masked array
            data = np.reshape(np.array(data[goodIx]), (1, 1, -1))
            t = np.array(time[goodIx])

            ta = timeArray.timeArray(t, 'epoch')
            meta = {}
            meta['location'] = stationNames[iSta]
            meta['instrument'] = 'model'
            meta['variable'] = var
            alongSLevel = False  # FIXME
            if alongSLevel:
                meta['slevel'] = kLevel
            else:
                zSign = 1 if zRelToSurf else -1  # zRelToSurf => depth below surface
                zTarget = 0.0 if var in VARS2D else staZ[iSta]
                msldepth = str(int(round(zSign * zTarget * 100)))
                meta['bracket'] = 'F' if zRelToSurf else 'A'
                meta['msldepth'] = msldepth
            meta['dataType'] = 'timeseries'
            z = np.mean(actualZ[iSta, :])
            x = staX[iSta]
            y = staY[iSta]
            dc = dataContainer(
                '', ta, x, y, z, data, fieldNameList.get(
                    var, [var]), coordSys='spcs', metaData=meta)
            dcs.append(dc)
        return dcs

    def extractVerticalProfile(
            self,
            startTime,
            endTime,
            var,
            staX,
            staY,
            stationNames=None):
        """Extracts vertical profiles for the given time period."""
        # TODO TEST
        stacks = self.getStacks(startTime, endTime)
        time, vals, zcoords = self.getVerticalProfileForStacks(
            stacks, var, staX, staY, stationNames)

        # time  : array_like (nTime,)
        # vals  : array_like (nPoints, nVert, nTime)
        # zcoords : array_like (nPoints, nVert, nTime)
        # build dataContainer for each station
        dcs = []
        for iSta in range(len(staX)):
            staName = '' if stationNames is None else stationNames[iSta]
            staStr = '{x:f} {y:f} {name:s}'.format(
                x=staX[iSta], y=staY[iSta], name=staName)
            # remove time steps with all bad values
            goodIxTime = np.logical_and(
                ~np.all(
                    vals[
                        iSta, :, :].mask, axis=0), ~np.all(
                    zcoords[
                        iSta, :, :].mask, axis=0))
            goodIxTime = np.nonzero(goodIxTime)[0]
            v = vals[iSta, :, :][:, goodIxTime]
            z = zcoords[iSta, :, :][:, goodIxTime]
            t = time[goodIxTime]
            print '1', z.shape
            # remove masked (below bottom) part (nVert,nTime) ->
            # (nGoodVert,nTime)
            goodIxVert = np.logical_and(~np.any(v.mask, axis=1),
                                        ~np.any(z.mask, axis=1))

            if not goodIxVert.any() or not goodIxTime.any():
                # all bad data
                print not goodIxVert.any(), not goodIxTime.any()
                print 'all bad data: ', staStr
                continue
            v = v[goodIxVert, :]
            z = z[goodIxVert, :]
            if v.mask.any() or z.mask.any() or t.mask.any():
                print v.mask.any(), z.mask.any(), t.mask.any()
                print v.mask
                print z.mask
                print t.mask
                raise Exception('bad values remain: ' + staStr)
            # to (nGoodVert,1,nTime)
            data = v[:, None, :]
            ta = timeArray.timeArray(np.array(t), 'epoch')
            # to (nGoodVert,nTime)
            nZ = z.shape[0]
            x = staX[iSta] * np.ones((nZ,))
            y = staY[iSta] * np.ones((nZ,))
            meta = {}
            meta['location'] = stationNames[iSta]
            meta['instrument'] = 'model'
            meta['bracket'] = 'A'
            meta['variable'] = var
            meta['dataType'] = 'profile'
            dc = dataContainer(
                '', ta, x, y, z, data, fieldNameList.get(
                    var, [var]), coordSys='spcs', metaData=meta)
            dcs.append(dc)
        return dcs

    def extractTransect(self, startTime, endTime, var, staX, staY, transName):
        """Extracts a transect for the given (x,y) points and time range."""
        stacks = self.getStacks(startTime, endTime)
        staX = np.array(staX)
        staY = np.array(staY)
        time, vals, zcoords = self.getVerticalProfileForStacks(
            stacks, var, staX, staY)

        # time  : array_like (nTime,)
        # vals  : array_like (nPoints, nVert, nTime)
        # zcoords : array_like (nPoints, nVert, nTime)

        # reorganize transect in trasect format
        data = []
        X = []
        Y = []
        Z = []
        for iSta in range(len(staX)):
            staStr = '{ix:d} {x:f} {y:f}'.format(
                ix=iSta, x=staX[iSta], y=staY[iSta])

            z = zcoords[iSta, :, :]
            v = vals[iSta, :, :]
            # remove masked (below bottom) part (nVert,nTime) ->
            # (nGoodVert,nTime)
            goodIxVert = np.logical_and(~np.all(v.mask, axis=1),
                                        ~np.all(z.mask, axis=1))
            if not goodIxVert.any():
                # all bad data
                print 'all bad data: ', staStr
                continue
            z = z[goodIxVert, :].filled(np.nan)
            v = v[goodIxVert, :].filled(np.nan)

            x = staX[iSta] * np.ones((z.shape[0],))
            y = staY[iSta] * np.ones((z.shape[0],))
            data.append(v)  # (nVert,nTime)
            X.append(x)  # (nVert,)
            Y.append(y)  # (nVert,)
            Z.append(z)  # (nVert,nTime)
        # concatenate
        X = np.concatenate(tuple(X), axis=0)
        Y = np.concatenate(tuple(Y), axis=0)
        Z = np.concatenate(tuple(Z), axis=0)
        data = np.concatenate(tuple(data), axis=0)
        # reshape for dataContainer
        # from (nVert*nSta,nTime) to (nVert*nSta,1,nTime)
        data = data[:, None, :]

        # build dataContainer
        ta = timeArray.timeArray(time, 'epoch')
        meta = {}
        meta['location'] = transName
        meta['instrument'] = 'model'
        meta['bracket'] = 'A'
        meta['variable'] = var
        meta['dataType'] = 'transect'
        dc = dataContainer('', ta, X, Y, Z, data, fieldNameList.get(
            var, [var]), coordSys='spcs', metaData=meta, acceptNaNs=True)
        return dc

    def extractTransectForBPFile(
            self,
            startTime,
            endTime,
            var,
            bpFile,
            transName):
        """Extracts a transect for the given build point ASCII file and time range."""
        bpObj = buildPoints.BuildPoint()
        bpObj.readFileFromDisk(bpFile)
        return self.extractTransect(startTime, endTime, var,
                                    bpObj.getX(), bpObj.getY(), transName)

    def extractTrack(self):
        """Extracts a (x,y,z,t) track for the given time range."""
        raise NotImplementedError('This feature has not been implemented yet.')

    def extractSlab(
            self,
            startTime,
            endTime,
            name,
            var,
            z=None,
            k=None,
            zRelToSurf=None):
        """Extracts a horiontal slice for the given time range."""
        stacks = self.getStacks(startTime, endTime)
        time, vals, zcoords = self.getSlabForStacks(
            stacks, var, z, k, zRelToSurf)
        ta = timeArray.timeArray(time, 'epoch')
        data = vals[:, None, :]
        data = data.filled(np.nan)
        connectivity = self.faceNodes
        x = self.nodeX
        y = self.nodeY
        msldepth = ''
        if k is not None:
            msldepth = 'slev' + str(k)
            z = zcoords[:, 0]  # FIXME include time in z coords ?
        else:
            zSign = 1 if zRelativeToSurf else -1  # zRelToSurf => depth below surface
            msldepth = str(int(round(zSign * z * 100)))
            zArray = z * np.ones_like(x)

        # make meshContainer
        meta = {}
        meta['dataType'] = 'slab'
        meta['location'] = name
        meta['instrument'] = 'model'
        meta['variable'] = var
        if k is not None:
            meta['slevel'] = k
        else:
            meta['bracket'] = 'F' if zRelativeToSurf else 'A'
            meta['msldepth'] = msldepth
        mc = meshContainer.meshContainer(
            '', ta, x, y, z, data, connectivity, fieldNameList[var],
            coordSys='spcs', metaData=meta)
        return mc

    def extractForXYZ(self, dataDir, var, startTime, endTime, x, y, z=None,
                      stationNames=None, profile=False, zRelToSurf=False,
                      stacks=None, verbose=False):
        """
        Extracts time series for given variable from stations defined by x,y,z.
        If profile=True, will extract profiles instead (z is ignored).
        """
        varStr, fileTypeStr = splitVarToFileType(var)
        if not profile and z is None:
            raise Exception('z coordinates must be provided')
        try:
            if profile:
                dcs = self.extractVerticalProfile(
                    startTime, endTime, varStr, x, y, stationNames, stacks=stacks)
            else:
                dcs = self.extractTimeSeries(
                    startTime, endTime, varStr, x, y, stationNames, z,
                    zRelToSurf=zRelToSurf)
            print ' * extracted'
        except Exception as e:
            print ' * extraction failed'
            traceback.print_exc(file=sys.stdout)
            dcs = []
        for dc in dcs:
            if profile:
                print ' '.join([dc.getMetaData('location'), dc.getMetaData('variable')])
            else:
                print ' '.join([dc.getMetaData('location'), dc.getMetaData('bracket'),
                                dc.getMetaData('msldepth'), dc.getMetaData('variable')])
        return dcs

    def extractForStations(self, dataDir, var, stationFile, startTime, endTime,
                           profile=False, stacks=None, verbose=False):
        """
        Extracts time series for given variable from stations defined in stationFile.
        """
        # read station file, allow duplicate stations (with different depth)
        if not os.path.isfile(stationFile):
            raise Exception('File does not exist: ' + stationFile)

        if profile:
            csvReader = csvStationFile.csvStationFile()
            csvReader.readFromFile(stationFile)
            tuples = csvReader.getTuples()  # all entries (loc,x,y)
            stationNames = [t[0] for t in tuples]
            x = np.array([t[1] for t in tuples])
            y = np.array([t[2] for t in tuples])
            z = None
            print ' *** extracting profiles for stations *** '
            for i, s in enumerate(stationNames):
                print s, x[i], y[i]
            return self.extractForXYZ(
                dataDir,
                var,
                startTime,
                endTime,
                x,
                y,
                z,
                stationNames,
                profile,
                False,
                stacks=stacks,
                verbose=verbose)

        # not profile, depths defined in stationFile
        csvReader = csvStationFile.csvStationFileWithDepth()
        csvReader.readFromFile(stationFile)
        tuples = csvReader.getTuples()  # all entries (loc,x,y,z,zType,var)
        stationNames = [t[0] for t in tuples]
        x = np.array([t[1] for t in tuples])
        y = np.array([t[2] for t in tuples])
        z = np.array([t[3] for t in tuples])
        zRelToSurf = np.array([t[4] == 'depth' for t in tuples], dtype=bool)

        dcList = []
        for zIsDepth in [True, False]:
            # filter for z coordinate/depth cases
            ix = np.nonzero(zRelToSurf == zIsDepth)[0]
            if len(ix) == 0:
                continue
            x_filt = x[ix]
            y_filt = y[ix]
            z_filt = z[ix]
            stationNames_filt = [stationNames[i] for i in ix]
            print ' *** extracting for stations *** '
            for i, s in enumerate(stationNames_filt):
                print s, x_filt[i], y_filt[i], z_filt[i], zRelToSurf[i]
            dcs = self.extractForXYZ(dataDir, var, startTime, endTime,
                                     x_filt, y_filt, z_filt, stationNames_filt,
                                     profile, zIsDepth, stacks=stacks,
                                     verbose=verbose)
            dcList.extend(dcs)

        return dcList


def splitVarToFileType(var):
    """Splits 'varname.ext' to 'var','ext'.
    If returns None as extension if cannot split."""
    if len(var.split('.')) == 2:
        return var.split('.')
    else:
        return var, None
