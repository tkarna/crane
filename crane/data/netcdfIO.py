"""
Methods for reading/writing dataContainer/meshContainer data from/to netCDF files.

Tuomas Karna 2013-01-23
"""

#-------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------

import os
import numpy as np
from netCDF4 import Dataset as NetCDFFile

import crane
from crane.data import timeArray

#-------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------

# variable strings used in netCDF files
ncVariableNames = {'elev': 'elevation',
                   'temp': 'water_temperature',
                   'salt': 'water_salinity',
                   'pres': 'water_pressure',
                   'cond': 'water_electrical_conductivity',
                   }
# inverse dictionary
ncVariableNamesInv = dict((ncVariableNames[k], k) for k in ncVariableNames)

ncVariableUnits = {'elev': 'm',
                   'temp': 'C',
                   'salt': 'psu',
                   'pres': 'Pa',
                   'cond': 'S/m',
                   }

# time format description for netCDF format
ncTimeUnits = {'epoch': 'seconds since 1970-01-01 00:00:00-00',  # UTC
               'corie': 'days since 1995-12-31 00:00:00 PST',
               'simulation': 'seconds since simulation started',
               }
# inverse dictionary
ncTimeUnitsInv = dict((ncTimeUnits[k], k) for k in ncTimeUnits)


def createDirectory(path):
    if path == '':
        return
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise Exception('file with same name exists', path)
    else:
        os.makedirs(path)
    return path

#-------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------


class netcdfIO(object):
    """A main class for reading/writing data/meshContainer data in netCDF format."""

    def __init__(self, filename):
        """Init instance with the given filename (for reading/writing)."""
        self.filename = filename

    def readHeader(self, startTime=None, endTime=None, verbose=True):
        """Reads file header.
        """
        # if verbose: print 'reading header '+self.filename
        if not os.path.isfile(self.filename):
            raise Exception('file not found {0:s}'.format(self.filename))
        ncfile = NetCDFFile(self.filename, 'r')
        return netcdfIO.readNetCDFFile(
            ncfile, startTime, endTime, headerOnly=True)

    def read(
            self,
            startTime=None,
            endTime=None,
            includeEnd=False,
            verbose=True):
        """Reads nc file header and contents."""
        if verbose:
            print 'reading ' + self.filename
        if not os.path.isfile(self.filename):
            raise Exception('file not found {0:s}'.format(self.filename))
        ncfile = NetCDFFile(self.filename, 'r')
        return netcdfIO.readNetCDFFile(
            ncfile, startTime, endTime, includeEnd=includeEnd)

    def readToDataContainer(
            self,
            startTime=None,
            endTime=None,
            includeEnd=False,
            verbose=True,
            dtype=None):
        """Reads file and returns data in data/meshContainer."""
        out = self.read(startTime, endTime, includeEnd, verbose=verbose)
        descr, ta, x, y, z, data, fns, cSys, conn, meta, bnds = out
        if conn is None:
            # NOTE import here to avoid cyclic import issues
            import crane.data.dataContainer as dataContainer
            dc = dataContainer.dataContainer(descr, ta,
                                             x, y, z, data,
                                             fns, cSys,
                                             acceptNaNs=True,
                                             dtype=dtype)
        else:
            import crane.data.meshContainer as meshContainer
            dc = meshContainer.meshContainer(descr, ta,
                                             x, y, z,
                                             data, conn, fns,
                                             cSys, dtype=dtype)
            if bnds:
                for bnd in bnds:
                    dc.addBoundary(bnd)
        dc.setMetaData(meta)
        return dc

    def saveDataContainer(
            self,
            dataContainer,
            dtype=np.float64,
            overwrite=False,
            compress=False,
            digits=None):
        """Saves dataContainer in netCDF format.
        dtype sets data type to any numpy format (except float16)
        overwrite=True forces writing over existing file, otherwise error is raised
        compress=True sets lossless zlib compression for netcdf4 dataset
        digits rounds data to given decimals before storing
        """
        print 'writing to ' + self.filename
        if os.path.exists(self.filename):
            if os.path.isfile(self.filename):
                if not overwrite:
                    raise IOError(
                        'file already exists: {0:s}'.format(
                            self.filename))
                else:
                    print ' file exists, overwriting ...'
                    os.remove(self.filename)
            else:
                raise IOError(
                    'unable to create file: {0:s}'.format(
                        self.filename))
        createDirectory(os.path.split(self.filename)[0])
        ncfile = NetCDFFile(self.filename, 'w')
        self.populateNetCDFObject(
            ncfile, dataContainer, dtype, compress, digits)
        # close the file.
        ncfile.close()

    @staticmethod
    def readNetCDFFile(ncfile, startTime=None, endTime=None, headerOnly=False,
                       includeEnd=False, dtype=np.float32):
        """Helper function that parses and manipulates netCDF file data"""

        # Workaround for 'bool' object has no attribute 'any' error
        # (missing_value in ncfile should not be a string)
        # v = file.variables[variable]
        # v.set_auto_maskandscale(False)
        description = ''
        # read metaData
        metaDataNames = ncfile.ncattrs()
        metaData = dict((k, ncfile.getncattr(k)) for k in metaDataNames)
        if 'description' in metaData.keys():
            # try to guess metaData from description string (deprecated)
            # grays.160.A.model.temp
            description = metaData['description']
            words = metaData['description'].split('.')
            #words = metaData.pop('description').split('.')
            if len(words) == 5:
                metaData.setdefault('location', words[0])
                msldepth = words[1]
                metaData.setdefault('bracket', words[2])
                metaData.setdefault('instrument', words[3])
                metaData.setdefault('variable', words[4])
                if msldepth == 'prof':
                    metaData.setdefault('dataType', 'profile')
                elif msldepth == 'profiler':
                    metaData.setdefault('dataType', 'profiler')
                elif msldepth[:4] == 'slev':
                    metaData.setdefault('dataType', 'slab')
                elif msldepth == 'tran':
                    metaData.setdefault('dataType', 'transect')
                elif msldepth == 'trck':
                    metaData.setdefault('dataType', 'track')
                else:
                    metaData.setdefault('dataType', 'timeseries')
                    metaData.setdefault('msldepth', msldepth)
        nTime = len(ncfile.dimensions['time'])
        nX = len(ncfile.dimensions['x']
                 ) if 'x' in ncfile.dimensions.keys() else 1
        nY = len(ncfile.dimensions['y']
                 ) if 'y' in ncfile.dimensions.keys() else 1
        nZ = len(ncfile.dimensions['z']
                 ) if 'z' in ncfile.dimensions.keys() else 1
        nData = len(ncfile.dimensions[
            'data']) if 'data' in ncfile.dimensions.keys() else nX
        vars = ncfile.variables
        timeVar = vars.pop('time')
        timeFormat = ncTimeUnitsInv[timeVar.units]
        ta = timeArray.timeArray(timeVar[:], timeFormat, acceptDuplicates=True)
        # infer correct time indices
        timeIx = ta.getRangeIndices(startTime, endTime, includeEnd=includeEnd)
        if len(timeIx) == 0:
            raise Exception(
                'No time stamps found for the given bounds {0:s} - {1:s}'.format(
                    str(startTime), str(endTime)))
        ta = timeArray.timeArray(
            timeVar[timeIx],
            timeFormat,
            acceptDuplicates=True)
        if headerOnly:
            return metaData, ta  # nX,nY,nZ ??
        coordSys = ''
        if 'coordSys' in metaData:
            coordSys = metaData.pop('coordSys')
        x = y = z = 0
        if 'x' in ncfile.variables:
            xVar = vars.pop('x')
            xDependsOnTime = False if len(xVar.shape) == 1 else True
            if xDependsOnTime:
                x = xVar[:, timeIx]
            else:
                x = xVar[:]
        if 'y' in ncfile.variables:
            yVar = vars.pop('y')
            yDependsOnTime = False if len(yVar.shape) == 1 else True
            if yDependsOnTime:
                y = yVar[:, timeIx]
            else:
                y = yVar[:]
        if 'z' in ncfile.variables:
            zVar = vars.pop('z')
            zDependsOnTime = False if len(zVar.shape) == 1 else True
            if zDependsOnTime:
                z = zVar[:, timeIx]
            else:
                z = zVar[:]
        if dtype is not None:
            x = x.astype(dtype)
            y = y.astype(dtype)
            z = z.astype(dtype)
        fieldNames = []
        # possible mesh connectivity (for meshContainer)
        connVar = vars.pop('connectivity', None)
        connectivity = connVar[:] if connVar else None

        # possible mesh boundary information
        bndVars = [vname for vname in vars if vname.find('bnd_nodes_') == 0]
        if bndVars:
            boundaries = []
            for bndVar in sorted(bndVars):
                v = vars.pop(bndVar)
                bnd = crane.data.meshContainer.meshBoundary(
                    v.getncattr('type'), v.getncattr('tag'), v[:])
                boundaries.append(bnd)
        else:
            boundaries = None

        # accept only variables with correct size
        goodVars = [
            v for v in vars if vars[v].shape == (
                nData, nTime) or vars[v].shape == (
                nTime,)]
        data = np.zeros((nData, len(goodVars), len(timeIx)), dtype=dtype)
        for i, var in enumerate(goodVars):
            fieldNames.append(ncVariableNamesInv.get(var, var))
            if vars[var].shape == (nData, nTime):
                data[:, i, :] = vars[var][:, timeIx]
            elif vars[var].shape == (nTime,):
                data[:, i, :] = vars[var][timeIx]

        return description, ta, x, y, z, data, fieldNames, coordSys, connectivity, metaData, boundaries

    def populateNetCDFObject(
            self,
            ncfile,
            dc,
            dtype=None,
            compress=False,
            digits=None):
        """Helper function that creates all necessary variables and assings the values"""
        # wrappers to handle older scipy
        def createDimension(file, tag, dim):
            try:  # for scipy >= 0.9
                return file.createDimension(tag, dim)
            except:
                return file.create_dimension(tag, dim)

        def createVariable(file, tag, typ, dim, zlib=False, digits=None):
            try:  # for scipy >= 0.9
                if digits is not None:
                    return file.createVariable(tag, typ, dim, zlib=zlib,
                                               least_significant_digit=digits)
                else:
                    return file.createVariable(tag, typ, dim, zlib=compress)
            except:
                if digits is not None:
                    return file.create_variable(tag, typ, dim, zlib=zlib,
                                                least_significant_digit=digits)
                else:
                    return file.create_variable(tag, typ, dim, zlib=compress)
        if dc.description:
            ncfile.description = dc.description
        # add all metaData
        ncfile.setncatts(dc.metaData)
        # create the x and y dimensions.
        nTime = len(dc.time)
        createDimension(ncfile, 'time', nTime)
        nXYZ = dc.x.shape[0]
        createDimension(ncfile, 'x', nXYZ)
        createDimension(ncfile, 'y', nXYZ)
        createDimension(ncfile, 'z', nXYZ)
        nData = dc.data.shape[0]
        dataAtNodes = nData == nXYZ
        if not dataAtNodes:
            createDimension(ncfile, 'data', nData)
        # create the variable
        # first argument is name of variable, second is datatype, third is
        # a tuple with the names of dimensions.
        # NOTE time must be double precision to avoid round-off in seconds
        timeVar = createVariable(
            ncfile, 'time', np.dtype(
                np.float64).char, ('time', ), compress)
        timeVar[:] = dc.time.asEpoch().array.astype(np.float64)
        timeVar.units = ncTimeUnits['epoch']
        if dtype is None:
            x = dc.x
            y = dc.y
            z = dc.z
        else:
            x = dc.x.astype(dtype)
            y = dc.y.astype(dtype)
            z = dc.z.astype(dtype)
        if dc.xDependsOnTime:
            xVar = createVariable(
                ncfile, 'x', x.dtype.char, ('x', 'time', ), compress)
            xVar[:, :] = x
        else:
            xVar = createVariable(ncfile, 'x', x.dtype.char, ('x',), compress)
            xVar[:] = x
        if dc.yDependsOnTime:
            yVar = createVariable(
                ncfile, 'y', y.dtype.char, ('y', 'time', ), compress)
            yVar[:, :] = y
        else:
            yVar = createVariable(ncfile, 'y', y.dtype.char, ('y', ), compress)
            yVar[:] = y
        if dc.zDependsOnTime:
            zVar = createVariable(
                ncfile, 'z', y.dtype.char, ('z', 'time', ), compress)
            zVar[:, :] = z
        else:
            zVar = createVariable(ncfile, 'z', z.dtype.char, ('z', ), compress)
            zVar[:] = z
        coordSysLabels = {
            'spcs': (
                'SPCS 3601 x coordinate', 'SPCS 3601 y coordinate'), 'lonlat': (
                'longitude', 'latitude'), 'utm': (
                'UTM x coordinate', 'UTM y coordinate')}
        coordSysUnits = {'spcs': ('m', 'm'),
                         'lonlat': ('degrees east', 'degrees north'),
                         'utm': ('m', 'm')}
        coordSys = dc.coordSys
        if not coordSys:
            coordSys = 'spcs'
        ncfile.setncatts({'coordSys': coordSys})
        if coordSys in coordSysUnits:
            xVar.units, yVar.units = coordSysUnits[coordSys]
        if coordSys in coordSysLabels:
            xVar.long_name, yVar.long_name = coordSysLabels[coordSys]
        zVar.units = 'm'
        zVar.long_name = 'elevation versus NGVD29 datum'
        if (hasattr(dc, 'metaData') and dc.metaData is not None and
                'bracket' in dc.metaData and dc.metaData['bracket'] == 'F'):
            zVar.long_name = 'meters below free surface'

        dataVars = []
        for i, var in enumerate(dc.fieldNames):
            d_array = dc.data[
                :, i, :] if dtype is None else dc.data[
                :, i, :].astype(dtype)
            if dataAtNodes:
                dv = createVariable(
                    ncfile, ncVariableNames.get(var, var),
                    d_array.dtype.char, ('x', 'time',),
                    compress, digits)
            else:
                dv = createVariable(
                    ncfile, ncVariableNames.get(var, var),
                    d_array.dtype.char, ('data', 'time',),
                    compress, digits)
            dataVars.append(dv)
            dv[:, :] = d_array
            if var in ncVariableUnits:
                dv.units = ncVariableUnits[var]

        if hasattr(dc, 'connectivity'):
            # add mesh connectivity
            nElem, nNodes = dc.connectivity.shape
            createDimension(ncfile, 'nElements', nElem)
            createDimension(ncfile, 'nElemNodes', nNodes)
            connVar = createVariable(ncfile, 'connectivity', np.dtype(
                np.int64).char, ('nElements', 'nElemNodes',), compress)
            connVar[:, :] = dc.connectivity.astype(np.int64)

        if hasattr(dc, 'boundaries'):
            # add mesh boundary information
            for i, bnd in enumerate(dc.boundaries, start=1):
                bndDimStr = 'nBnd' + str(i)
                bndNameStr = 'bnd_nodes_' + str(i)
                # add array of bnd nodes
                createDimension(ncfile, bndDimStr, bnd.nodes.shape[0])
                bndVar = createVariable(
                    ncfile, bndNameStr, np.dtype(
                        np.int64).char, (bndDimStr,), compress)
                #bndVar = createVariable(ncfile,'bnd_nodes_'+bnd.tag,np.dtype(dtype).char,(bndDimStr,))
                bndVar[:] = bnd.nodes.astype(np.int64)
                # add string tags
                attrDict = {'type': bnd.type, 'tag': bnd.tag}
                bndVar.setncatts(attrDict)
