#!/usr/bin/python
"""
Extract station data using the efficient SELFE extract_mod python module.

Basic python usage:
  ee = extractStation(dataDir,var,profile=False,modelCoordSys='spcs')
  ee.setStations( stationNames, x,y,z, zRelativeToSurf=False )
  dcs = ee.extractDates( startTime, endTime )

Or:
  extractForOfferings( dataDir, var, offerings, startTime, endTime, profile=False,
                          modelCoordSys='spcs', stationFile=None )

Examples for commandline interface:

# extract elevation from /some/run/outputs and store netCDF files to outputDir
python extractStation.py -v elev -d /some/run/outputs -o outputDir -s 2011-05-11 -e 2011-07-20

# extract vertical profile of temperature from /some/run/outputs and store netCDF files to outputDir
python extractStation.py -v temp -p -d /some/run/outputs -o outputDir -s 2011-05-11 -e 2011-07-20

# same as above except for a model that is UTM coordinates (e.g. deb28 grid)
python extractStation.py -v elev -d /some/run/outputs -o outputDir -s 2011-05-11 -e 2011-07-20 -c utm

# working example
python extractStation.py -d /home/workspace/ccalmr42/karnat/runs/db28dev/run03/outputs/ -v elev -o tmp_net -s 2002-5-17 -e 2002-5-19 -c utm

Tuomas Karna 2012-10-10
"""
import os
import sys
import datetime
import glob
import time as timeMod
# set Ctrl-C to default signal (terminates fortran routines immediately)
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import numpy as np

from crane.data import dataContainer
from crane.data import timeArray
from crane.data import loadHindcastStations
from crane.files import csvStationFile
from crane.physicalVariableDefs import addTracers

# TODO extracting legacy selfe format is now obsolete?
#import extract_mod

#------------------------------------------------------------------------------
# Constants
#------------------------------------------------------------------------------

VALID_MIN = -89


def extractForXYZ(dataDir, var, stationFile, startTime, endTime, profile=False,
                  modelCoordSys='spcs', tracers=None, ntracers=0):
    """
    Extracts all stations and depths defined in the stationFile.

    Args:
      dataDir -- (str) path to model outputs directory
      var     -- (str) variable to extract, e.g. 'elev'. 'temp'
      stationFile   -- (str) a stations.csv file for reading x,y,z coordinates
      startTime -- (datetime) first time stamp of extraction
      endTime   -- (datetime)  last time stamp of extraction
      profile -- (bool) if true, extracts vertical profile instead of a value at z
      modelCoordSys -- (str) either 'spcs' (default) or 'utm'
      tracers -- (str) Tracer model type {'sed', 'oxy', 'generic'}
      ntracers -- (int) Number of tracers if using sed or generic tracer models

    Returns:
      dcList    -- (list of dataContainer) all dataContainers in a list
    """

    if profile:
        csvReader = csvStationFile.csvStationFile()
        csvReader.readFromFile(stationFile)
        tuples = csvReader.getTuples()  # all entries (loc,x,y)
        stationNames = [t[0] for t in tuples]
        x = np.array([t[1] for t in tuples])
        y = np.array([t[2] for t in tuples])
    else:
        csvReader = csvStationFile.csvStationFileWithDepth()
        csvReader.readFromFile(stationFile)
        tuples = csvReader.getTuples()  # all entries (loc,x,y,z,zType,var)
        stationNames = [t[0] for t in tuples]
        x = np.array([t[1] for t in tuples])
        y = np.array([t[2] for t in tuples])
        z = np.array([t[3] for t in tuples])
        zRelToSurf = np.array([t[4] == 'depth' for t in tuples], dtype=bool)

    if tracers is not None:
        addTracers(tracers, varList, numTracers=ntracers)

    # extract
    dcList = []
    ee = extractStation(
        dataDir,
        var,
        profile=profile,
        modelCoordSys=modelCoordSys)
    if profile:
        # discard depth information, use only station name
        print ' *** extracting profiles for stations *** '
        for i, s in enumerate(stationNames):
            print s, x[i], y[i]

        # execute
        ee.setStations(stationNames, x, y)
        dcs = ee.extractDates(startTime, endTime)
        dcList.extend(dcs)

    else:
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

            # execute
            ee.setStations(
                stationNames_filt,
                x_filt,
                y_filt,
                z_filt,
                zRelativeToSurf=zIsDepth)
            try:
                dcs = ee.extractDates(startTime, endTime)
                print ' * extracted'
            except Exception as e:
                print ' * extraction failed'
                print e
                dcs = []
            for dc in dcs:
                print (dc.getMetaData('location'), dc.getMetaData('variable'),
                       dc.getMetaData('msldepth'))
            dcList.extend(dcs)

    return dcList


def extractForOfferings(
        dataDir,
        var,
        offerings,
        startTime,
        endTime,
        profile=False,
        modelCoordSys='spcs',
        stationFile=None):
    """
    Extracts all stations and depths as in the list of offerings.
    If offerings list is not provided, it will be fetched from the database.

    Args:
      dataDir -- (str) path to model outputs directory
      var     -- (str) variable to extract, e.g. 'elev'. 'temp'
      offerings -- (list of str) list of offerins from the database
                   each string is station.msldepth.bracket.instrument[.var]
      startTime -- (datetime) first time stamp of extraction
      endTime   -- (datetime)  last time stamp of extraction
      profile -- (bool) if true, extracts vertical profile instead of a value at z
      modelCoordSys -- (str) either 'spcs' (default) or 'utm'
      stationFile   -- (str) a stations.csv file for reading coordinates

    Returns:
      dcList    -- (list of dataContainer) all dataContainers in a list
    """

    # read station file
    staReader = csvStationFile.csvStationFile()
    staReader.readFromFile(stationFile)

    # screen possible duplicates in the offerings (e.g. instrument can be
    # ignored)
    uniqueOffrs = {}
    for o in offerings:
        if var in ['elev']:
            # 2D variable, depth makes no difference
            key = (o['location'])
        else:
            key = (o['location'], o['msldepth'], o['bracket'])
        uniqueOffrs[key] = o
    offerings = uniqueOffrs.values()

    # extract
    dcList = []
    ee = extractStation(
        dataDir,
        var,
        profile=profile,
        modelCoordSys=modelCoordSys)
    if profile:
        # discard msldepth information, use only station name
        stationNames = set()
        for o in offerings:
            stationNames.add(o['location'])
        # omit stations that are missing from sta file
        stationNames = list(set(staReader.getStations()).intersection(
            set(stationNames)))
        x = np.array([staReader.getX(s) for s in stationNames])
        y = np.array([staReader.getY(s) for s in stationNames])
        print ' *** extracting profiles for stations *** '
        for s in stationNames:
            print s

        # execute
        ee.setStations(stationNames, x, y)
        dcs = ee.extractDates(startTime, endTime)
        dcList.extend(dcs)

    else:
        # split offerings to free surface ones ( bracket = 'F' ) and others
        bracketF = [o for o in offerings if o['bracket'] == 'F']
        bracketA = [o for o in offerings if o['bracket'] != 'F']
        for offerings, zRelToSurf in [(bracketF, True), (bracketA, False)]:
            if not offerings:
                continue
            bracketStr = 'F' if zRelToSurf else 'A'
            print ' *** extracting for offerings: {0:s} bracket *** '.format(bracketStr)

            for o in offerings:
                print tuple(o[k] for k in ['location', 'msldepth', 'bracket', 'variable'])
            stationNames = [o['location'] for o in offerings]
            # omit stations that are missing from sta file
            stationNames = list(set(staReader.getStations()).intersection(
                set(stationNames)))
            offerings = [o for o in offerings if o['location'] in stationNames]
            # station of each offering (may contain duplicates)
            stationNames = [o['location'] for o in offerings]
            x = np.array([staReader.getX(o['location']) for o in offerings])
            y = np.array([staReader.getY(o['location']) for o in offerings])
            zSign = 1 if zRelToSurf else -1  # zRelToSurf => depth below surface
            z = np.array([zSign * float(o['msldepth']) / 100.0
                          for o in offerings])

            # execute
            ee.setStations(stationNames, x, y, z, zRelativeToSurf=zRelToSurf)
            try:
                dcs = ee.extractDates(startTime, endTime)
                print ' * extracted'
            except Exception as e:
                print ' * extraction failed'
                print e
                dcs = []
            for dc in dcs:
                print ' '.join([dc.getMetaData('location'), dc.getMetaData('bracket'),
                                dc.getMetaData('msldepth'), dc.getMetaData('variable')])
            dcList.extend(dcs)

    return dcList

#-------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------


class selfeExtractor(object):
    """A wrapper object for calling extract_mod module. extract_mod module should not be visible outside this object."""

    def __init__(
            self,
            dataDir,
            fileName,
            transect=False,
            slab=False,
            track=False):
        """Create new instance"""
        # this shallow copy is simply a coding convenience
        self._mod = extract_mod.extract_mod
        # specify a bad value to be used in the extracted data
        self._mod.fill_in = np.nan

        self.fileName = fileName
        self.dataDir = dataDir
        self.transect = transect
        self.slab = slab
        self.track = track
        # specify coordinate system of binary model output
        self._mod.ics = 1  # 1 for cartesian coordinates, 2 for latlon in binary output
        self._mod.itransect = int(transect)
        self.startTime = None
        self.headerIsRead = False
        self.hasCoordinates = False

    def setTransectMode(self):
        self.transect = True
        self._mod.itransect = 1

    def setSlabMode(self):
        self.slab = True

    def setTrackMode(self):
        self.track = True
        self._mod.ixy_or_xyz = 2  # x,y,z depend on time (not only x,y)

    def readHeader(self, stack):
        # initialize
        fname = os.path.join(self.dataDir, '%d_%s' % (stack, self.fileName))
        if not os.path.isfile(fname):
            fnames = []
            pattern = os.path.join(self.dataDir, '?_%s' % (self.fileName))
            fnames.extend(glob.glob(pattern))
            pattern = os.path.join(self.dataDir, '??_%s' % (self.fileName))
            fnames.extend(glob.glob(pattern))
            pattern = os.path.join(self.dataDir, '???_%s' % (self.fileName))
            fnames.extend(glob.glob(pattern))
            if len(fnames) == 0:
                raise IOError('File not found: ' + fname)
            fname = fnames[0]
        rtn = extract_mod.readheader(fname)
        if rtn != 0:
            raise IOError('Error reading header: ' + fname)
        # convert start time specified in header of output file to epoch time
        self.startTime = datetime.datetime.strptime(
            ''.join(self._mod.start_time[: 23]).strip(),
            '%m/%d/%Y %H:%M:%S %Z')
        self.headerIsRead = True

    def getConnectivityArray(self):
        """Return connectivity array of the 2D horizontal grid"""
        if not self.headerIsRead:
            raise Exception(
                'file header must be read before calling this function')
        return self._mod.nm[:, :3].copy() - 1

    def getMeshNodeCoords(self):
        """Returns nodes of the 2D horizontal grid"""
        if not self.headerIsRead:
            raise Exception(
                'file header must be read before calling this function')
        return self._mod.x.copy(), self._mod.y.copy()

    def getNumberOfComponents(self):
        """Returns number of components of the field (scalar 1, vector 2,3)"""
        if not self.headerIsRead:
            raise Exception(
                'file header must be read before calling this function')
        return self._mod.ivs

    def getStartTime(self):
        """Returns simulation start time in datetime format"""
        if not self.headerIsRead:
            raise Exception(
                'file header must be read before calling this function')
        return self.startTime

        """Returns array of wet(=1) / dry(=0) elements"""

    def allocateArrays(self):
        # preallocate storage for extracted data in the fortran module
        self._mod.outtime = np.ones(self._mod.nrec) * np.nan
        if self.slab:
            self._mod.varout3 = np.ones(
                [3, self._mod.np, self._mod.nrec],
                order='F') * np.nan
        elif self.track:
            # for track ( (dim1,dim2,z),nvrt,nxy )
            # ixy_or_xyz=1 extract all z; ixy_or_xyz=2 extract given z
            self._mod.varout2 = np.ones(
                [3, self._mod.nvrt, self._mod.nxy],
                order='F') * np.nan
        else:  # station or transect
            if self._mod.itransect == 1 and self._mod.i23d == 3:
                self._mod.varout = np.ones(
                    [3, self._mod.nvrt, self._mod.nxy, self._mod.nrec],
                    order='F') * np.nan
            else:
                # allocate a degenerate dimension for depth to be consistent
                # with the profile reader
                self._mod.varout = np.ones(
                    [3, 1, self._mod.nxy, self._mod.nrec],
                    order='F') * np.nan

    def setSlab(self, z=None, k=None, zRelativeToSurf=False):
        '''assigns parameters for extracting slabs. If k is given, extracts along given S coordinate level. If z is given, vertical coordinates are interpolated to z. If zRelativeToSurf is True, z>0 is depth below surface, otherwise it is the z coordinate (z<0 under datum).'''
        if not self.slab:
            raise Exception('this function is only available for slab mode')
        self._mod.ialong_s = int(bool(z is None))
        if self._mod.ialong_s:
            self._mod.klev0 = k  # S-level to extract ( only for ialong_s=1 )
        else:
            self._mod.zout = z  # z-coord to extract
        self._mod.ifs = int(zRelativeToSurf)
        self.hasCoordinates = True

    def setTrack(self, x, y, z, t, zRelativeToSurf=False):
        '''assigns x,y, and z locations for data extraction in the fortran global variables.
        All the arrays must be of same length.
        '''
        if not self.track:
            raise Exception('this function is only available for track mode')
        nxy = len(x)
        if len(y) != nxy or len(z) != nxy or len(t) != nxy:
            print x.shape, y.shape, z.shape, t.shape
            raise Exception('all input arrays must be of the same length')
        self._mod.x00 = x
        self._mod.y00 = y
        self._mod.z00 = z
        self._mod.t00 = t
        self._mod.nxy = len(x)
        self._mod.ifs = int(zRelativeToSurf)
        self.hasCoordinates = True

    def setCoordinates(self, x, y, z=None, zRelativeToSurf=False):
        '''assigns x,y, and z locations for data extraction in the fortran global variables. If no z values are provided, set global variable specifying profile/transect output.
        Automatically checks if the x,y coordinates can be found in the grid, and removes the outliers (silently).'''
        if self.slab or self.track:
            raise Exception(
                'this function is not available for slab or track modes')
        self._mod.nxy = len(x)
        self._mod.x00 = x
        self._mod.y00 = y
        # specify coordinate system of requested z coordinates
        # 0 (input z are z coordinates) or 1 (input z are relative to free
        # surface;
        self._mod.ifs = int(zRelativeToSurf)
        if self._mod.itransect == 0:  # 1 for depth profile at specified point, 0 for timeseries at single depth
            self._mod.z00 = z
        self.hasCoordinates = True
        goodCoordinates = self.checkCoordinates()
        # find_parents allocates and fills some arrays, redo to get correct
        if not np.all(goodCoordinates == False):
            self.discardBadCoordinates()
            self.checkCoordinates()
        return goodCoordinates

    def checkCoordinates(self):
        """Check if all the horizontal coordinates are in the grid"""
        # find parent element for each location at which data will be extracted
        if not self.headerIsRead:
            self.readHeader(1)
        rtn = extract_mod.find_parents()
        self.goodCoordinates = self._mod.stations.copy().astype(bool)
        return self.goodCoordinates

    def discardBadCoordinates(self):
        self._mod.x00 = self._mod.x00[self.goodCoordinates].copy()
        self._mod.y00 = self._mod.y00[self.goodCoordinates].copy()
        self._mod.nxy = len(self._mod.x00)
        if self._mod.itransect == 0:
            self._mod.z00 = self._mod.z00[self.goodCoordinates].copy()

    def getExtractedStationZCoords(self):
        """For each station, returns the z coord where the data was extracted.
        May differ from the requested z if under bottom or over free surface."""
        actualZ = self._mod.stazcoord.copy()
        for i, z in enumerate(self._mod.z00):
            if actualZ[i] == 99:
                # z has not been modified in extraction routine
                actualZ[i] = z
        return actualZ

    def extractTrack(self):
        """Extracts track (x,y,z,t) from the data directory"""
        if not self.headerIsRead:
            self.readHeader(1)
        self.allocateArrays()
        t0 = timeMod.clock()
        sys.stdout.write('extracting track... ')
        sys.stdout.flush()
        rtn = extract_mod.find_parents()
        if rtn != 0:
            raise IOError(
                'Error finding parents ' +
                self.dataDir +
                ' ' +
                self.fileName)
        baseDir = self.dataDir
        if baseDir[-1] != '/':
            baseDir += '/'
        rtn = extract_mod.find_xyzt(baseDir, self.fileName)
        if rtn != 0:
            raise IOError(
                'Error reading data from ' +
                self.dataDir +
                ' ' +
                self.fileName)
        sys.stdout.write(' duration %.2f s\n' % (timeMod.clock() - t0))
        # varout2(1:ivs,1,it)
        nComponents = self.getNumberOfComponents()
        return self._mod.varout2[:nComponents, :1, :].copy()

    def extract(self, stacks):
        """Exctracts time series from the given stacks"""
        if self.track:
            raise Exception('this method is not callable in track mode')
        # read first header to get setup right
        # TODO handle case where stacks[0] does not exist
        self.readHeader(stacks[0])
        if not self.hasCoordinates:
            raise Exception('coordinates or slab not set')
        self.allocateArrays()
        # loop through time period of interest
        time = []
        data = []
        if not self.slab and np.all(self.goodCoordinates == False):
            # if no coordinates are found in the grid
            return time, data
        # check if all stacks exist
        verifiedFiles = []
        for i in stacks:
            fname = os.path.join(self.dataDir, '%d_%s' % (i, self.fileName))
            if not os.path.isfile(fname):
                print 'File not found, skipping: ' + fname
                #raise IOError( 'File not found: '+fname )
            else:
                verifiedFiles.append(fname)
        if not verifiedFiles:
            raise IOError(
                'No output files found ',
                stacks,
                self.dataDir,
                self.fileName)
        t0 = timeMod.clock()
        sys.stdout.write('extracting ... ')
        sys.stdout.flush()
        for fname in verifiedFiles:
            # read data from file
            #sys.stdout.write( ' ... '+fname+'\n' )
            # sys.stdout.flush()
            if self.slab:
                rtn = extract_mod.readslab(fname)
            else:
                rtn = extract_mod.readdata(fname)
            if rtn != 0:
                raise IOError('Error reading data from ' + fname)

            nComponents = self.getNumberOfComponents()
            t = self._mod.outtime[:self._mod.nrec].copy()
            if self.slab:
                # For slabs
                # varout3(1:ivs,1:nxy,1:nrec)
                d = self._mod.varout3[:nComponents, :, :].copy()
            elif self.transect:
                # For transect input (itransect=1; for 3D variables only),
                # varout(1:ivs,1:nvrt,1:nxy,1:nrec)
                # is the final output with times given by outtime(1:nrec) (in sec).
                # The vertical structure (i.e. z coordinates) is given by
                # varout(3,1:nvrt,1:nxy,1:nrec), where nvrt is the total # of
                # vertical levels.
                d = self._mod.varout[:, :, :, :].copy()
            else:  # station/profile
                # For time series input (itransect=0),
                # varout(1:ivs,1:1,1:nxy,1:nrec) is the output with
                # times given by outtime(1:nrec) (in sec), where nrec is # of time steps within
                # each stack, nxy is the # of points in bp_file, ivs=1 indicates
                # scalar (1) output; ivs=2 indicates vector outputs (e.g. 1=u; 2=v).
                # TODO handle >1 components
                d = self._mod.varout[0, 0, :, :].copy()
            # HACK workaround for f22 time arrays, force increasing time
            if len(time) > 0 and time[-1][-1] > t[0]:
                t += time[-1][-1]
            time.append(t)
            data.append(d)
        sys.stdout.write(' duration %.2f s\n' % (timeMod.clock() - t0))

        # merge arrays
        if self.slab:
            time = np.concatenate(tuple(time), axis=2)
            data = np.concatenate(tuple(data), axis=2)
        elif self.transect:
            time = np.concatenate(tuple(time), axis=3)
            data = np.concatenate(tuple(data), axis=3)
        else:
            # TODO handle >1 components
            time = np.concatenate(tuple(time), axis=1)
            data = np.concatenate(tuple(data), axis=1)
        return time, data


class extractBase(object):
    """Abstract base class for all extraction modes: station, profile, transect, slab, track."""

    def __init__(self, dataDir, fieldName, firstStack=1, modelCoordSys='spcs'):
        self.fieldName = fieldName
        self.dataDir = dataDir
        self.modelCoordSys = modelCoordSys
        self.firstStack = firstStack
        self.extractor = selfeExtractor(
            dataDir, fieldNameToFilename[fieldName])

    def extract(self, stacks):
        """Extracts data for the given stacks (list of integers). Returns [t,d],
        where t is a time array (nTime,) and d is data array (shape can vary). """
        raise NotImplementedError(
            'This method must be implemented in the derived class')

    def extractDates(
            self,
            startTime,
            endTime,
            simulationStartTime=None,
            headerStack=None):
        """Extracts data for given date range. startTime and endTime are datetime objects.

        The routine will try to deduce the simulation start date by reading the header of stack 1. To override this behavior, the user can provide the correct simulationStartTime object or another headerStack number.

        NOTE: startTime,endTime is converted to midnight of the given date, and
        the days in between are considered."""

        # deduce correct stacks
        if headerStack is None:
            headerStack = self.firstStack
        if simulationStartTime is None:
            if not self.extractor.headerIsRead:
                self.extractor.readHeader(headerStack)
            simulationStartTime = self.extractor.getStartTime()
        startStack = (startTime - simulationStartTime).days + 1
        endStack = (endTime - simulationStartTime).days
        endStack = max(startStack, endStack)
        if startStack < 0:
            print startStack, simulationStartTime
            raise Exception(
                'Negative start day: requested extraction date earlier than simulation start date.')
        stacks = range(startStack, endStack + 1)
        return self.extract(stacks)


class extractStation(extractBase):
    """Higher level extraction object for stations and profiles"""

    def __init__(
            self,
            dataDir,
            fieldName,
            firstStack=1,
            profile=False,
            modelCoordSys='spcs'):
        extractBase.__init__(
            self,
            dataDir,
            fieldName,
            firstStack,
            modelCoordSys)
        self.profile = profile
        if profile:
            self.extractor.setTransectMode()

    def setStations(
            self,
            stationNames,
            staX,
            staY,
            staZ=None,
            zRelativeToSurf=False):
        '''Sets station coordinates. Removes stations outside the grid.'''
        if len(stationNames) != len(staX):
            raise Exception(
                'stationNames and station coordinates do not match')
        self.zRelativeToSurf = zRelativeToSurf
        # This automatically discards outside coordinates in the extractor
        stationsInGrid = self.setCoordinates(staX, staY, staZ, zRelativeToSurf)
        badStations = [s for i, s in enumerate(
            stationNames) if not stationsInGrid[i]]
        for station in badStations:
            print 'Station out of grid, removing: ' + station
        self.stationNames = [s for s in stationNames if s not in badStations]
        if not self.stationNames:
            print 'Warning: no stations found in grid.'
        mask = stationsInGrid.copy().astype(bool)
        self.staX = staX[mask]  # always in spcs coords
        self.staY = staY[mask]
        if self.profile:
            self.staZ = None
        else:
            self.staZ = staZ[mask]

    def setCoordinates(self, staX, staY, staZ, zRelativeToSurf):
        x = staX.copy()
        y = staY.copy()
        if self.modelCoordSys == 'utm':
            # convert stations to utm for extraction
            for i in range(len(staX)):
                x[i], y[i] = spcs2utm(x[i], y[i])
        return self.extractor.setCoordinates(x, y, staZ, zRelativeToSurf)

    def extract(self, stacks):
        """Extract data for given stacks. Returns a list of dataContainers, one for each station."""
        t, d = self.extractor.extract(stacks)
        if t == []:
            return []
        ta = timeArray.timeArray(
            t, 'simulation', self.extractor.startTime).asEpoch()
        var = self.fieldName
        dcList = []
        if not self.profile:
            extractedZ = self.extractor.getExtractedStationZCoords()
        for i, station in enumerate(self.stationNames):
            # Reshaping is determined by data type (scalar vs. vector & time
            # series vs. transect)
            nComponents = self.extractor.getNumberOfComponents()
            if not self.profile:
                goodIx = np.isfinite(t) * np.isfinite(d[i, :])
                ti = t[goodIx]
                di = d[i, goodIx]
                if len(ti) == 0:
                    print 'all bad data', station, len(t)
                    continue
                datai = np.reshape(di, (1, nComponents, len(ti)))
                x = self.staX[i]
                y = self.staY[i]
                z = self.staZ[i]
                zSign = 1 if self.zRelativeToSurf else -1  # zRelToSurf => depth below surface
                msldepth = str(int(round(zSign * self.staZ[i] * 100)))
                if abs(z - extractedZ[i]) > 1e-6:
                    print 'warning: station z coordinate has changed: ', station, z, extractedZ[i], extractedZ[i] - z
            else:
                z = d[2, :, i, :]  # (dim1, dim2, z), nvrt, nxy, ntime
                # remove bad values (below bottom)
                goodIxZ = np.logical_not(np.isnan(np.sum(z, axis=1)))
                if np.any(goodIxZ == False):
                    # check if nan mask is changing in time
                    maxNaNIx = 0
                    minNaNIx = 1e9
                    for kkk in range(z.shape[1]):
                        nanIx = np.nonzero(np.isnan(z[:, kkk]))[0]
                        maxNaNIx = max(nanIx.max(), maxNaNIx)
                        minNaNIx = min(nanIx.max(), minNaNIx)
                    if minNaNIx != maxNaNIx:
                        print 'zCoord nan mask changes in time:', minNaNIx, maxNaNIx

                if np.all(goodIxZ == False):
                    continue
                z = z[goodIxZ, :]
                # x,y static in time
                x = self.staX[i] * np.ones((z.shape[0],))
                y = self.staY[i] * np.ones((z.shape[0],))

                ti = t
                # from (dim,vrt,xy,time) to (dim,xyz,time)
                datai = np.reshape(
                    d[: nComponents, goodIxZ, i, :],
                    (nComponents, len(z),
                     len(ti)))
                # from (dim,xyz,time) to (xyz,dim,time)
                datai = datai.swapaxes(0, 1)
                msldepth = 'prof'

            tai = timeArray.timeArray(
                ti, 'simulation', self.extractor.startTime).asEpoch()
            # if suspected bad values, print warning
            hasBadValues = np.isnan(datai).any() or np.isinf(
                datai).any() or np.any(datai < VALID_MIN)
            if hasBadValues:
                print 'Warning: bad values in', station, msldepth
            meta = {}
            meta['location'] = station
            meta['instrument'] = 'model'
            meta['variable'] = var
            # TODO check bracket for profile ??
            if self.profile:
                meta['bracket'] = 'A'
                meta['dataType'] = 'profile'
            else:
                meta['bracket'] = 'F' if self.zRelativeToSurf else 'A'
                meta['msldepth'] = msldepth
                meta['dataType'] = 'timeseries'
            dc = dataContainer.dataContainer(
                '', tai, x, y, z, datai, fieldNameList[var],
                coordSys='spcs', metaData=meta)

            dcList.append(dc)
        return dcList

#-------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------


def parseCommandLine():
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option(
        '-r',
        '--runTag',
        action='store',
        type='string',
        dest='runTag',
        help='Run tag, used as a label in post-proc.')
    parser.add_option(
        '-d',
        '--dataDirectory',
        action='store',
        type='string',
        dest='dataDir',
        help='directory where model outputs are stored')
    parser.add_option(
        '-C',
        '--read-netcdf',
        action='store_true',
        dest='readNetcdf',
        help='Extract from SELFE netcdf output files instead of SELFE binary files (default %default)',
        default=False)
    parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
    parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
    parser.add_option('', '--stacks', action='store', type='string',
                      dest='stackStr', help='range of output files to read '
                      '(e.g 1,14) if start,end not given')
    parser.add_option(
        '-v',
        '--variable',
        action='store',
        type='string',
        dest='varList',
        help='variable(s) to extract: elev,temp,salt, ...\nTo use specific output file define extrension, e.g. salt.70')
    parser.add_option(
        '-t',
        '--stationFile',
        action='store',
        type='string',
        dest='stationFile',
        help='text file (*.cvs) containing stations '
        'and horizontal coordinates (in spcs coordinates)',
        default=None)
    parser.add_option('-c', '--modelCoordSys', action='store', type='string',
                      dest='modelCoordSys', default='spcs',
                      help='horizontal coordinate system used in model: '
                      'spcs or utm (Default: %default)')
    parser.add_option(
        '-o',
        '--outDirectory',
        action='store',
        type='string',
        dest='outDir',
        help='base directory for netCDF file tree '
        '(optional)')
    parser.add_option(
        '-p',
        '--profile',
        action='store_true',
        dest='profile',
        help='extract vertical profile instead of value at given z level (default %default)',
        default=False)
    parser.add_option(
        '-n',
        '--no-offerings',
        action='store_true',
        dest='noOfferings',
        help='Do not extract based on offering strings. If set, stationFile must contain (x,y,z) coordinates (default %default)',
        default=False)
    parser.add_option(
        '',
        '--save-in-tree',
        action='store_true',
        dest='saveInTree',
        help='saves extracted data in file tree with monthly files instead of a single file (default %default)',
        default=False)
    parser.add_option(
        '-A',
        '--all-stations',
        action='store_true',
        dest='allOfferings',
        help='Do not extract based on available offerings, but for all stations and depths in the database (default %default)',
        default=False)
    parser.add_option(
        '-T', '--tracerModel', action='store', type='string',
        dest='tracerModel',
        help='Enable extraction of tracers: sed, oxy, generic. Must '
        'supply number of tracers for \'sed\' and \'generic\' '
        'models via the -N switch. \'oxy\' model provides tracers: '
        '\'NO3\',\'NH4\',\'phy\',\'zoo\',\'det\' and \'oxy\'.', default=None)
    parser.add_option(
        '-N',
        '--numTracers',
        action='store',
        type='int',
        dest='numTracers',
        help='Tracer number to extract for \'sed\' and \'generic\' models',
        default=None)
    parser.add_option(
        '',
        '--decimals',
        action='store',
        type='int',
        dest='digits',
        help='Round extracted data to given decimal precision to save disk space',
        default=None)

    (options, args) = parser.parse_args()

    dataDir = options.dataDir
    varList = options.varList.split(',') if options.varList else None
    stationFile = options.stationFile
    modelCoordSys = options.modelCoordSys
    outDir = options.outDir
    startStr = options.startStr
    endStr = options.endStr
    stackStr = options.stackStr
    readNetcdf = options.readNetcdf
    profile = options.profile
    noOfferings = options.noOfferings
    allOfferings = options.allOfferings
    runTag = options.runTag
    saveInTree = options.saveInTree
    tracerModel = options.tracerModel
    numTracers = options.numTracers
    digits = options.digits

    if not dataDir:
        parser.print_help()
        parser.error('dataDir  undefined')
    if not varList and tracerModel is None:
        parser.print_help()
        parser.error('variable undefined')
    # if not outDir :
        # parser.print_help()
        #parser.error('outDir   undefined')
    if startStr is None and stackStr is None:
        parser.print_help()
        parser.error('startStr undefined')
    if endStr is None and stackStr is None:
        parser.print_help()
        parser.error('endStr   undefined')
    if noOfferings and not stationFile:
        parser.print_help()
        parser.error(
            'stationFile must be provided, if offerings are not used.')
    if not runTag:
        parser.print_help()
        parser.error('runTag  undefined')
    if tracerModel:
        if not numTracers and tracerModel.split('.')[0] in ['sed', 'generic']:
            parser.print_help()
            parser.error(
                'numTracers must be provided if sed or generic tracer models are used.')
        addTracers(tracerModel, varList, numTracers=numTracers)

    if stackStr is not None:
        limits = [int(v) for v in stackStr.split(',')]
        stacks = np.arange(limits[0], limits[1] + 1)
    else:
        stacks = None

    if startStr is not None:
        startTime = datetime.datetime.strptime(startStr, '%Y-%m-%d')
        endTime = datetime.datetime.strptime(endStr, '%Y-%m-%d')
    else:
        startTime = endTime = None

    print 'Parsed options:'
    if stackStr is None:
        print ' - time range:', str(startTime), '->', str(endTime)
    else:
        print ' - stacks:', str(stacks[0]), '->', str(stacks[-1])
    print ' - dataDir:', dataDir
    print ' - SELFE output format:', 'netCDF' if readNetcdf else 'binary'
    print ' - runTag:', runTag
    if outDir:
        print ' - output dir:', outDir
    print ' - variables:', varList
    if profile:
        print ' - extracting vertical profiles'
    if noOfferings and stationFile is not None:
        print ' - reading stations from:', stationFile
    sys.stdout.flush()

    for var in varList:
        if noOfferings:
            if readNetcdf:
                from crane.data.ncExtract import extractForStations as extractNetCDF
                dcs = extractNetCDF(
                    dataDir,
                    var,
                    stationFile,
                    startTime,
                    endTime,
                    profile,
                    stacks=stacks)
            else:
                dcs = extractForXYZ(
                    dataDir,
                    var,
                    stationFile,
                    startTime,
                    endTime,
                    profile,
                    modelCoordSys)
        else:
            # get all available offerings
            import crane.data.netcdfCacheInterface as netcdfDB
            if allOfferings:
                print(
                    'fetching {0:s} offerings from the database for all stations'.format(var))
                offerings = netcdfDB.getAllOfferings([var.split('.')[0]])
            else:
                print(
                    'fetching {0:s} offerings from the database for the time period'.format(var))
                offerings = netcdfDB.getAvailableOfferings(
                    startTime, endTime, [var.split('.')[0]])
            if len(offerings) == 0:
                print(
                    'No offerings received, skipping variable {0:s}'.format(var))
            if readNetcdf:
                from crane.data.ncExtract import extractForOfferings as extractNetCDF
                dcs = extractNetCDF(
                    dataDir,
                    var,
                    offerings,
                    startTime,
                    endTime,
                    profile,
                    stationFile)
            else:
                dcs = extractForOfferings(
                    dataDir,
                    var,
                    offerings,
                    startTime,
                    endTime,
                    profile,
                    modelCoordSys,
                    stationFile)
        for dc in dcs:
            dc.setMetaData('tag', runTag)
        import crane.data.dirTreeManager as dtm
        if saveInTree:
            rule = 'monthlyFile'
        else:
            rule = 'singleFile'
        dtm.saveDataContainerInTree(
            dcs,
            rootPath=outDir,
            rule=rule,
            dtype=np.float32,
            overwrite=True,
            compress=True,
            digits=digits)

if __name__ == '__main__':
    parseCommandLine()
