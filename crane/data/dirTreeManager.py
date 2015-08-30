"""
Implementation of routines for reading and writing netcdf files from/to
a directory tree.

The two high level routines users normally need to call are:
getDataContainer
saveDataContainerInTree

Currently supported file trees are:

singleFileTree, with file naming rule:
    tag/data/location_variable_msldepth_startTime-endTime.nc

monthlyFileTree, with file naming rule:
    tag/data/stations/location/offering/yyyymm.nc

Any prefix path to these file trees can be set with rootPath parameter.

Tuomas Karna 2015-08-08
"""
import os
import numpy as np
import datetime
from dateutil import rrule, relativedelta
from glob import glob
import traceback
import sys

from crane.utility import createDirectory
from crane.data import netcdfIO

ruleAlias = {'singleFile': ['singleFile', 'old', ],
             'monthlyFile': ['monthlyFile', 'default', ],
             }


def detectTreeRule(verbose=False, **kwargs):
    """
    Test if any file reader is able to find matching files,
    return first good one.
    """
    rules = [singleFileTree(verbose=verbose),
             monthlyFileTree(verbose=verbose)]
    for r in rules:
        files = r.findMatchingFiles(**kwargs)
        if len(files) > 0:
            if verbose:
                print '* directory structure detected:', r.__class__.__name__
            return r
    argsStr = ' '.join(['{0:s}={1:s}'.format(k, repr(kwargs[k])) for k in kwargs])
    raise Exception('Could not detect rule: No matching files found: '+argsStr)


def getDataContainer(rootPath=None, rule=None, dataType=None, tag=None,
                     location=None, variable=None,
                     startTime=None, endTime=None,
                     msldepth=None, slevel=None, verbose=False):
    """
    Reads a dataContainer from file tree given the meta data. If multiple
    files match the meta data, the first one is returned.

    Parameters
    ----------
    rootPath : string, optional
        Set the root path of the file tree. By default the current
        directory is used.
    rule : string or fileTree object, optional
        Specify which tree abstraction to use. Possible values are
        'singleFile' or 'monthlyFile'. If unspecified, first
        method that finds a matching file is used.
    dataType : string
        dataType of the dataContainer to read. Possible values are
        'timeseries', 'profile', 'track', 'transect', 'slab', etc.
    tag : string, optional
    location : string, optional
    variable : string, optional
    msldepth : string, optional
    slevel : string, optional
        Meta data of the data to look for.
    startTime : datetime, optional
    endTime : datetime, optional
        start and end times of the data. If the data set on disk is
        longer, than the requested range, only the matching time steps
        are read. Saves time in case of long, large data sets.
    verbose : bool, optional
        Print information on stdin

    Returns
    -------
    dc : dataContainer or meshContainer
        First data set that matches the query. In case of meshed data
        returns a meshContainer object
    """
    return getAllDataContainers(rootPath, rule, dataType, tag, location,
                                variable, startTime, endTime, msldepth,
                                slevel, verbose)[0]


def getAllDataContainers(rootPath=None, rule=None, dataType=None, tag=None,
                         location=None, variable=None,
                         startTime=None, endTime=None,
                         msldepth=None, slevel=None, verbose=False):
    """
    Reads all matching dataContainers from file tree given the meta data.

    Parameters
    ----------
    rootPath : string, optional
        Set the root path of the file tree. By default the current
        directory is used.
    rule : string or fileTree object, optional
        Specify which tree abstraction to use. Possible values are
        'singleFile' or 'monthlyFile'. If unspecified, first
        method that finds a matching file is used.
    dataType : string
        dataType of the dataContainer to read. Possible values are
        'timeseries', 'profile', 'track', 'transect', 'slab', etc.
    tag : string, optional
    location : string, optional
    variable : string, optional
    msldepth : string, optional
    slevel : string, optional
        Meta data of the data to look for.
    startTime : datetime, optional
    endTime : datetime, optional
        start and end times of the data. If the data set on disk is
        longer, than the requested range, only the matching time steps
        are read. Saves time in case of long, large data sets.
    verbose : bool, optional
        Print information on stdin

    Returns
    -------
    dcs : list of dataContainer or meshContainer
        All data sest that matche the query. In case of meshed data
        returns a meshContainer object
    """
    keys = ['dataType', 'tag', 'location', 'variable',
            'startTime', 'endTime', 'msldepth', 'slevel',
            'rootPath']
    vals = [dataType, tag, location, variable,
            startTime, endTime, msldepth, slevel, rootPath]
    kwargs = dict(zip(keys, vals))
    if rule is None:
        tree = detectTreeRule(verbose=verbose, **kwargs)
    elif isinstance(rule, str):
        if rule in ruleAlias['singleFile']:
            tree = singleFileTree(verbose=verbose)
        elif rule in ruleAlias['monthlyFile']:
            tree = monthlyFileTree(verbose=verbose)
        else:
            raise Exception('Unknown dir tree rule: '+rule)
    else:
        # assume that user provided the tree object
        tree = rule

    argsStr = ' '.join(['{0:s}={1:s}'.format(k, repr(kwargs[k])) for k in kwargs])
    dcList = tree.readFiles(**kwargs)
    if len(dcList) == 0:
        raise Exception('Reading data from tree failed: ' + argsStr)
    return dcList


def saveDataContainerInTree(dcs, rootPath=None, rule=None, dtype=np.float64,
                            overwrite=False, compress=False, digits=None):
    """Saves dataContainer(s) in a tree with the given fileTreeRule.

    Parameters
    ----------
    dcs : dataContainer or meshContainer or a list of thereof
        single object or list of objects to save.
    rootPath : string, optional
        Set the root path of the file tree. By default the current
        directory is used.
    rule : string or fileTree object, optional
        Specify which tree abstraction to use. Possible values are
        'singleFile' or 'monthlyFile'. If unspecified, first
        method that finds a matching file is used.
    dtype : numpy dtype, optional
        data type of the arrays to write to disk. If unspecified
    overwrite : bool, optional
        Overwrite existing files. (%default = False)
    compress : bool, optional
        Compress netcdf files. (%default = False)
    digits : int, optional
        If set, truncates data precision to given number of significant digits.
    """
    if not isinstance(dcs, list):
        dcs = [dcs]

    if rule is None:
        rule = 'singleFile'
    if isinstance(rule, str):
        if rule in ruleAlias['singleFile']:
            tree = singleFileTree()
        elif rule in ruleAlias['monthlyFile']:
            tree = monthlyFileTree()
        else:
            raise Exception('Unknown dir tree rule: '+rule)
    else:
        raise Exception('rule must be a string')

    for dc in dcs:
        tree.saveDataContainer(dc, overwrite=overwrite, dtype=dtype,
                               compress=compress, digits=digits,
                               rootPath=rootPath)


class fileTree(object):
    """
    A base class for all file tree objects.
    """

    def readFiles(self, *args, **kwargs):
        raise NotImplementedError('This method must be implemented in derived class')

    def saveDataContainer(self, *args, **kwargs):
        raise NotImplementedError('This method must be implemented in derived class')


class singleFileTree(fileTree):
    """A class for storing dataContainers in a single netcdf file.

    Naming rule:
    tag/data/location_variable_msldepth_startTime-endTime.nc
    tag/data/datatype/location_variable_0_startTime-endTime.nc
    e.g.
    obs/data/saturn01_salt_1950_2012-05-01_2012-05-16.nc
    run127/data/transect/ncAUVeast_salt_0_2012-05-01_2012-05-19.nc
    """
    # list of metadata keys
    metaKeys = ['startTime', 'endTime', 'dataType', 'tag', 'location',
                'variable', 'msldepth', 'slevel', 'bracket']

    def __init__(self, verbose=False):
        self.verbose = verbose

    def generateFileName(self, startTime, endTime, dataType, tag, location,
                         variable, msldepth=None, slevel=None, rootPath=None,
                         **kwargs):
        """Generates a filename for the given metadata"""
        if isinstance(startTime, datetime.datetime):
            startStr = startTime.strftime('%Y-%m-%d')
        else:  # assume str
            startStr = startTime
        if isinstance(startTime, datetime.datetime):
            endStr = endTime.strftime('%Y-%m-%d')
        else:  # assume str
            endStr = endTime
        if dataType == 'timeseries':
            pattern = '{tag:s}/data/{loc:s}_{var:s}_{dep:s}_{st:s}_{et:s}.nc'
            f = pattern.format(tag=tag, loc=location, dep=msldepth,
                               var=variable, st=startStr, et=endStr)
        else:  # any other data type
            if slevel is None or slevel == '*':
                dep = '0'
            else:
                dep = 's' + str(slevel)
            pattern = '{tag:s}/data/{typ:s}/{loc:s}_{var:s}_{dep:s}_{st:s}_{et:s}.nc'
            f = pattern.format(tag=tag, typ=dataType, loc=location,
                               dep=dep, var=variable, st=startStr, et=endStr)
        if rootPath is not None:
            f = os.path.join(rootPath, f)
        return f

    def generateSearchPattern(self, **kwargs):
        """
        Generates a pattern for seaching files.
        A wildcard '*' is used for all missing inputs.
        """
        kw = dict.fromkeys(self.metaKeys)
        # replace None with user input
        kw.update(kwargs)
        for arg in kw.keys():
            if arg not in ['rootPath'] and kw[arg] is None:
                # replace None with arg or wildcard string
                kw[arg] = '*'
        pathPattern = self.generateFileName(**kw)
        return pathPattern

    def saveDataContainer(self, dc, dtype=np.float64, overwrite=False,
                          rootPath=None, compress=False, digits=None):
        m = dc.getMetaData()
        st = dc.time.getDatetime(0)
        et = dc.time.getDatetime(-1)
        m['startTime'] = st
        m['endTime'] = et
        m['rootPath'] = rootPath
        filename = self.generateFileName(**m)
        netcdfIO.netcdfIO(filename).saveDataContainer(dc, dtype, overwrite,
                                             compress, digits)

    def findMatchingFiles(self, **kwargs):
        """Returns a list of all files that match the query."""
        st = kwargs.get('startTime')
        et = kwargs.get('endTime')
        m = dict(kwargs)
        m['startTime'] = None
        m['endTime'] = None
        pattern = self.generateSearchPattern(**m)
        if self.verbose:
            print 'searching for ' + pattern
        files = glob(pattern)
        return files

    def readFiles(self, **kwargs):
        """Tries to find files matching the given keywords. Undefined keywords are
        treated as wildcards '*' in the file search. Returns all matching
        dataContainers in a list.
        """
        rootPath = kwargs.get('rootPath')
        st = kwargs.get('startTime')
        et = kwargs.get('endTime')
        files = self.findMatchingFiles(**kwargs)
        dcList = []
        for f in files:
            if os.path.isfile(f):
                nc = netcdfIO.netcdfIO(f)
                try:
                    desc, ta = nc.readHeader(verbose=self.verbose)
                    if (st is not None and et is not None and not ta.overlaps(st, et)):
                        continue  # time out of requested range
                    dc = nc.readToDataContainer(st, et, verbose=self.verbose)
                    dcList.append(dc)
                except Exception as e:
                    print 'Error while reading file: {0:s}'.format(f)
                    print e
                    traceback.print_exc(file=sys.stdout)
        return dcList


class monthlyFileTree(fileTree):
    """A class for storing dataContainers in monthly files.

    Naming rule:
    tag/data/stations/location/offering/yyyymm.nc
    tag/dataDir/dataType/location/variable/yyyymm.nc
    e.g.
    run27/data/stations/saturn01/saturn01.0.F.var/201205.nc
    run127/data/transect/nc_fine/salt/201205.nc
     """
    # list of metadata keys
    metaKeys = ['startTime', 'endTime', 'dataType', 'tag', 'location',
                'variable', 'msldepth', 'slevel', 'bracket']

    def __init__(self, verbose=False):
        self.verbose = verbose

    def generateFileName(self, startTime, endTime, dataType, tag,
                         location, variable, msldepth=None, bracket=None,
                         slevel=None, rootPath=None, **kwargs):
        """Returns a path to a file for the given metadata."""

        if not dataType:
            raise Exception('dataType string must be given')
        # subdir for each datatype
        typeDirs = {'timeseries': 'stations',
                    'transect': 'transect',
                    'slab': 'slab',
                    'track': 'track',
                    'auv': 'auv',
                    'profile': 'profile'}
        subdir = typeDirs.get(dataType, dataType)
        if dataType == 'timeseries':
            if msldepth is None:
                raise Exception('msldepth required')
            if bracket is None:
                raise Exception('bracket required')
            offer = '.'.join([location, msldepth, bracket, variable])
            directory = os.path.join(tag, 'data', subdir, location, offer)
        else:
            if msldepth is not None and msldepth != '*':
                variable += '.'+msldepth
            elif slevel is not None and slevel != '*':
                variable += '.s'+str(slevel)
            directory = os.path.join(tag, 'data', subdir, location, variable)
        if isinstance(startTime, str):
            fname = '{0:s}.nc'.format(startTime)
        else:
            fname = '{y:d}{m:02d}.nc'.format(y=startTime.year,
                                             m=startTime.month)
        if rootPath is not None:
            directory = os.path.join(rootPath, directory)
        fullpath = os.path.join(directory, fname)
        return fullpath

    def generateSearchPattern(self, **kwargs):
        """
        Generates a pattern for seaching files.
        A wildcard '*' is used for all missing inputs.
        """
        kw = dict.fromkeys(self.metaKeys)
        # replace None with user input
        kw.update(kwargs)
        # TODO handle case where start/end months differ?
        if kw['startTime'] is None or kw['endTime'] is None:
            kw['startTime'] = None
            kw['endTime'] = None
        for arg in kw.keys():
            if arg not in ['rootPath'] and kw[arg] is None:
                # replace None with arg or wildcard string
                kw[arg] = '*'
        pathPattern = self.generateFileName(**kw)
        return pathPattern

    def getMonthWindows(self, startTime, endTime):
        """Returns an interator for time windows."""
        windows = []
        for dt in rrule.rrule(rrule.MONTHLY, until=endTime,
                              dtstart=datetime.datetime(startTime.year,startTime.month,1)):
            windows.append((dt, dt+relativedelta.relativedelta(months=1)))
        return windows

    def saveDataContainer(self, dc, dtype=np.float64, overwrite=False,
                          rootPath=None, compress=False, digits=None):
        # get metadata from dataContainer
        meta = dc.getMetaData()
        dataType = meta.get('dataType')
        tag = meta.get('tag')
        loc = meta.get('location')
        var = meta.get('variable')
        msldepth = meta.get('msldepth')
        slevel = meta.get('slevel')
        bracket = meta.get('bracket')
        # figure out range of months
        startTime = dc.time.getDatetime(0)
        endTime = dc.time.getDatetime(-1)
        windows = self.getMonthWindows(startTime, endTime)
        # for each month
        for win_start, win_end in windows:
            # generate file name
            fn = self.generateFileName(win_start, win_end, dataType,
                                       tag, loc, var,
                                       msldepth, bracket, slevel,
                                       rootPath=rootPath)
            print 'saving to', fn
            createDirectory(os.path.split(fn)[0])
            # crop data dataContainer
            try:
                cropped_dc = dc.timeWindow(win_start, win_end, includeEnd=False)
            except:
                print 'Given time window out of range'
                continue
            #  save
            netcdfIO.netcdfIO(fn).saveDataContainer(cropped_dc, dtype,
                                                    overwrite, compress,
                                                    digits)

    def getLocalFilename(self, filename, rootPath=None):
        """Removes base path from the file name"""
        localPath = filename
        if rootPath is not None and len(rootPath) > 0:
            if rootPath[-1] != '/':
                rootPath += '/'
            localPath = filename.replace(rootPath, '')
        return localPath

    def parseFilename(self, filename, rootPath=None):
        """Parses metadata from filename
        """
        msldepth = None
        slevel = None
        meta = {}
        localFile = self.getLocalFilename(filename, rootPath)
        words = localFile.split('/')
        tag = words[0]
        fname = os.path.splitext(words[-1])[0]
        dataType = words[2]
        if dataType == 'stations':
            dataType = 'timeseries'
        location = words[3]
        if dataType == 'timeseries':
            offer = words[4]
            loc, msldepth, bracket, variable = offer.split('.')
        else:
            variable = words[4]
        w = variable.split('.')
        if len(w) > 1:
            variable = w[0]
            if w[1].startswith('s'):
                slevel = int(w[1][1:])
            else:
                msldepth = w[1]
        startTime = datetime.datetime.strptime(fname, '%Y%m')

        meta['tag'] = tag
        meta['dataType'] = dataType
        meta['location'] = location
        meta['variable'] = variable
        meta['startTime'] = startTime
        meta['endTime'] = startTime + relativedelta.relativedelta(months=1)
        if msldepth is not None:
            meta['msldepth'] = msldepth
        if slevel is not None:
            meta['slevel'] = slevel
        if dataType == 'timeseries':
            meta['bracket'] = bracket
        return meta

    def getAvailableSamples(self, files, rootPath=None):
        """Traverses through files and returns metadata for all disctict data
        samples, ingnoring time stamps.

        I.e. these two files belong to the same data set:
        run27/data/stations/saturn01/saturn01.1950.A.salt/201205.nc
        run27/data/stations/saturn01/saturn01.1950.A.salt/201206.nc
        Whereas these are different data sets
        run27/data/stations/saturn01/saturn01.1950.A.salt/201205.nc
        run27/data/stations/saturn03/saturn03.1300.A.salt/201206.nc
        """
        # store values in set of tuples to omit duplicates
        sampleSet = set()
        for f in files:
            meta = self.parseFilename(f, rootPath)
            meta['startTime'] = None
            meta['endTime'] = None
            vals = tuple(meta.get(k) for k in self.metaKeys)
            sampleSet.add(vals)
        # morph back to a list of dicts
        outputs = []
        for vals in sorted(list(sampleSet)):
            outputs.append(dict(zip(self.metaKeys, vals)))
        return outputs

    def findMatchingFiles(self, **kwargs):
        """Returns a list of all files that match query."""
        st = kwargs.get('startTime')
        et = kwargs.get('endTime')
        rootPath = kwargs.get('rootPath')
        m = dict(kwargs)
        m['startTime'] = None
        m['endTime'] = None
        pattern = self.generateSearchPattern(**m)
        if self.verbose:
            print 'searching for ' + pattern
        files = sorted(glob(pattern))
        # exclude files that are out of time range
        for f in files:
            meta = self.parseFilename(f, rootPath)
            file_st = meta['startTime']
            file_et = meta['endTime']
            # skip file if its end time is less than query start time
            if (st is not None and file_et < st):
                files.remove(f)
            # skip file if its start time is greater than query end
            if (et is not None and file_st > et):
                files.remove(f)
        return files

    def readFiles(self, **kwargs):
        """Tries to find files matching the given keywords. Undefined keywords are
        treated as wildcards '*' in the file search. Returns all matching
        dataContainers in a list.
        """
        st = kwargs.get('startTime')
        et = kwargs.get('endTime')
        rootPath = kwargs.get('rootPath')
        files = self.findMatchingFiles(**kwargs)
        samples = self.getAvailableSamples(files, rootPath)
        # discard all samples that do not match the query (after parsing the filename)
        dcList = []
        for s in samples:
            # search for files that are mergeable, i.e. belong to the same data set
            s['startTime'] = st
            s['endTime'] = et
            s['rootPath'] = rootPath
            sampleFiles = self.findMatchingFiles(**s)
            # read files and merge them if they fit the query range
            outputDC = None
            for f in sampleFiles:
                if os.path.isfile(f):
                    nc = netcdfIO.netcdfIO(f)
                    try:
                        desc, ta = nc.readHeader(verbose=self.verbose)
                        file_st = ta.getDatetime(0)
                        file_et = ta.getDatetime(-1)
                        # skip file if its end time is less than query start time
                        if (st is not None and file_et < st):
                            continue
                        # skip file if its start time is greater than query end
                        if (et is not None and file_st > et):
                            continue
                        dc = nc.readToDataContainer(st, et, verbose=self.verbose)
                        if outputDC is None:
                            outputDC = dc
                        else:
                            # merge (assuming that data sets are compatible)
                            outputDC.mergeTemporal(dc, testSanity=False)
                    except Exception as e:
                        print 'Error while reading file: {0:s}'.format(f)
                        print e
                        traceback.print_exc(file=sys.stdout)
            if outputDC is not None:
                dcList.append(outputDC)
        return dcList
