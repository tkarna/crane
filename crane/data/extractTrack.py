#!/usr/bin/python
"""
Extract track (xyzt) data using the efficient SELFE extract_mod python module.

Examples:


Tuomas Karna 2012-11-29
"""

import numpy as np
import time as timeMod
import os
import sys
import datetime
import subprocess as sub

from crane.data import dataContainer
from crane.data import timeArray
from crane.data import extractStation
from crane.files import buildPoints
from crane.physicalVariableDefs import addTracers

#-------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------


def extractForDataContainer(
        dataDir,
        trackDC,
        var=None,
        name=None,
        bpFile=None):
    """Extracts a track based on x,y,z,time in the given dataContaner."""
    # sanity check that trackDC is suitable
    if not trackDC.isTrack():
        raise Exception(
            'given dataContainer does not contain track information')
    # track metadata
    bracket = trackDC.getMetaData('bracket')
    if not name:
        name = trackDC.getMetaData('location')
    if not var:
        var = trackDC.getMetaData('variable')
    x = trackDC.x.flatten()
    y = trackDC.y.flatten()
    z = trackDC.z.flatten()

    # Overwrite x,y coordinate if using alternative coordinate system (open
    # channels)
    if bpFile is not None:
        print 'Overwriting xy from buildpoint'
        bp = buildPoints.BuildPoint()
        bp.readFileFromDisk(bpFile)
        x = bp.getX()
        y = bp.getY()

    # expand scalar coordinates to time-dep track
    nx = max(max(len(x), len(y)), len(z))
    if len(x) == 1:
        x = np.tile(x, (nx,))
    if len(y) == 1:
        y = np.tile(y, (nx,))
    if len(z) == 1:
        z = np.tile(z, (nx,))
    ta = trackDC.time
    zRelativeToSurf = True if bracket == 'F' else False
    ee = extractTrack(
        dataDir,
        var,
        firstStack=1,
        modelCoordSys=trackDC.coordSys)
    print ' * extracting track', trackDC.getMetaData('location'), trackDC.getMetaData('variable')
    ee.setTrack(name, x, y, z, ta, zRelativeToSurf)
    dc = ee.extract()
    return dc


#------------------------------------------------------------------------------
# Classes
#------------------------------------------------------------------------------

class extractTrack(extractStation.extractBase):
    """A higher lever extraction object for tracks"""

    def __init__(self, dataDir, fieldName, firstStack=1, modelCoordSys='spcs'):
        super(
            extractTrack,
            self).__init__(
            dataDir,
            fieldName,
            firstStack,
            modelCoordSys)
        self.extractor.setTrackMode()
        self.name = None
        self.extractor.readHeader(firstStack)

    def setTrack(self, name, x, y, z, time, zRelativeToSurf=False):
        '''Set track coordinates parameters. Time is given as a timeArray object.'''
        self.name = name
        self.zRelativeToSurf = zRelativeToSurf
        # convert time to simulation time ( seconds since the start )
        startTime = self.extractor.getStartTime()
        t = time.asSimulation(startTime).array
        self.extractor.setTrack(x, y, z, t, zRelativeToSurf)
        self.x = x
        self.y = y
        self.z = z
        self.time = time

    def extract(self):
        """Extract data from data directory. Returns the data in a dataContainer."""
        data = self.extractor.extractTrack()  # (dim,1,nTime)
        data = data.swapaxes(0, 1)  # (1,dim,nTime)
        x = self.x[None, :]
        y = self.y[None, :]
        z = self.z[None, :]
        var = self.fieldName

        # TODO remove bad values (nan mask?)
        data[data < extractStation.VALID_MIN] = np.nan
        # TODO add support for removing/keeping dry elements

        # if suspected bad values, print warning
        hasBadValues = np.isnan(data).any() or np.isinf(
            data).any() or np.any(data < extractStation.VALID_MIN)
        if hasBadValues:
            print 'Warning: bad values in', self.name
        meta = {}
        meta['dataType'] = 'track'
        meta['location'] = self.name
        meta['instrument'] = 'model'
        meta['bracket'] = 'F' if self.zRelativeToSurf else 'A'
        meta['variable'] = var
        dc = dataContainer.dataContainer(
            '', self.time, x, y, z, data, fieldNameList[var],
            coordSys='spcs', metaData=meta, acceptNaNs=True)
        return dc

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
    parser.add_option(
        '-v',
        '--variable',
        action='store',
        type='string',
        dest='varList',
        help='variable(s) to extract: elev,temp,salt, ...\nTo use specific output file define extrension, e.g. salt.70')
    parser.add_option('-n', '--name', action='store', type='string',
                      dest='name', help='name of the track for identification')
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
    # parser.add_option('-s', '--start', action='store', type='string',
    # dest='startStr', help='Date to start processing')
    # parser.add_option('-e', '--end', action='store', type='string',
    # dest='endStr', help='Date to end processing')
    parser.add_option('', '--stacks', action='store', type='string',
                      dest='stackStr', help='range of output files to read '
                      '(e.g 1,14) if start,end not given')
    parser.add_option(
        '-i',
        '--inputTrack',
        action='store',
        type='string',
        dest='ncFile',
        help='netCDF file (dataContainer) of the track to extract, e.g. auv mission')
    parser.add_option(
        '',
        '--decimals',
        action='store',
        type='int',
        dest='digits',
        help='Round extracted data to given decimal precision to save disk space',
        default=None)
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
        help='Number of tracers to support for \'sed\' and \'generic\' models',
        default=None)
    parser.add_option(
        '-b',
        '--buildPoint',
        action='store',
        type='str',
        dest='bpFile',
        help='BuildPoint file to read x,y location to extract track from model',
        default=None)

    (options, args) = parser.parse_args()

    dataDir = options.dataDir
    varList = options.varList.split(',') if options.varList else None
    name = options.name
    modelCoordSys = options.modelCoordSys
    outDir = options.outDir
    #startStr      = options.startStr
    #endStr        = options.endStr
    stackStr = options.stackStr
    readNetcdf = options.readNetcdf
    runTag = options.runTag
    ncFile = options.ncFile
    runTag = options.runTag
    digits = options.digits
    tracerModel = options.tracerModel
    numTracers = options.numTracers
    bpFile = options.bpFile

    if not dataDir:
        parser.print_help()
        parser.error('dataDir  undefined')
    if not varList and tracerModel is None:
        parser.print_help()
        parser.error('variable undefined')
    # if not outDir :
        # parser.print_help()
        #parser.error('outDir   undefined')
    # if not name :
        # parser.print_help()
        #parser.error('track name undefined')
    if not ncFile:
        parser.print_help()
        parser.error('input netcdf file undefined')
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

    print 'Parsed options:'
    if name:
        print ' - name ', name
    print ' - input track', ncFile
    if stackStr is not None:
        print ' - stacks:', str(stacks[0]), '->', str(stacks[-1])
    # if startTime :
        # print ' - time range:',str(startTime),'->', str(endTime)
    print ' - dataDir', dataDir
    print ' - SELFE output format:', 'netCDF' if readNetcdf else 'binary'
    print ' - output dir', outDir
    print ' - variables ', varList
    print ' - model coord system', modelCoordSys
    print ' - runTag ', runTag
    if bpFile:
        print ' - bpFile ', bpFile
    sys.stdout.flush()

    import crane.data.dirTreeManager as dtm
    rule = 'singleFile'
    if readNetcdf:
        from crane.data.ncExtract import extractTrackForDataContainer as extractNetCDF

    if tracerModel == 'sed':
        trackDC = dataContainer.dataContainer.loadFromNetCDF(ncFile)
        if readNetcdf:
            dc = extractNetCDF(dataDir, trackDC, 'trcr_1', name, stacks=stacks)
        else:
            dc = extractForDataContainer(
                dataDir, trackDC, 'sed_1', name, bpFile)
        dc.setMetaData('tag', runTag)
        dc.setMetaData('variable', 'sed_1')
        dc.fieldNames = ['sed_1']
        dtm.saveDataContainerInTree(
            dc,
            rootPath=outDir,
            rule=rule,
            dtype=np.float32,
            overwrite=True)
        # Combines sediment files into one
        if numTracers > 1:
            for sed_class in range(2, numTracers + 1):
                if readNetcdf:
                    tmp = extractNetCDF(
                        dataDir, trackDC, 'trcr_%d' %
                        sed_class, name, bpFile, stacks=stacks)
                else:
                    tmp = extractForDataContainer(
                        dataDir, trackDC, 'sed_%d' %
                        sed_class, name, bpFile)
                dc.data = dc.data + tmp.data
                tmp.setMetaData('tag', runTag)
                tmp.setMetaData('variable', 'sed_%d' % sed_class)
                tmp.fieldNames = ['sed_%d' % sed_class]
                dtm.saveDataContainerInTree(
                    tmp,
                    rootPath=outDir,
                    rule=rule,
                    dtype=np.float32,
                    overwrite=True)
        dc.setMetaData('tag', runTag)
        dc.setMetaData('variable', 'sed')
        dc.fieldNames = ['sed']
        dtm.saveDataContainerInTree(
            dc,
            rootPath=outDir,
            rule=rule,
            dtype=np.float32,
            overwrite=True)
    else:
        for var in varList:
            trackDC = dataContainer.dataContainer.loadFromNetCDF(ncFile)
            if readNetcdf:
                dc = extractNetCDF(
                    dataDir, trackDC, var, name, bpFile, stacks=stacks)
            else:
                dc = extractForDataContainer(
                    dataDir, trackDC, var, name, bpFile)
            dc.setMetaData('tag', runTag)
            dtm.saveDataContainerInTree(
                dc,
                rootPath=outDir,
                rule=rule,
                dtype=np.float32,
                overwrite=True)

if __name__ == '__main__':
    parseCommandLine()
