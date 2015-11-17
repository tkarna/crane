#!/usr/bin/python
"""
Extract slab data using the efficient SELFE extract_mod python module.

Examples:

# extract salt,hvel, -d data dir, -t defines transect bp file, -n transect name string, -o output dir, -s -e time range
python extractSlab.py -d /home/tuomas/workspace/cmop/selfe/runs/channel_boxtest/outputs_dihv_spool10/ -v salt,hvel -n test -o tmp -s 2010-6-14 -e 2010-6-14 -z -4

Tuomas Karna 2012-11-16
"""

import numpy as np
import time as timeMod
import os
import sys
import datetime
import subprocess as sub

from crane.data import meshContainer
from crane.data import timeArray
from crane.data import extractStation
from crane.physicalVariableDefs import addTracers

#-------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------


def extractSlabForLevel(
        dataDir,
        varList,
        startTime,
        endTime,
        name,
        zCoord=None,
        kLevel=None,
        zRelativeToSurf=False,
        modelCoordSys='spcs'):

    mcs = []
    for var in varList:
        ee = extractSlab(dataDir, var, firstStack=1,
                         modelCoordSys=modelCoordSys)
        ee.setSlab(name, zCoord, kLevel, zRelativeToSurf)
        mc = ee.extractDates(startTime, endTime)
        mcs.append(mc)
    return mcs

#-------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------


class extractSlab(extractStation.extractBase):
    """A higher lever extraction object for slabs"""

    def __init__(self, dataDir, fieldName, firstStack=1, modelCoordSys='spcs'):
        super(
            extractSlab,
            self).__init__(
            dataDir,
            fieldName,
            firstStack,
            modelCoordSys)
        self.extractor.setSlabMode()
        self.name = None

    def setSlab(self, name, z=None, k=None, zRelativeToSurf=False):
        """Set slab parameters."""
        self.name = name
        self.extractor.setSlab(z, k, zRelativeToSurf)
        self.zRelativeToSurf = zRelativeToSurf
        self.alongSCoord = z is None
        self.z = z
        self.k = k

    def extract(self, stacks):
        """Extract data for given stacks. Returns the slab data in meshContainer."""
        t, data = self.extractor.extract(stacks)
        if t == []:
            return []
        var = self.fieldName
        connectivity = self.extractor.getConnectivityArray()
        x, y = self.extractor.getMeshNodeCoords()

        # TODO remove bad values (nan mask?)
        data[data < extractStation.VALID_MIN] = np.nan
        # TODO add support for removing/keeping dry elements
        data = data.swapaxes(0, 1)  # from (dim,np,time) to (np,dim,time)
        ta = timeArray.timeArray(
            t, 'simulation', self.extractor.startTime).asEpoch()
        msldepth = ''
        if self.alongSCoord:
            msldepth = 'slev' + str(self.k)
            z = 0.0 * np.ones_like(x)
        else:
            zSign = 1 if self.zRelativeToSurf else -1  # zRelToSurf => depth below surface
            msldepth = str(int(round(zSign * self.z * 100)))
            z = self.z * np.ones_like(x)

        # if suspected bad values, print warning
        hasBadValues = np.isnan(data).any() or np.isinf(
            data).any() or np.any(data < VALID_MIN)
        if hasBadValues:
            print 'Warning: bad values in', self.name
        meta = {}
        meta['dataType'] = 'slab'
        meta['location'] = self.name
        meta['instrument'] = 'model'
        meta['variable'] = var
        if self.alongSCoord:
            meta['slevel'] = self.k
        else:
            meta['bracket'] = 'F' if self.zRelativeToSurf else 'A'
            meta['msldepth'] = msldepth
        mc = meshContainer.meshContainer(
            '', ta, x, y, z, data, connectivity, fieldNameList[var],
            coordSys='spcs', metaData=meta)
        return mc

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
        '-n',
        '--name',
        action='store',
        type='string',
        dest='name',
        help='name of the slab for identification (default %default)',
        default='slab')
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
        '-z',
        '--zCoord',
        action='store',
        type='float',
        dest='zCoord',
        help='vertical position of slab: z coordinate (z<0 below datum)')
    parser.add_option(
        '-k',
        '--kLevel',
        action='store',
        type='int',
        dest='kLevel',
        help='vertical position of slab: S level, 1 means bottom level, -1 surface')
    parser.add_option(
        '-S',
        '--zRelToSurf',
        action='store_true',
        dest='zRelativeToSurf',
        help='z coordinate (z>0) is depth below free surface (default %default)',
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
        help='Number of tracers to support for \'sed\' and \'generic\' models',
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
    name = options.name
    modelCoordSys = options.modelCoordSys
    outDir = options.outDir
    startStr = options.startStr
    endStr = options.endStr
    stackStr = options.stackStr
    zCoord = options.zCoord
    kLevel = options.kLevel
    zRelativeToSurf = options.zRelativeToSurf
    runTag = options.runTag
    readNetcdf = options.readNetcdf
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
    # if not name :
        # parser.print_help()
        #parser.error('slab name undefined')
    if zCoord is None and kLevel is None:
        parser.print_help()
        parser.error('zCoord or kLevel must be defined')
    if not runTag:
        parser.print_help()
        parser.error('runTag  undefined')

    if zRelativeToSurf and zCoord < 0:
        parser.error('zCoord must be >0 for zRelToSurf option')

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
    if name:
        print ' - name ', name
    if stackStr is None:
        print ' - time range:', str(startTime), '->', str(endTime)
    else:
        print ' - stacks:', str(stacks[0]), '->', str(stacks[-1])
    print ' - dataDir', dataDir
    print ' - SELFE output format:', 'netCDF' if readNetcdf else 'binary'
    print ' - output dir', outDir
    print ' - variables ', varList
    print ' - model coord system', modelCoordSys
    if zCoord:
        print ' - z coordinate', zCoord
        if zRelativeToSurf:
            print ' - z relative to surface'
    else:
        print ' - S coordinate level', kLevel
    sys.stdout.flush()

    dcs = []
    if readNetcdf:
        from crane.data.ncExtract import extractSlabForLevel as extractNetCDF
        dcs = extractNetCDF(
            dataDir,
            varList,
            startTime,
            endTime,
            name,
            zCoord,
            kLevel,
            zRelativeToSurf,
            stacks=stacks)
    else:
        dcs = extractSlabForLevel(dataDir, varList, startTime, endTime,
                                  name, zCoord, kLevel, zRelativeToSurf,
                                  modelCoordSys=modelCoordSys)
    for dc in dcs:
        dc.setMetaData('tag', runTag)
    import crane.data.dirTreeManager as dtm
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
