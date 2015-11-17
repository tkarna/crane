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
import sys
import datetime
# set Ctrl-C to default signal (terminates fortran routines immediately)
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import numpy as np

from crane.physicalVariableDefs import addTracers
import crane.data.slimNcExtract as SNE

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
    print ' - runTag:', runTag
    if outDir:
        print ' - output dir:', outDir
    print ' - variables:', varList
    if profile:
        print ' - extracting vertical profiles'
    if noOfferings and stationFile is not None:
        print ' - reading stations from:', stationFile
    sys.stdout.flush()

    se = SNE.slimExtract(dataDir, varList[0], verbose=True)
    for var in varList:
        dcs = se.extractForStations(
            dataDir,
            var,
            stationFile,
            startTime,
            endTime,
            profile=profile,
            stacks=stacks,
            verbose=True)
        # if noOfferings :
        #  dcs = se.extractForStations(dataDir, var, startTime,endTime,var,staX,staY,stationNames, staZ=None, k=None, zRelToSurf=None)
        # else :
        #  # get all available offerings
        #  import crane.data.netcdfCacheInterface as netcdfDB
        #  if allOfferings :
        #    print('fetching {0:s} offerings from the database for all stations'.format(var))
        #    offerings = netcdfDB.getAllOfferings( [var.split('.')[0]] )
        #  else :
        #    print('fetching {0:s} offerings from the database for the time period'.format(var))
        #    offerings = netcdfDB.getAvailableOfferings( startTime, endTime, [var.split('.')[0]] )
        #  if len(offerings) == 0:
        #      print('No offerings received, skipping variable {0:s}'.format(var))
        #  dcs = SNE.slimExtract( dataDir, var, offerings, startTime, endTime,
        #                        profile, stationFile)
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
