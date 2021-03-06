#!/usr/bin/env python
"""
Script for generating bottom friction Cd plots.

Tuomas Karna 2013-03-01
"""

import os
import sys
import numpy as np
import datetime
from scipy.interpolate import interp1d

from crane.data import timeArray
from crane.data import dataContainer
from crane.data import dirTreeManager
from crane.plotting import transectPlot

TRANSECT_NAME = 'mainChannel'
TRANSECT_BPFILE = '/home/workspace/users/pturner/db29/processing.dev/scripts.working/intrusion_length.bp'


def compute_intrusion_length(x, salt, salt_threshold, min_gap_length=5.0, verbose=False):
    """Computes salt intrusion length."""
    dx_along = np.hstack(([0], np.diff(x)))
    # Detect gaps, i.e. nodes where S is below the threshold
    binary_ix = salt > salt_threshold
    if binary_ix.all():
        if verbose:
            print('whole transect is within threshold')
        x_threshold = x[-1]
        return x_threshold

    # pad beginning assuming that ocean is always above threshold
    # this ensures that there's at least one gap start
    binary_ix = np.hstack(([True], binary_ix))
    # detect gaps
    transition = np.diff(binary_ix.astype(int))
    # first node below threshold
    gap_start = np.nonzero(transition == -1)[0]
    # first node above threshold
    gap_end = np.nonzero(transition == 1)[0]
    # compute the length of each gap
    gap_len = np.zeros_like(gap_start)
    for i in range(len(gap_start)):
        if i < len(gap_end):
            gap_len[i] = np.sum(dx_along[gap_start[i]:gap_end[i]])
            if verbose:
                print('gap {:} len: {:}'.format(i, gap_len[i]))
        else:
            # gap end missing, assume infinite gap
            gap_len[i] = 1.0e10
    # index of SIL node is the gap_start_ix - 1
    # for the first gap whose length is above min_gap_length
    good_gaps = np.nonzero(gap_len > min_gap_length)[0]
    if len(good_gaps) == 0:
        if verbose:
            print('no suitable gap found, assume maximal SIL')
        x_threshold = x[-1]
        return x_threshold

    ix_sil_node = max(gap_start[good_gaps[0]] - 1, 0)
    if verbose:
        print('good gap found: {:}'.format(ix_sil_node))
    if ix_sil_node == 0:
        if verbose:
            print('gap found in the beginning, SIL = 0')
        x_threshold = x[0]
        return x_threshold
    # interpolate distance
    func = interp1d(salt[[ix_sil_node + 1, ix_sil_node]],
                    x[[ix_sil_node + 1, ix_sil_node]])
    x_threshold = func(salt_threshold)

    if verbose:
        for i in range(len(x)):
            arrow = '<-- SIL' if i == ix_sil_node else ''
            print('{:4d} {:12.3f} {:12.3f} {:}'.format(i, x[i], salt[i], arrow))

    return x_threshold


def computeSaltIntrusion(transectDC, salt_threshold_list, min_gap_length=5.0):
    print 'Generating salt intrusion length... ',
    sys.stdout.flush()
    ntime = len(transectDC.time)
    nthresholds = len(salt_threshold_list)
    sil = [np.zeros((ntime,)) for i in xrange(nthresholds)]
    # for each time step
    for it in range(ntime):
        # convert transect to array
        x_along, Z, salt, time, uniqueXYCoords = transectPlot.generateTransectFromDataContainer(
            transectDC, it)
        x_along = x_along[0, :]
        # compute max salt in each column
        bottom_salt = salt.max(axis=0)
        for i, salt_threshold in enumerate(salt_threshold_list):
            gap_len = min_gap_length*1000.0  # convert to meters
            x_sil = compute_intrusion_length(x_along, bottom_salt,
                                             salt_threshold,
                                             min_gap_length=gap_len)
            sil[i][it] = x_sil / 1000.  # convert to km

    output_dc_list = []
    for i, salt_threshold in enumerate(salt_threshold_list):
        # generate dataContainer
        meta = {}
        meta['location'] = transectDC.getMetaData('location')
        meta['instrument'] = 'model'
        meta['variable'] = '''sil_%d''' % (salt_threshold,)
        meta['dataType'] = 'sil'
        meta['tag'] = transectDC.getMetaData('tag')
        data = sil[i][None, None, :]
        silDC = dataContainer.dataContainer(
            '', transectDC.time, 0, 0, 0, data, [
                meta['variable']], coordSys='', metaData=meta)
        # compute daily max
        dailyMax = []
        dailyTime = []
        dayBegin = silDC.time.getDatetime(0).replace(
            minute=0, hour=0, second=0, microsecond=0)
        while dayBegin < silDC.time.getDatetime(-1):
            dayEnd = dayBegin + datetime.timedelta(hours=24.8)
            try:
                dailySIL = silDC.timeWindow(dayBegin, dayEnd)
                m = dailySIL.data.max()
                dailyMax.append(m)
                dailyMax.append(m)
                dailyTime.append(dailySIL.time.array[0])
                dailyTime.append(dailySIL.time.array[-1])
            except Exception as e:
                print 'Cannot compute SIL: day missing, skipping'
                print e
            dayBegin = dayEnd
        data = np.array(dailyMax)[None, None, :]
        ta = timeArray.timeArray(np.array(dailyTime), 'epoch')
        meta = {}
        meta['location'] = transectDC.getMetaData('location')
        meta['instrument'] = 'model'
        meta['variable'] = '''max_sil_%d''' % (salt_threshold,)
        meta['dataType'] = 'sil'
        meta['tag'] = transectDC.getMetaData('tag')
        dailyMaxDC = dataContainer.dataContainer(
            '', ta, 0, 0, 0, data, [meta['variable']],
            coordSys='', metaData=meta)
        output_dc_list.append(silDC)
        output_dc_list.append(dailyMaxDC)
    return output_dc_list

#-------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------


def parseCommandLine():
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-r', '--runTag', action='store', type='string',
                      dest='runTag',
                      help='Run tag, used as a label in post-proc.')
    parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
    parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
    parser.add_option('-d', '--dataDirectory', action='store',
                      type='string', dest='dataDir', help='directory where model outputs are stored')
    parser.add_option('-C', '--read-netcdf', action='store_true',
                      dest='readNetcdf',
                      help='Extract from SELFE netcdf output files instead of SELFE binary files (default %default)',
                      default=False)
    parser.add_option('-v', '--salinityThreshold', action='store',
                      type='string', dest='salinityThreshold',
                      help='PSU value for determining presence of salinity along the transect line (default %default). Multiple thresholds can be given in a list e.g. 1.0,5.0',
                      default='1.0')
    parser.add_option('--minGapThreshold', action='store', type='float',
                      dest='minGapThreshold',
                      help='Ignore gaps in salinity transect that are shorter than this value (in km) (default %default)',
                      default='5.0')
    parser.add_option('-t', '--buildPointFile', action='store',
                      type='string', dest='bpFile',
                      help='text file (*.bp) containing the transect (x,y) coordinates (optional)',
                      default=None)
    parser.add_option('', '--variable', action='store', type='string',
                      dest='var',
                      help='Define alternative output file process, e.g. salt.70 (default %default)',
                      default='salt')
    parser.add_option('-c', '--modelCoordSys', action='store', type='string',
                      dest='modelCoordSys', default='spcs',
                      help='horizontal coordinate system used in model: '
                      'spcs or utm (Default: %default)')

    (options, args) = parser.parse_args()

    runTag = options.runTag
    startStr = options.startStr
    endStr = options.endStr
    dataDir = options.dataDir
    bpFile = options.bpFile
    modelCoordSys = options.modelCoordSys
    salinityThreshold = options.salinityThreshold
    readNetcdf = options.readNetcdf
    var = options.var
    minGapThreshold = options.minGapThreshold

    salt_threshold_list = [float(vStr)
                           for vStr in salinityThreshold.split(',')]

    if not bpFile:
        # TODO move to shared location
        bpFile = TRANSECT_BPFILE

    if not dataDir:
        parser.print_help()
        parser.error('dataDir  undefined')
    if not startStr:
        parser.print_help()
        parser.error('startStr undefined')
    if not endStr:
        parser.print_help()
        parser.error('endStr   undefined')
    if not runTag:
        parser.print_help()
        parser.error('runTag  undefined')

    startTime = datetime.datetime.strptime(startStr, '%Y-%m-%d')
    endTime = datetime.datetime.strptime(endStr, '%Y-%m-%d')

    print 'Parsed options:'
    print ' - time range:', str(startTime), '->', str(endTime)
    print ' - salinity threshold(s):', salt_threshold_list
    print ' - transect gap threshold:', minGapThreshold
    print ' - dataDir', dataDir
    print ' - variable', var
    print ' - SELFE output format:', 'netCDF' if readNetcdf else 'binary'
    print ' - runTag', runTag
    print ' - transect file', bpFile
    print ' - model coord system', modelCoordSys

    # Extract
    varList = [var]
    name = TRANSECT_NAME
    dcs = []
    if readNetcdf:
        from crane.data.ncExtract import extractTransectForBPFile
        dcs = extractTransectForBPFile(bpFile, dataDir, varList,
                                       startTime, endTime, name)
    else:
        from crane.data.extractTransect import extractTransectForBPFile
        dcs = extractTransectForBPFile(bpFile, dataDir, varList,
                                       startTime, endTime, name=name,
                                       modelCoordSys=modelCoordSys)
    for dc in dcs:
        dc.setMetaData('tag', runTag)

    # compute SIL
    silDCs = computeSaltIntrusion(dcs[0], salt_threshold_list, min_gap_length=minGapThreshold)
    for dc in silDCs:
        print dc

    rule = 'monthlyFile'
    dirTreeManager.saveDataContainerInTree(silDCs, rule=rule, dtype=np.float32,
                                           overwrite=True)

if __name__ == '__main__':
    parseCommandLine()
