"""
Computes and plots residual salinity transport profile from
velocity and salinity profiles.

Tuomas Karna 2013-11-08
"""

import string
import datetime
import numpy as np
from scipy.interpolate import interp1d, griddata
from scipy.spatial import cKDTree
import matplotlib.gridspec as gridspec

from crane import matplotlib
from crane import plt

from crane.data import dataContainer
from crane.data import timeArray
from crane.data import dirTreeManager
from crane.data import pca
from crane.plotting import profilePlot
from crane.physicalVariableDefs import VARS
from crane.physicalVariableDefs import UNITS
from crane.utility import createDirectory
from crane.utility import saveFigure
from crane.utility import parseTimeStr

fontsize = 14
matplotlib.rcParams['font.size'] = fontsize


def cleanPCA(pc_u):
    """Rotates PCA velocity so that flood (in bottom layer) is positive"""
    z = pc_u.z.mean(axis=1)
    ix = z > -10.0
    mvel = pc_u.data[ix, 0, :].mean()
    # print 'mean surf vel',mvel
    if mvel > 0:
        # print 'flipping'
        pc_u.data *= -1
    return pc_u


def doPCAprof(dc, fraction=0.9):
    """Computes PCA directions for the given uv dataContainer"""
    A = np.swapaxes(dc.data[:, :2, :], 1, 2)
    A = np.reshape(A, (-1, dc.data.shape[1]))
    # print "PCA:", p.npc
    # print "% variance:", p.sumvariance * 100
    # print 'dir:',p.Vt
    B = p.pc()
    dc2 = dc.copy()
    dc2.data = np.reshape(B, (dc.data.shape[0], 1, dc.data.shape[2]))
    dc2.fieldNames = dc.fieldNames[:dc2.data.shape[1]]
    dc2.setMetaData('variable', 'pcavel')
    dc2 = cleanPCA(dc2)
    return dc2


def computeMeanProfile(dc):
    """Computes temporal mean"""
    dc2 = dc.copy()
    dc2.data = dc2.data.mean(axis=2)[:, :, None]
    dc2.z = dc2.z.mean(axis=1)
    dc2.zDependsOnTime = False
    dc2.time = timeArray.timeArray(np.array([dc.time.array.mean()]), 'epoch')
    return dc2


def getTidalSaltTransport(salt, hvel):
    """Computes tidal component <u'S'> of salt transport
      s' = salt - <salt>
      <.> is tidal average
    """
    sMean = computeMeanProfile(salt)
    sPrime = salt.copy()
    sPrime.data -= np.tile(sMean.data, (1, 1, sPrime.data.shape[2]))
    uMean = computeMeanProfile(hvel)
    uPrime = hvel.copy()
    uPrime.data -= np.tile(uMean.data, (1, 1, uPrime.data.shape[2]))
    # prod = s'*u'
    prod = sPrime.copy()
    prod.data = sPrime.data * uPrime.data
    # res = <s'*u'>
    res = computeMeanProfile(prod)
    res.setMetaData('label', '<u\'S\'>')
    res.setMetaData('variable', 'sTranspTide')
    return res


def getMeanSaltTransport(salt, hvel):
    """Computes mean component <u><S> of salt transport
      s' = salt - <salt>
      <.> is tidal average
    """
    s = computeMeanProfile(salt)
    u = computeMeanProfile(hvel)
    res = s.copy()
    res.data = s.data * u.data
    res.setMetaData('label', '<u><S>')
    res.setMetaData('variable', 'sTranspMean')
    return res


def getSaltTransport(salt, hvel):
    """Computes salt transport <uS> profile, where
      s' = salt - <salt>
      <.> is tidal average
    """
    res = salt.copy()
    res.data = salt.data * hvel.data
    res.setMetaData('label', '<uS>')
    res.setMetaData('variable', 'sTransport')
    res = computeMeanProfile(res)
    return res

# TODO add to global list of variables ?
VARS['sTransport'] = 'Salt Transport'
UNITS['sTransport'] = 'psu m s-1'
VARS['sTranspTide'] = 'Salt Transport'
UNITS['sTranspTide'] = 'psu m s-1'
VARS['sTranspMean'] = 'Salt Transport'
UNITS['sTranspMean'] = 'psu m s-1'


def makeSubPlot(ax, dcs, startTime, endTime):
    legendProp = {'loc': 'upper left', 'bbox_to_anchor': (1.01, 1.00)}
    var = dcs[0].getMetaData('variable')
    loc = dcs[0].getMetaData('location')
    dia = profilePlot.verticalProfilePlotDC(
        xunit=UNITS.get(var, '-'),
        xlabel=VARS.get(var, var))
    dia.setAxes(ax)
    colorDict = {}
    for dc in dcs:
        tag = dc.getMetaData('tag')
        if tag not in colorDict:
            colorDict[tag] = dia.colorSequence[len(colorDict)]
    for dc in dcs:
        #dc.z = -dc.z
        var = dc.getMetaData('variable')
        tag = dc.getMetaData('tag')
        label = tag
        if dc.hasMetaData('label'):
            label += ' ' + dc.getMetaData('label')
        if var == 'sTransport':
            lStyle = 'solid'
            lWidth = 1.5
        elif var == 'sTranspTide':
            lStyle = '-.'
            lWidth = 1.0
        elif var == 'sTranspMean':
            lStyle = '--'
            lWidth = 1.0
        dia.addSample(
            dc,
            label,
            color=colorDict[tag],
            linestyle=lStyle,
            linewidth=lWidth)
    dia.showLegend(**legendProp)
    start_str = startTime.strftime('%b %d')
    end_str = endTime.strftime('%b %d')
    # , size='large')
    dia.addTitle(' '.join([loc.upper(), start_str, '-', end_str]))

#------------------------------------------------------------------------------
# Main routine
#------------------------------------------------------------------------------


def computeAndPlotSaltTransport(tags, locs, imgDir, startTime, endTime):

    if isinstance(locs, str):
        locs = [locs]

    rule = 'monthlyFile'

    fPrefix = 'resProf'
    filetype = 'png'

    # generate data, store in dict:
    # data['stationName']=[total,mean,tidalTransp]
    data = {}
    for loc in locs:
        dcs = []
        for tag in tags:
            saltPr = dirTreeManager.getDataContainer(
                tag=tag, dataType='profile', location=loc, variable='salt',
                startTime=startTime, endTime=endTime)
            hvelUV = dirTreeManager.getDataContainer(
                tag=tag, dataType='profile', location=loc, variable='hvel',
                startTime=startTime, endTime=endTime)
            hvelPr = doPCAprof(hvelUV.copy())
            dc = getSaltTransport(saltPr, hvelPr)
            dcs.append(dc)
            dc = getMeanSaltTransport(saltPr, hvelPr)
            dcs.append(dc)
            dc = getTidalSaltTransport(saltPr, hvelPr)
            dcs.append(dc)
        data[loc] = dcs

    def addChar(ax, char):
        if char:
            ax.text(0.0, 1.02, char + ')', fontsize=fontsize + 4,
                    verticalalignment='bottom', horizontalalignment='left',
                    transform=ax.transAxes)
    chList = string.ascii_lowercase
    chIter = iter(chList)

    nPlots = len(locs)
    fig = plt.figure(figsize=(5, 3.5 * nPlots))
    Grid = gridspec.GridSpec(nPlots, 1)
    axList = []

    for i in range(nPlots):
        if i == 0:
            axList.append(fig.add_subplot(Grid[i:i + 1, :]))
        else:
            axList.append(fig.add_subplot(Grid[i:i + 1, :], sharex=axList[0]))
    Grid.tight_layout(fig, w_pad=0.0, h_pad=0.7)

    for ax, loc in zip(axList, locs):
        makeSubPlot(ax, data[loc], startTime, endTime)
        if nPlots > 1:
            addChar(ax, chIter.next())

    for ax in axList[:-1]:
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), visible=False)

    locStr = '-'.join(locs)
    varStr = 'saltTransport'
    tagStr = '-'.join(list(set(dc.getMetaData('tag') for dc in dcs)))
    dateStr = startTime.strftime('%Y-%m-%d')
    dateStr += '_' + endTime.strftime('%Y-%m-%d')
    fname = '_'.join([fPrefix, locStr, tagStr, varStr, dateStr])
    saveFigure(imgDir, fname, filetype, verbose=True, dpi=100, bbox_tight=True)
    plt.close()

#------------------------------------------------------------------------------
# Command line interface
#------------------------------------------------------------------------------


def parseCommandLine():

    from optparse import OptionParser
    usage = 'Usage: %prog [options] runTag1 runTag2 ...'
    parser = OptionParser(usage=usage)
    parser.add_option(
        '-a',
        '--stations',
        action='store',
        type='string',
        dest='stations',
        help='list of stations to plot: jetta,saturn01')
    parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Start time of temporal average')
    parser.add_option(
        '-d',
        '--n-days',
        action='store',
        type='float',
        default=2.0,
        dest='days',
        help='Number of tidal days over which temporal average is computed (default:%default)')
    parser.add_option(
        '-o',
        '--output-directory',
        action='store',
        type='string',
        dest='imgDir',
        help='base directory where generated figures are stored')

    (opt, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        parser.error('runTags missing')

    tags = args

    if not opt.stations:
        parser.print_help()
        parser.error('stations undefined')
    if not opt.startStr:
        parser.print_help()
        parser.error('start time undefined')
    if not opt.imgDir:
        parser.print_help()
        parser.error('output directory undefined')

    if opt.startStr:
        startTime = parseTimeStr(opt.startStr)
    endTime = startTime + datetime.timedelta(seconds=opt.days * 2 * 44714.0)

    stations = opt.stations.split(',')

    print 'Parsed options:'
    print ' - tags:', tags
    print ' - stations:', stations
    print ' - time range:', startTime, endTime
    print ' - output directory:', opt.imgDir

    if opt.imgDir:
        createDirectory(opt.imgDir)

    computeAndPlotSaltTransport(tags, stations, opt.imgDir, startTime, endTime)

if __name__ == '__main__':
    parseCommandLine()
