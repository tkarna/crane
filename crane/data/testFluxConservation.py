"""
Plots volume/tracer mass conservation for each region in flux partitioning.
Also plots fluxes and fluxes corrected to match volume changes.

Tuomas Karna 2014-10-29
"""
import os
import numpy as np

from crane import plt

from crane.data import timeArray
from crane.data import dirTreeManager
from crane.data import extractStation
from crane.plotting import timeSeriesPlot
from crane.utility import createDirectory
from crane.utility import saveFigure
from crane.physicalVariableDefs import VARS
from crane.physicalVariableDefs import UNITS
from crane.physicalVariableDefs import addTracers

def processFluxes(runTag, var, location, imgDir):
    fluxname = 'flux'
    volname = 'volume'
    fluxname = var+fluxname
    rule = 'singleFile'
    fluxes = dirTreeManager.getDataContainer(rule=rule, tag=runTag, dataType='flux',
                                  location=location,variable=fluxname)
    corrfluxes = dirTreeManager.getDataContainer(rule=rule, tag=runTag, dataType='flux',
                                      location=location,
                                      variable='corr'+fluxname)
    vol = dirTreeManager.getDataContainer(rule=rule, tag=runTag, dataType='flux',
                               location=location, variable=volname)
    if var != 'volume':
        trcrvolname = 'mean'+var
        trcr = dirTreeManager.getDataContainer(rule=rule, tag=runTag, dataType='flux',
                                    location=location, variable=trcrvolname)
        volume = trcr.copy()
        volume.data = vol.data*trcr.data
    else:
        volume = vol

    print volume
    print fluxes
    print corrfluxes

    nFluxes = corrfluxes.data.shape[1]
    nBoxes = volume.data.shape[1]

    for i in range(nFluxes):
        scalar = max(np.abs(corrfluxes.data[0,i,:]).max(), 1.0)
        err = np.abs((fluxes.data[0,i,:]-corrfluxes.data[0,i,:])/scalar)
        print fluxes.fieldNames[i], 'rel. difference', np.abs(err).max()

    def computeVolError(fluxDC, volDC):
        # [from,to] box ids for each flux time series
        fluxMap = []
        for f in fluxDC.fieldNames:
            words = f.split(' ')
            fromBox = int(words[2])
            toBox = int(words[4])
            fluxMap.append([fromBox, toBox])
        fluxMap = np.array(fluxMap)
        # [boxId] for each volume
        volMap = []
        for f in volDC.fieldNames:
            words = f.split(' ')
            iBox = int(words[1])
            volMap.append(iBox)
        volMap = np.array(volMap)

        # Matrix that computes the total inflow to each box
        # volchange = dt*A*fluxes
        nBoxes = len(volMap)
        nFluxes = fluxMap.shape[0]
        A = np.zeros((nBoxes,nFluxes))
        for i in xrange(nBoxes):
            iBox = volMap[i]
            incoming_ix = fluxMap[:,1] == iBox
            outgoing_ix = fluxMap[:,0] == iBox
            A[i,incoming_ix] = +1.0
            A[i,outgoing_ix] = -1.0

        dt = np.diff(volDC.time.array)[0]
        fluxes = fluxDC.data[0,:,:]
        vol = volDC.data[0,:,:]
        volchange = np.diff(vol)
        # omit the first flux value, flux[1] corresponds to vol[1]-vol[0] diff
        fluxes = fluxes[:,1:]
        volchange_from_fluxes = np.dot(dt*A, fluxes)
        volError = volchange - volchange_from_fluxes

        # dump time series in dataContainers for plotting
        newtime = timeArray.timeArray(volDC.time.array[1:], 'epoch')
        volchangeDC = volDC.copy()
        volchangeDC.time = newtime
        volchangeDC.data = volchange[None, :, :]
        volchangeDC.setMetaData('variable','volchange')
        volchangeDC2 = volDC.copy()
        volchangeDC2.time = newtime
        volchangeDC2.data = volchange_from_fluxes[None, :, :]
        volchangeDC2.setMetaData('variable','volchange')
        volErrorDC = volDC.copy()
        volErrorDC.time = newtime
        volErrorDC.data = volError[None, :, :]
        volErrorDC.setMetaData('variable','volerror')
        return volchangeDC, volchangeDC2, volErrorDC

    volchangeDC, volchange_flx_DC, volErrorDC = computeVolError(fluxes, volume)
    volchangeDC_r, volchange_flx_DC_r, volErrorDC_r = computeVolError(corrfluxes,
                                                                    volume)

    # plot volumes and volume conservation
    for i in range(nBoxes):
        dia = timeSeriesPlot.stackTimeSeriesPlotDC(figwidth=15)
        var = volume.getMetaData('variable')
        volunit = 'm3'
        vollabel = 'volume'
        if var != 'volume':
            primaryVar = var.replace('mean','')
            tracerUnit = UNITS.get(primaryVar,'')
            vollabel = primaryVar+' volume'
            if tracerUnit:
                volunit = tracerUnit+' '+volunit
        varStr = var.replace('flux',' flux')
        dia.addPlot('vol', ylabel=vollabel, unit=volunit, lw=1.5)
        dia.addSample('vol',volume.extractField(i), label='exact', color='k')
        ts = volchange_flx_DC.data[0,i,:]
        data = volume.data[0,i,1] + np.cumsum(ts[:-1])-ts[0]
        tmpDC = volchange_flx_DC.copy()
        tmpDC.data = data
        dia.addSample('vol',tmpDC, label='flux', color='b')
        ts = volchange_flx_DC_r.data[0,i,:]
        data = volume.data[0,i,1] + np.cumsum(ts[:-1])-ts[0]
        tmpDC = volchange_flx_DC.copy()
        tmpDC.data = data
        dia.addSample('vol',tmpDC, label='corr. flux', color='r', linestyle='dashed')
        dia.showLegend()
        dia.addPlot('dvol', ylabel='Vol. change', unit=volunit, lw=1.5)
        dia.addSample('dvol',volchange_flx_DC.extractField(i),
                    label='flux', color='b')
        dia.addSample('dvol',volchange_flx_DC_r.extractField(i),
                    label='corr. flux', color='r', linestyle='dashed')
        dia.addPlot('volerr', ylabel='Vol. error', unit=volunit, lw=1.5)
        dia.addSample('volerr',volErrorDC.extractField(i),
                    label='flux', color='b')
        dia.addSample('volerr',volErrorDC_r.extractField(i),
                    label='corr. flux', color='r', linestyle='dashed')
        dia.addTitle('Box '+volume.fieldNames[i].split(' ')[-1])

        # save fig
        loc = volume.getMetaData('location')
        boxStr = 'volume-'+str(i)
        if var != 'volume':
            boxStr = primaryVar+boxStr
        sT = str(volume.time.getDatetime(0).date())
        eT = str(volume.time.getDatetime(-1).date())
        fname = '_'.join([loc,boxStr,sT,eT])
        saveFigure(imgDir,fname,'png',verbose=True, bbox_tight=True)
        plt.close()

    # plot non-zero fluxes
    for i in range(nFluxes):
        if np.abs(corrfluxes.data[0,i,:]).max() < 1.0:
            continue
        dia = timeSeriesPlot.stackTimeSeriesPlotDC(figwidth=15)
        var = fluxes.getMetaData('variable')
        primaryVar = var.replace('flux','')
        varStr = var.replace('flux',' flux')
        unit = 'm3/s'
        tracerUnit = UNITS.get(primaryVar,'')
        if tracerUnit:
            unit = tracerUnit+' '+unit
        dia.addPlot('flux', ylabel=varStr, unit=unit, lw=1.5)
        dia.addSample('flux',fluxes.extractField(i), label='flux', color='b')
        dia.addSample('flux',corrfluxes.extractField(i), label='corr. flux', color='r', linestyle='dashed')
        dia.showLegend()
        dia.addTitle(fluxes.fieldNames[i].replace('flux',varStr))

        # save fig
        loc = fluxes.getMetaData('location')
        fluxStr = fluxes.fieldNames[i]
        fluxStr = fluxStr.replace(' from ','-').replace(' to ','-')
        if var != 'volume':
            fluxStr = primaryVar+fluxStr
        sT = str(fluxes.time.getDatetime(0).date())
        eT = str(fluxes.time.getDatetime(-1).date())
        fname = '_'.join([loc,fluxStr,sT,eT])
        saveFigure(imgDir,fname,'png',verbose=True, bbox_tight=True)
        plt.close()

#-------------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------------
def parseCommandLine() :
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-r', '--runTag', action='store', type='string',
                        dest='runTag', help='Run tag, used as a label in post-proc.')
    parser.add_option('-o', '--imageDirectory', action='store', type='string',
                        dest='imgDir', help='(optional) directory where generated images are stored. If not specified, shows the image instead.')
    parser.add_option('-v', '--variable', action='store', type='string',
                      dest='var', help='variable to process, either \'volume\' or a tracer name e.g. \'salt\'')
    parser.add_option('-n', '--name', action='store', type='string',
                      dest='name', help='name of the region configuration (e.g. shoreaches)')
    #parser.add_option('-s', '--start', action='store', type='string',
                        #dest='startStr', help='Date to start processing')
    #parser.add_option('-e', '--end', action='store', type='string',
                        #dest='endStr', help='Date to end processing')
    parser.add_option('-T', '--tracerModel', action='store', type='string', dest='tracerModel',
                        help='Tracer model. Must supply number of tracers '
                        'for \'sed\' and \'generic\' models via the -N '
                        'switch.', default=None)
    parser.add_option('-N', '--numTracers', action='store', type='int', dest='numTracers',
                        help='Tracer number to extract for \'sed\' and \'generic\' models',
                        default=None)

    (options, args) = parser.parse_args()

    runTag = options.runTag
    imgDir = options.imgDir
    name = options.name
    var = options.var
    #startStr = options.startStr
    #endStr = options.endStr
    tracerModel = options.tracerModel
    numTracers = options.numTracers

    if not imgDir :
        parser.print_help()
        parser.error('imgDir  undefined')
    if not name :
        parser.print_help()
        parser.error('name  undefined')
    #if not startStr :
        #parser.print_help()
        #parser.error('stacks or startStr must be defined')
    #if not endStr :
        #parser.print_help()
        #parser.error('stacks or startStr must be defined')
    if not runTag :
        parser.print_help()
        parser.error('runTag  undefined')
    if tracerModel :
        if not numTracers and tracerModel.split('.')[0] in ['sed','generic']:
            parser.print_help()
            parser.error('numTracers must be provided if sed or generic tracer models are used.')
        addTracers( tracerModel, numTracers=numTracers)

    print 'Parsed options:'
    print ' - runTag:', runTag
    print ' - output dir', imgDir
    print ' - regions', name

    imgDir = createDirectory(imgDir)
    processFluxes(runTag, var, name, imgDir)

if __name__=='__main__' :
  parseCommandLine()
