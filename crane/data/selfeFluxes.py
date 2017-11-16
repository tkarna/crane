"""
Flux computation routine for SELFE netcdf outputs.

1)   Computes volumes and volume fluxes from dihv.71 files (always)
2a)  If tracers are defined (-v option), computes tracer fluxes from 2D
     depth and time integrated files di_trcr_1_flux.65.nc
2b)  Alternatively uses hvel.67 and trcr_1.70 files to approximate the tracer
     flux instead (option --use-hvel). This is less accurate.
3)   Corrects volume and tracer fluxes to match changes in volume/mass.
     For tracers this means that one needs to compute mean tracer in each region
     (from trcr_1.70.nc file) which takes some time. To omit use --no-correction
     option. This correction may be needed even for *.65 fluxes to account for
     losses due to wetting-drying.

Tuomas Karna 2014-10-10
"""

import os
import sys
import datetime
import traceback
import numpy as np
import numpy.linalg
import time as timeMod

from crane.files import gr3Interface
from crane.data import ncExtract
from crane.data import dataContainer
from crane.data import meshContainer
from crane.data import timeArray
from crane.data import gridUtils
from crane.data import extractStation
from crane.data import dirTreeManager
from crane.physicalVariableDefs import addTracers

# -----------------------------------------------------------------------------
#  classes
# -----------------------------------------------------------------------------


class interface(object):
    """
    interface represents all the edges between two regions.

       normal ^
              |  elem2
         n1---x---n2
                 elem1

    n1 - node1
    n2 - node2
    x  - edgenode

    Members:
    --------

    nodes : array_like (nEdges, 2)
            list of nodes in each edge of the interface
    elems : array_like (nEdges, 2)
            list of elements on both sides of the interface
            first element is in the side of lower region tag.
    edges : array_like (nEdges, )
            list of all edge nodes on the interface
    normal_x : array_like (nEdges, )
    normal_y : array_like (nEdges, )
            unit normal vector defined from the lower region to higher
    edge_len : array_like (nEdges, )
            length of each edge
    """

    def __init__(self, elems, meshSearchObj):
        """
        Parameters
        ----------
        elems : array_like (nEdges, 2)
            list of elements on both sides of the interface
            first element is in the side of lower region tag.
        meshSearchObj : meshSearch2d object
            mesh search object for a file that has edge information (hvel.67)
        """
        self.elems = elems
        edgesHi = meshSearchObj.elem2edge[elems[:, 0], :]
        edgesLo = meshSearchObj.elem2edge[elems[:, 1], :]

        # find the edge that appears in both
        edgesAll = np.concatenate((edgesLo, edgesHi), axis=1)
        b = np.sort(edgesAll, axis=1)
        btmp = b[:, :-1]
        edges = btmp[b[:, 1:] == b[:, :-1]]

        node_x = meshSearchObj.node_x
        node_y = meshSearchObj.node_y
        faceNodes = meshSearchObj.faceNodes

        nodes = meshSearchObj.edgeNodes[edges]
        # compute normal
        dx = (node_y[nodes[:, 1]] - node_y[nodes[:, 0]])
        dy = -(node_x[nodes[:, 1]] - node_x[nodes[:, 0]])
        edgeLen = np.hypot(dx, dy)
        dx /= edgeLen
        dy /= edgeLen
        # compute direction from high to low region
        elemsHi = elems[:, 0]
        elemsLo = elems[:, 1]
        cdx = np.mean(node_x[faceNodes[elemsLo, :]], axis=1) -\
            np.mean(node_x[faceNodes[elemsHi, :]], axis=1)
        cdy = np.mean(node_y[faceNodes[elemsLo, :]], axis=1) -\
            np.mean(node_y[faceNodes[elemsHi, :]], axis=1)
        # fix cases where normal is in wrong direction
        flipDir = cdx * dx + cdy * dy < 0
        nodes[flipDir, :] = nodes[flipDir, ::-1]
        dx[flipDir] *= -1
        dy[flipDir] *= -1

        self.nodes = nodes
        self.edgenodes = edges
        self.normal_x = dx
        self.normal_y = dy
        self.edge_len = edgeLen

        uni_nodes, inv_ix = np.unique(self.nodes, return_inverse=True)
        self.uniq_nodes = uni_nodes
        self.uniq_nodes_inv_ix = inv_ix


def buildInterfaces(elemToRegion, meshSearchObj):
    """
    Constructs interface data structures for each interface between regions.

    Parameters
    ----------
    elemToRegion : array_like (nElems, )
            Maps each element to a region (int)
    meshSearchObj : meshSearch2d object
            Mesh data structure object, must contain edge information

    Returns
    -------
    intefaces : dict of interface objects
            All interfaces between regions, indentified by
            (ibox_hi, ibox_lo) tuple of ints.
    """
    elem2neigh = meshSearchObj.elem2neigh
    neighRegions = elemToRegion[elem2neigh]

    # find all elements that lie in an interface between regions
    interfaceElements = neighRegions == np.tile(elemToRegion[:, None], (1, 3))
    # omit boundaries where neigh=-1
    interfaceElements[elem2neigh == -1] = True
    interfaceElements = np.sum(interfaceElements, axis=1) != 3
    interfaceElements = np.nonzero(interfaceElements)[0]

    # maps an interface (iReg_high, iReg_low) to elements in the regions
    # first element is in the higher region, the latter in the lower
    interfaceToElems = dict()
    for iE in interfaceElements:
        # interface is (iReg_high, iReg_low)
        for i in xrange(3):
            if elem2neigh[iE, i] == -1:
                continue
            if elemToRegion[iE] == neighRegions[iE, i]:
                continue
            all_regions = np.array([elemToRegion[iE], neighRegions[iE, i]])
            sortedIx = np.argsort(all_regions)[::-1]
            key = all_regions[sortedIx]
            elems = np.array([iE, elem2neigh[iE, i]])[sortedIx]
            interfaceToElems.setdefault(tuple(key), set()).add(tuple(elems))

    # convert from set of tuples to ndarray
    # interfaceToElems[(high, low)] = [[elem1_high, elem1_low],
    #                                  [elem2_high, elem2_low],
    #                                  ...                     ]
    for k in interfaceToElems:
        interfaceToElems[k] = np.array([list(pair) for pair in
                                        interfaceToElems[k]])

    # build interface objects
    # dictionary of all interfaces, keys are (iReg_high, iReg_low)
    interfaces = {}
    for k in sorted(interfaceToElems.keys()):
        face = interface(interfaceToElems[k], meshSearchObj)
        interfaces[k] = face
    return interfaces


def makeFluxDataContainer(times, volFlux, location, runTag,
                          variable='volume flux', prefix='flux'):
    """
    Creates a single dataContainer that contains all the flux time series
    """
    ta = timeArray.timeArray(times, 'epoch')
    x = np.array([0])
    y = np.array([0])
    z = np.array([0])
    nFluxes = len(volFlux)
    data = np.zeros((1, nFluxes, len(times)))
    fieldNames = []
    for i, k in enumerate(sorted(volFlux)):
        data[0, i, :] = volFlux[k]
        fieldNames.append(prefix + ' from ' + str(k[0]) + ' to ' + str(k[1]))
    meta = {}
    meta['location'] = location
    meta['instrument'] = 'model'
    meta['variable'] = variable
    meta['dataType'] = 'flux'
    meta['tag'] = runTag

    dc = dataContainer.dataContainer(
        '',
        ta,
        x,
        y,
        z,
        data,
        fieldNames,
        coordSys='spcs',
        metaData=meta)
    return dc


def makeVolumeDataContainer(times, volume, location, runTag,
                            variable='volume', prefix='vol'):
    """
    Creates a single dataContainer that contains all the volume time series
    """
    ta = timeArray.timeArray(times, 'epoch')
    x = np.array([0])
    y = np.array([0])
    z = np.array([0])
    nVolumes = len(volume)
    data = np.zeros((1, nVolumes, len(times)))
    fieldNames = []
    for i, k in enumerate(sorted(volume)):
        data[0, i, :] = volume[k]
        fieldNames.append(prefix + ' ' + str(k))
    meta = {}
    meta['location'] = location
    meta['instrument'] = 'model'
    meta['variable'] = variable
    meta['dataType'] = 'flux'
    meta['tag'] = runTag

    dc = dataContainer.dataContainer(
        '',
        ta,
        x,
        y,
        z,
        data,
        fieldNames,
        coordSys='spcs',
        metaData=meta)
    return dc


def correctFluxes(fluxDC, volDC):
    """
    Corrects fluxes so that cumsum matches the changes in volume.
    """
    var = fluxDC.getMetaData('variable')

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
    A = np.zeros((nBoxes, nFluxes))
    for i in xrange(nBoxes):
        iBox = volMap[i]
        incoming_ix = fluxMap[:, 1] == iBox
        outgoing_ix = fluxMap[:, 0] == iBox
        A[i, incoming_ix] = +1.0
        A[i, outgoing_ix] = -1.0
    # volchange and flux time series
    dt = np.diff(volDC.time.array)[0]
    fluxes = fluxDC.data[0, :, :]
    vol = volDC.data[0, :, :]
    volchange = np.diff(vol)
    # omit the first flux value, flux[1] corresponds to vol[1]-vol[0] diff
    fluxes = fluxes[:, 1:]
    volchange_from_fluxes = np.dot(dt * A, fluxes)
    # conservation error
    volError = volchange - volchange_from_fluxes

    # Assume that box 0 is the "background" and omit it from the correction
    A = A[1:, :]
    volError = volError[1:, :]

    # omit fluxes that are close to zero (e.g. salt in river)
    meanFlux = np.mean(np.abs(fluxes), axis=1)
    include_flux_ix = meanFlux > 1.0
    A = A[:, include_flux_ix]
    fluxes = fluxes[include_flux_ix, :]
    if ~(include_flux_ix.any()):
        raise Exception('Cannot correct fluxes: all fluxes are zero, ' + var)

    # solve equation
    # dt*A*(fluxes + correction) - volchange = 0
    # dt*A*fluxes + dt*A*correction - volchange = 0
    # volchange_from_fluxes + dt*A*correction - volchange = 0
    # correction = inv(A)/dt*(volchange - volchange_from_fluxes)

    scaleFluxes = True
    if scaleFluxes:
        # scale fluxes by the mean value to distribute correction
        # proportionally to each flux
        # NOTE this has an effect iff nFluxes > nBoxes
        fluxes_scaled = fluxes / meanFlux[include_flux_ix][:, np.newaxis]
        A_scaled = A * meanFlux[include_flux_ix][np.newaxis, :]
        correction, res, rank, s = numpy.linalg.lstsq(dt * A_scaled, volError)
        correction = correction * meanFlux[include_flux_ix][:, np.newaxis]
    else:
        correction, res, rank, s = numpy.linalg.lstsq(dt * A, volError)

    # dump to new dataContainer
    newFluxDC = fluxDC.copy()
    newFluxDC.data[0, include_flux_ix, 1:] = fluxes + correction
    newFluxDC.setMetaData('variable', 'corr' + var)
    # have a beer
    return newFluxDC

# -----------------------------------------------------------------------------
# main routine
# -----------------------------------------------------------------------------


def computeSelfeFluxes(path, regionFile, location, runTag, stacks=None,
                       startTime=None, endTime=None,
                       trcrVarList=[], useHVel=True, applyCorrection=False):

    verbose = False
    regionMC = gr3Interface.readGR3FileToMC(regionFile)

    regions = np.unique(regionMC.data).astype(int)
    nRegions = len(regions)
    print 'regions', regions

    elemToRegion = np.max(regionMC.data[regionMC.connectivity, 0, 0],
                          axis=1).astype(int)
    regionToElem = {}
    for iReg in regions:
        regionToElem[iReg] = np.nonzero(elemToRegion == iReg)[0]

    fname = ncExtract.getNCVariableName('dihv')
    dihvReader = ncExtract.selfeExtractBase(
        path, fname, fileTypeStr='71', verbose=verbose)
    # construct edge info and mesh search object
    edgeNodes = gridUtils.constructEdgeNodeArray(dihvReader.dataFile.faceNodes)
    meshSearchObj = gridUtils.meshSearch2d(dihvReader.dataFile.nodeX,
                                           dihvReader.dataFile.nodeY,
                                           dihvReader.dataFile.faceNodes,
                                           edgeNodes=edgeNodes)
    trcrReaders = {}
    if useHVel or applyCorrection:
        for v in trcrVarList:
            fname = ncExtract.getNCVariableName(v)
            trcrReaders[v] = ncExtract.selfeExtractBase(
                path, fname, fileTypeStr='70', verbose=verbose,
                meshSearchObj=meshSearchObj)
    if useHVel:
        fname = ncExtract.getNCVariableName('hvel')
        hvelReader = ncExtract.selfeExtractBase(path, fname, fileTypeStr='67',
                                                verbose=verbose,
                                                meshSearchObj=meshSearchObj)
    else:
        trcrFluxReaders = {}
        for v in trcrVarList:
            fname = ncExtract.getNCVariableName(v)
            fname = 'di_' + fname + '_flux'
            trcrFluxReaders[v] = ncExtract.selfeExtractBase(
                path, fname, fileTypeStr='65', verbose=verbose,
                meshSearchObj=meshSearchObj)

    if regionMC.connectivity.shape[0] != meshSearchObj.faceNodes.shape[0]:
        errstr = 'nb elems: '
        errstr += str(regionMC.connectivity.shape[0])
        errstr += ' != '
        errstr += str(meshSearchObj.faceNodes.shape[0])
        raise Exception('Region file doesn\'t match output mesh, ' + errstr)

    if stacks is None:
        stacks = dihvReader.dataFile.getStacks(
            startTime, endTime, wholeDays=True)

    interfaces = buildInterfaces(elemToRegion, meshSearchObj)

    nTime = dihvReader.dataFile.nTime
    nFluxes = len(interfaces)
    nStacks = len(stacks)

    # time series of volume flux across each interface, k=(iReg_high, iReg_low)
    volFlux = {}
    trcrFlux = {}
    for v in trcrVarList:
        trcrFlux[v] = {}

    # time series of volume of each region, key=iReg
    volume = {}
    meanTrcr = {}
    if applyCorrection:
        for v in trcrVarList:
            meanTrcr[v] = {}

    # allocate time series arrays
    times = np.zeros((nTime * nStacks,))
    # flux arrays
    arrayList = []
    arrayList.append(volFlux)
    for v in trcrFlux:
        arrayList.append(trcrFlux[v])
    for array in arrayList:
        for k in interfaces.keys():
            array[k] = np.zeros((nTime * nStacks,))
    # volume arrays
    arrayList = []
    arrayList.append(volume)
    for v in meanTrcr:
        arrayList.append(meanTrcr[v])
    for array in arrayList:
        for iReg in regions:
            array[iReg] = np.zeros((nTime * nStacks,))

    vCoords = dihvReader.dataFile.vCoords
    bathymetry = dihvReader.dataFile.bath
    bathMC = meshContainer.meshContainer('', timeArray.timeArray(np.array([0]), 'epoch'),
                                         dihvReader.dataFile.nodeX,
                                         dihvReader.dataFile.nodeY,
                                         np.zeros_like(dihvReader.dataFile.nodeX),
                                         bathymetry[:, None, None],
                                         dihvReader.dataFile.faceNodes, ['bathymetry'])
    tri_areas = bathMC.computeAreas()
    bathymetry_elem = (bathymetry[meshSearchObj.faceNodes[:, 0]] +
                       bathymetry[meshSearchObj.faceNodes[:, 1]] +
                       bathymetry[meshSearchObj.faceNodes[:, 2]]) / 3.0

    def compute_bottom_z(bathymetry, meshSearchObj, vCoords):
        """Compute bottom z coordinate for each column of prisms."""
        faceNodes = meshSearchObj.faceNodes
        eta0 = np.zeros_like(bathymetry)
        Z, kbottom, iwet = vCoords.computeVerticalCoordinates(eta0, bathymetry)
        k_bottom_elem = np.max(kbottom[faceNodes], axis=1)
        z_bottom_elem = np.zeros_like(k_bottom_elem, dtype=float)
        for iElem in xrange(len(k_bottom_elem)):
            k = k_bottom_elem[iElem]
            z_bottom_elem[iElem] = np.mean(Z[k, faceNodes[iElem, :]])
        return z_bottom_elem, k_bottom_elem

    z_bottom_elem, k_bottom_elem = compute_bottom_z(bathymetry, meshSearchObj,
                                                    vCoords)

    nFaces = meshSearchObj.nFaces
    nNodes = meshSearchObj.nNodes
    Z = np.zeros((vCoords.nvrt, nNodes))
    kbottom = np.zeros((nNodes,), dtype=int)
    iwet = np.zeros((nNodes,), dtype=int)

    Z_elem = np.zeros((vCoords.nvrt, nFaces))
    kbot_elem = np.zeros((nFaces,), dtype=int)
    iwet_elem = np.zeros((nFaces,), dtype=int)

    t0 = timeMod.clock()
    compVCoords = vCoords.computeVerticalCoordinates
    #@profile

    def excecute(Z, kbottom, Z_elem, kbot_elem):
        for iStack in xrange(len(stacks)):
            print 'stack', iStack
            try:
                dihvFile = dihvReader.getNCFile(iStack=stacks[iStack])
                etaDataSet = dihvFile.variables['elev']
                diuDataSet = dihvFile.variables['dihv_x']
                divDataSet = dihvFile.variables['dihv_y']
                trcrDataSet = {}
                if useHVel or applyCorrection:
                    for v in trcrVarList:
                        trcrFile = trcrReaders[v].getNCFile(
                            iStack=stacks[iStack])
                        ncVar = ncExtract.getNCVariableName(v)
                        trcrDataSet[v] = trcrFile.variables[ncVar]
                if useHVel and trcrVarList:
                    uvFile = hvelReader.getNCFile(iStack=stacks[iStack])
                    uDataSet = uvFile.variables['hvel_x']
                    vDataSet = uvFile.variables['hvel_y']
                else:
                    trcrFluxArray = {}
                    for v in trcrVarList:
                        fluxFile = trcrFluxReaders[
                            v].getNCFile(iStack=stacks[iStack])
                        ncVar = ncExtract.getNCVariableName(v)
                        ncVar = 'di_' + ncVar + '_flux'
                        trcrFluxArray[v] = fluxFile.variables[ncVar][:]

                time = dihvFile.getTime()
                times[iStack * nTime:iStack * nTime + nTime] = time[:nTime]

                etaAll = etaDataSet[:]
                diuAll = diuDataSet[:]
                divAll = divDataSet[:]

                for iTime in xrange(nTime):
                    eta = etaAll[iTime, :]  # etaDataSet[iTime,:]
                    diu = diuAll[iTime, :]  # diuDataSet[iTime,:]
                    div = divAll[iTime, :]  # divDataSet[iTime,:]

                    # compute dry mask and mask values for dry nodes/elems
                    faceNodes = meshSearchObj.faceNodes
                    idry_elem, idry = vCoords.computeDryElemMask(
                        eta, bathymetry, faceNodes)
                    eta[idry] = np.nan
                    eta_elem = (eta[faceNodes[:, 0]] +
                                eta[faceNodes[:, 1]] +
                                eta[faceNodes[:, 2]]) / 3.0
                    diu[idry] = 0
                    div[idry] = 0

                    if useHVel or applyCorrection:
                        trcr = {}
                        for v in trcrVarList:
                            trcr[v] = trcrDataSet[v][
                                iTime, :, :]  # layers, iElem
                    if applyCorrection:
                        Z_elem, kbot_elem, _ = compVCoords(eta_elem,
                                                           bathymetry_elem,
                                                           Z_elem, kbot_elem)
                    if useHVel and trcrVarList:
                        Z, kbottom, _ = compVCoords(eta, bathymetry, Z,
                                                    kbottom)

                    # compute volumes and total tracers
                    for iReg in regions:
                        elems = regionToElem[iReg]
                        # compute total volume
                        h = eta_elem[elems] - z_bottom_elem[elems]
                        goodIx = np.isfinite(h) & ~idry_elem[elems]
                        tot_vol = np.sum(h[goodIx] * tri_areas[elems[goodIx]])
                        volume[iReg][iStack * nTime + iTime] = tot_vol
                        if applyCorrection and trcrVarList:
                            # compute element sizes
                            tmp = Z_elem[:, elems]
                            h_elem = tmp[1:, :] - tmp[:-1, :]
                            #h_elem = Z_elem[1:, elems] - Z_elem[:-1, elems]
                            vol_elem = h_elem * tri_areas[elems]
                            for v in trcrVarList:
                                tot = np.sum(trcr[v][1:, elems] * vol_elem)
                                meanTrcr[v][iReg][
                                    iStack * nTime + iTime] = tot / tot_vol

                    # compute fluxes
                    for k in interfaces.keys():
                        face = interfaces[k]
                        idry_edge = idry[
                            face.nodes[:, 0]] | idry[
                            face.nodes[:, 1]]
                        diu_edge = (diu[face.nodes[:, 0]] +
                                    diu[face.nodes[:, 1]]) / 2
                        div_edge = (div[face.nodes[:, 0]] +
                                    div[face.nodes[:, 1]]) / 2
                        # flux from dihv files
                        diu_flux = (
                            diu_edge *
                            face.normal_x +
                            div_edge *
                            face.normal_y)
                        diu_flux[idry_edge] = 0
                        vol_flux = diu_flux

                        volFlux[k][
                            iStack *
                            nTime +
                            iTime] = np.sum(
                            vol_flux *
                            face.edge_len)

                        if useHVel and trcrVarList:
                            # compute tracer fluxes for hvel and tracer values
                            # compute normal u trough each prism quad
                            u_edge = uDataSet[iTime, :, face.edgenodes]
                            v_edge = vDataSet[iTime, :, face.edgenodes]
                            # normal velocity
                            un = u_edge * face.normal_x + v_edge * face.normal_y
                            un[:, idry_edge] = 0
                            # compute flux
                            Z_edge = np.mean(Z[:, face.nodes], axis=2)
                            tot_depth = Z_edge[-1, :] - Z_edge[0, :]
                            quad_height = Z_edge[1:, :] - Z_edge[:-1, :]
                            quad_height = quad_height.filled(0)

                            un_quad = 0.5 * (un[1:, :] + un[:-1, :])
                            upwind_ix = (un_quad > 0).astype(int)

                            for v in trcrVarList:
                                trcr_vals = trcr[v][1:, face.elems]
                                trcr_upwind = trcr_vals[:, :, 0] * upwind_ix +\
                                    (1 - upwind_ix) * trcr_vals[:, :, 1]
                                # trcr_upwind = trcr_vals.mean(axis=2) # as in
                                # fortran
                                trcr_flux = np.sum(
                                    trcr_upwind * un_quad * quad_height, axis=0)
                                trcrFlux[v][k][
                                    iStack *
                                    nTime +
                                    iTime] = np.sum(
                                    trcr_flux *
                                    face.edge_len)
                        else:
                            for v in trcrVarList:
                                # use depth integrated flux field
                                trcr_flux = trcrFluxArray[
                                    v][iTime, face.edgenodes]
                                # correct sign (di_flux points from low elem to
                                # hi)
                                invSign = face.elems[:, 0] > face.elems[:, 1]
                                trcr_flux[invSign] *= -1
                                #trcr_flux[idry_edge] = 0
                                trcrFlux[v][k][
                                    iStack * nTime + iTime] = np.sum(trcr_flux)
            except Exception as e:
                print 'failed:'
                traceback.print_exc(file=sys.stdout)

    excecute(Z, kbottom, Z_elem, kbot_elem)
    print 'duration', timeMod.clock() - t0, 's'

    dcs = []
    fluxDC = makeFluxDataContainer(times, volFlux, location, runTag,
                                   'volumeflux', 'flux')
    print fluxDC
    dcs.append(fluxDC)
    volDC = makeVolumeDataContainer(times, volume, location, runTag,
                                    'volume', 'vol')
    print volDC
    dcs.append(volDC)

    if applyCorrection:
        try:
            corrFluxDC = correctFluxes(fluxDC, volDC)
            print corrFluxDC
            dcs.append(corrFluxDC)
        except Exception as e:
            print 'Flux correction failed:'
            traceback.print_exc(file=sys.stdout)

    trcrFluxDCs = []
    for v in trcrFlux:
        trcrFluxDC = makeFluxDataContainer(times, trcrFlux[v], location,
                                           runTag, v + 'flux', 'flux')
        print trcrFluxDC
        dcs.append(trcrFluxDC)
        trcrFluxDCs.append(trcrFluxDC)

    for i, v in enumerate(meanTrcr):
        meanTrcrDC = makeVolumeDataContainer(times, meanTrcr[v], location,
                                             runTag, 'mean' + v, 'mean' + v)
        print meanTrcrDC
        dcs.append(meanTrcrDC)

        if applyCorrection:
            try:
                totTrcrDC = meanTrcrDC.copy()
                totTrcrDC.data = meanTrcrDC.data * volDC.data
                corrTrcrFluxDC = correctFluxes(trcrFluxDCs[i], totTrcrDC)
                print corrTrcrFluxDC
                dcs.append(corrTrcrFluxDC)
            except Exception as e:
                print 'Flux correction failed:'
                traceback.print_exc(file=sys.stdout)

    return dcs

#-------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------


def parseCommandLine():
    from optparse import OptionParser

    usage = """%prog [options]
    Computes volume/tracer fluxes from SELFE netcdf outputs:
1)   Computes region volumes and volume fluxes from dihv.71 files (always)
2a)  If tracers are defined (-v option), computes tracer fluxes from 2D
     depth and time integrated files di_trcr_1_flux.65.nc
2b)  Alternatively uses hvel.67 and trcr_1.70 files to approximate the tracer
     flux instead (option --use-hvel). This is less accurate.
3)   Corrects volume and tracer fluxes to match changes in volume/mass.
     For tracers this means that one needs to compute mean tracer in each region
     (from trcr_1.70.nc file) which takes some time. To omit use --no-correction
     option. This correction may be needed even for *.65 fluxes to account for
     losses due to wetting-drying.

     Units are:
     volume: [m3], volume flux: [m3/s]
     mean tracer: [trcr_unit], tracer flux: [trcr_unit m3/s]
"""
    parser = OptionParser(usage=usage)
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
    parser.add_option('-S', '--stacks', action='store', type='string',
                      dest='stackStr', help='range of output files to read '
                      '(e.g 1,14) if start,end not given')
    parser.add_option(
        '-t',
        '--regionFile',
        action='store',
        type='string',
        dest='regionFile',
        help='gr3 file defining the regions as different integers',
        default=None)
    parser.add_option(
        '-n',
        '--name',
        action='store',
        type='string',
        dest='name',
        help='name of the region configuration for identification (e.g. myRegions01)')
    parser.add_option(
        '',
        '--save-in-tree',
        action='store_true',
        dest='saveInTree',
        help='saves extracted data in file tree with monthly files instead of a single file (default %default)',
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
        '-v',
        '--tracer-variables',
        action='store',
        type='string',
        dest='trcrVars',
        help='tracer variable(s) to process, e.g. salt,NO3,NH4. '
        'Will compute both fluxes across regions and mean values in each region.')
    parser.add_option(
        '',
        '--use-hvel',
        action='store_true',
        dest='useHVel',
        help='Compute tracer fluxes from trcr.70.nc and hvel.67.nc data, instead of di_trcr_flux.65.nc',
        default=False)
    parser.add_option(
        '', '--no-correction', action='store_false', dest='applyCorrection',
        help='Do not correct fluxes to match changes in volume. '
        'In this case mean tracers are not computed, which saves '
        ' (a lot) of cpu time.', default=True)

    (options, args) = parser.parse_args()

    runTag = options.runTag
    dataDir = options.dataDir
    startStr = options.startStr
    endStr = options.endStr
    stackStr = options.stackStr
    regionFile = options.regionFile
    saveInTree = options.saveInTree
    name = options.name
    tracerModel = options.tracerModel
    numTracers = options.numTracers
    trcrVars = options.trcrVars
    useHVel = options.useHVel
    applyCorrection = options.applyCorrection

    def error(status, msg):
        parser.print_help()
        sys.stderr.write('\nerror: ' + msg + '\n')
        sys.exit(status)

    if not dataDir:
        error(2, 'dataDir  undefined')
    if not startStr and not stackStr:
        error(2, 'stacks or startStr must be defined')
    if not endStr and not stackStr:
        error(2, 'stacks or startStr must be defined')
    if not runTag:
        error(2, 'runTag  undefined')
    if not name:
        error(2, 'name  undefined')
    if not regionFile:
        error(2, 'regionFile  undefined')
    if tracerModel:
        if not numTracers and tracerModel.split('.')[0] in ['sed', 'generic']:
            parser.print_help()
            error(
                2, 'numTracers must be provided if sed or generic tracer models are used.')
        extraTrcrFiles = addTracers(tracerModel, numTracers=numTracers)
    if startStr and endStr:
        startTime = datetime.datetime.strptime(startStr, '%Y-%m-%d')
        endTime = datetime.datetime.strptime(endStr, '%Y-%m-%d')
    else:
        startTime = None
        endTime = None
    if stackStr:
        limits = [int(v) for v in stackStr.split(',')]
        stacks = np.arange(limits[0], limits[1] + 1)
    else:
        stacks = None
    if trcrVars is not None:
        trcrVarList = [v for v in trcrVars.split(',')]
    else:
        trcrVarList = []

    print 'Parsed options:'
    if stacks is None:
        print ' - time range:', str(startTime), '->', str(endTime)
    else:
        print ' - stacks:', stacks
    print ' - runTag:', runTag
    print ' - dataDir:', dataDir
    print ' - regions:', name, regionFile
    if useHVel:
        print ' - tracer fluxes: from hvel.67.nc and trcr.70.nc data'
    else:
        print ' - tracer fluxes: from di_trcr_flux.65.nc fields'
    if trcrVarList:
        print ' - tracers:', trcrVarList
    if trcrVarList and tracerModel:
        print ' - tracer model:', tracerModel

    dcs = computeSelfeFluxes(dataDir, regionFile, name, runTag, stacks=stacks,
                             startTime=startTime, endTime=endTime,
                             trcrVarList=trcrVarList, useHVel=useHVel,
                             applyCorrection=applyCorrection)

    if saveInTree:
        rule = 'monthlyFile'
    else:
        rule = 'singleFile'
    dirTreeManager.saveDataContainerInTree(dcs, rule=rule, dtype=np.float32,
                                           overwrite=True, compress=True)

if __name__ == '__main__':
    parseCommandLine()
