"""
Computes Geyer-MacCready estuary classification parameters from a cross-section
transect.

Tuomas Karna 2016-08-01
"""
import os
import numpy as np
import datetime
import sys
from optparse import OptionParser
import glob
from scipy.interpolate import interp1d, griddata, splev, splrep

import crane.data.dirTreeManager as dtm
import crane.data.dataContainer as dataContainer
import crane.data.timeArray as timeArray
from crane.data.timeSeriesFilters import removeTides, runningMax
from crane.plotting.transectPlot import dcToTransect


def computeNormalVelocity(X, Y, Z, U, V, Xdelta, Ydelta):
    """Computes velocity normal to the transect"""
    Umean = 0.5*(U[:, :-1, :]+U[:, 1:, :])
    Vmean = 0.5*(V[:, :-1, :]+V[:, 1:, :])
    Xmean = 0.5*(Y[:, :-1]+X[:, 1:])
    Ymean = 0.5*(Y[:, :-1]+Y[:, 1:])
    Zmean = 0.5*(Z[:, :-1, :]+Z[:, 1:, :])
    nZ, nX, nTime = Umean.shape
    Unormal = Umean*Ydelta - Vmean*Xdelta
    Vnormal = Umean*Xdelta + Vmean*Ydelta # tangential
    UNormalMean = Unormal.mean()
    print 'Mean normal velocity', UNormalMean
    if UNormalMean > 0 :
        # ensure that mean (residual) is negative, i.e. seaward
        Unormal *= -1
    return Xmean, Ymean, Zmean, Unormal, Vnormal


def computeStratification(S):
    """Computes stratification as max bottom value - min surface value"""
    S_max_bot = np.nanmax(S[0, :, :], axis=0)
    S_min_surf = np.nanmin(S[-1, :, :], axis=0)
    return S_max_bot - S_min_surf


def computeDepth(Zm):
    """Compute transect depth from z coordinates"""
    Zdelta = np.diff(Zm, axis=0)
    dZ = np.zeros_like(Zm)
    dZ[1:, :, :] = 0.5*Zdelta
    dZ[:-1, :, :] = 0.5*Zdelta
    depth = Zm[-1, :, :]-Zm[0, :, :]
    dZ[ ~np.isfinite(dZ) ] = 0.0
    depth[ ~np.isfinite(depth) ] = 0.0
    return depth, dZ


def computeSectionalAverage(Un, depth, dZ, sectLen):
    """Averages field Un over the transect area"""
    transectLen = np.sum(sectLen)
    depthAv = np.sum(Un*sectLen*dZ, axis=0)/depth
    depthAv[~np.isfinite(depthAv)] = 0
    sectAv = np.sum(depthAv, axis=0)/transectLen
    return sectAv


def computeSectionalMaximum(Un, depth, dZ, sectLen):
    """Computes maximum value of field Un over the transect"""
    return np.max(np.max(Un, axis=0), axis=0)


def computeSectionalMaxUbf(U, V, depth, dZ, sectLen):
    """
    Computes max. bottom friction velocit over the transect

    NOTE uses the same bottom friction formula as used in SELFE DB33
    """
    Umag = np.sqrt(U**2 + V**2)
    botU = np.abs(Umag[1, :, :])
    h = dZ[0, :, :]*2
    roughArr = np.zeros_like(h)*1e-4
    roughArr[np.mean(depth, axis=1)>12.0,:] = 1e-6
    Cd = 1/(2.5*np.log(h/roughArr))**2
    botU = np.sqrt(Cd*botU**2)
    botU = np.nanmax(botU,axis=0)
    return botU


def makeDataContainer(t, d, runTag, location, variable):
    ta = timeArray.timeArray(t,'epoch').asEpoch()
    data = d[None,None,:]
    meta = {}
    meta['location'] = location
    meta['instrument'] = 'model'
    meta['variable'] = variable
    meta['dataType'] = 'gmcparams'
    meta['msldepth'] = '0'
    meta['tag'] = runTag
    dc = dataContainer.dataContainer('', ta, 0, 0, 0, data,
                                     [variable], coordSys='', metaData=meta)
    return dc


def reshapeTransect(dc):
    """Converts dataContainer transect into a (XY,Z,U,T) arrays"""
    T, X, Y, Z, U = dcToTransect(dc, iComp=0)
    if dc.data.shape[1] > 1:
        T, X, Y, Z, V = dcToTransect(dc, iComp=1)
    else:
        V = None
    x = X[0,:]
    y = Y[0,:]
    Xdiff = np.diff(x)
    Ydiff = np.diff(y)
    Xdelta = Xdiff
    Ydelta = Ydiff
    sectLen = np.sqrt(Xdelta**2+Ydelta**2)
    Xdelta /= sectLen
    Ydelta /= sectLen
    Xdelta = Xdelta[None,:,None]
    Ydelta = Ydelta[None,:,None]
    sectLen = sectLen[None,:,None]
    return T, X, Y, Z, U, V, Xdelta, Ydelta, sectLen


def generateGMcData(hvelDC, saltDC):
    """
    Computes Geyer-MacCready parameters from the given transects

    The parameters are:

    M = sqrt( u_bf**2 / (omega*N_0*H**2 ) )
    Fr_f = U_R / sqrt( Beta*g*s_ocean*H )

    where

    U_R : residual seaward flow speed (~Q_R/A_section)
    Beta ~ 7.7e-4 : haline contraction coefficient, rho = rho_0*(1+ Beta*salt)
    g ~ 9.81 : gravitation acceleration
    s_ocean ~34.0 : ocean salinity; max salt value in system
    u_bf : bottom friction velocity
    omega = 2*pi*1/T_M2: tidal frequency
    N_0 = sqrt(Beta*g*s_ocean/H): maximal buoyancy frequency
    H : datum depth

    If saltDC is not provided, assumes constant stratification, i.e.
    stratification = s_ocean
    """

    # constants
    Beta = 7.7e-4
    g = 9.81
    s_ocean = 34.0
    Cd = 1.0e-3 #2.5e-3
    T_M2 = 44714.
    omega = 2*np.pi*1/T_M2

    T, X, Y, Z, U, V, Xdelta, Ydelta, sectLen = reshapeTransect(hvelDC)

    # compute normal velocity (nVert, nPoint, nTime)
    Xm, Ym, Zm, Un, Vn = computeNormalVelocity(X, Y, Z, U, V, Xdelta, Ydelta)

    depth, dZ = computeDepth(Zm)
    # compute sectionally averaged velocity
    sectAvU = computeSectionalAverage(Un, depth, dZ, sectLen)
    sectMaxU = computeSectionalMaximum(Un, depth, dZ, sectLen)
    sectMaxUbf = computeSectionalMaxUbf(Un, Vn, depth, dZ, sectLen)

    sectAvUDC = makeDataContainer(T, sectAvU, hvelDC.getMetaData('tag'),
                                  hvelDC.getMetaData('location'), 'sectAvU')
    sectMaxUDC = makeDataContainer(T, sectMaxU, hvelDC.getMetaData('tag'),
                                   hvelDC.getMetaData('location'), 'sectMaxU')
    sectMaxUbfDC = makeDataContainer(T, sectMaxUbf, hvelDC.getMetaData('tag'),
                                     hvelDC.getMetaData('location'), 'sectMaxUbf')



    # bottom friction velocity (max over tidal day)
    ubf = runningMax(sectMaxUbfDC, T=2*44714.)
    ubf = removeTides(ubf)
    ubf.setMetaData('variable', 'ubf')
    ubf.fieldNames=['ubf']

    # generate common time stamps
    new_time = ubf.time  # this is the shortest

    # compute tidal average => residual velocity U_R
    ur = removeTides(sectAvUDC).interpolateInTime(new_time)
    ur.data *= -1  # flip sign
    ur.setMetaData('variable', 'ur')
    ur.fieldNames=['ur']

    if saltDC is not None:
        T, X, Y, Z, S, _, Xdelta, Ydelta, sectLen = reshapeTransect(saltDC)
        strat = computeStratification(S)
        stratDC = makeDataContainer(T, strat, saltDC.getMetaData('tag'),
                                    saltDC.getMetaData('location'), 'strat')
        stratDC = removeTides(stratDC).interpolateInTime(new_time)
        deltaS = stratDC.data
    else:
        deltaS = np.ones_like(ubf.data)*s_ocean

    # compute H (constant in time)
    transectLen = np.sum(sectLen)
    D = depth.copy()
    D[~np.isfinite(D)]=0
    # datum depth (does not necessarily correspond to actual datum!)
    H = np.mean(np.sum(D*sectLen[0,:,:],axis=0)/transectLen)

    # depth time series
    depth_dc = ur.copy()
    depth_dc.data[:] = H
    depth_dc.setMetaData('variable', 'depth')
    depth_dc.fieldNames=['depth']

    # denominator of Fr, wave celecity
    uc = ur.copy()
    uc.data[:] = np.sqrt(Beta*g*deltaS*H)
    uc.setMetaData('variable', 'uc')
    uc.fieldNames=['uc']

    # Fr_f
    Fr_f = ur.copy()
    Fr_f.data /= uc.data[:]
    Fr_f.setMetaData('variable', 'froude')
    Fr_f.fieldNames=['froude']

    # max buoyancy frequency
    N_0 = np.sqrt(Beta*g*deltaS/H)

    # stratification
    strat_dc = ubf.copy()
    strat_dc.data[:] = deltaS
    strat_dc.setMetaData('variable', 'strat')
    strat_dc.fieldNames=['strat']

    # buoyancy frequency
    buoy_dc = ubf.copy()
    buoy_dc.data[:] = N_0
    buoy_dc.setMetaData('variable', 'buoy')
    buoy_dc.fieldNames=['buoy']

    # denominator of M
    um = ubf.copy()
    um.data[:] = np.sqrt(omega*N_0*H**2)
    um.setMetaData('variable', 'um')
    um.fieldNames=['um']

    # M
    M = ubf.copy()
    M.data /= um.data[:]
    M.setMetaData('variable', 'mixing')
    M.fieldNames=['mixing']

    def printMinMax(dc, name):
        if isinstance(dc, dataContainer.dataContainer):
            d = dc.data.ravel()
        else:
            d = dc
        print name, d.min(), d.mean(), d.max()

    print M
    print Fr_f

    return M, Fr_f, ubf, um, ur, uc, depth_dc, strat_dc, buoy_dc


# -----------------
# user interface
# -----------------

def parse_commandline():
    import argparse

    parser = argparse.ArgumentParser(description='Computes Geyer-MacCready estuary classification parameters from cross-channel transect.')

    parser.add_argument('-v', required=True, dest='hvelfile',
                        help='Horizontal velocity (hvel) transect file to process')
    parser.add_argument('-s', dest='saltfile',
                        help='Salinity (salt) transect file to process')
    parser.add_argument('--constant-strat', action='store_true',
                        help='Use a constant stratification (s_ocean - 0.0), '
                        'instead of computing instantaneous stratification '
                        'from the salinity transect. '
                        'NOTE: This is default use in many applications.')


    args = parser.parse_args()


    if args.constant_strat:
        salt_dc = None
    else:
        assert args.saltfile is not None, 'salt transect file must be provided'
        salt_dc = dataContainer.dataContainer.loadFromNetCDF(args.saltfile)
    hvel_dc = dataContainer.dataContainer.loadFromNetCDF(args.hvelfile)

    dc_list = generateGMcData(hvel_dc, salt_dc)

    rule='singleFile'
    dtm.saveDataContainerInTree(list(dc_list), rule=rule, overwrite=True)

if __name__ == '__main__':
    parse_commandline()
