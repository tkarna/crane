#!/usr/bin/env python
"""
Collection of methods for manipulating mesh data.

Tuomas Karna 2013-02-15
"""
import os
import sys
import time as timeMod

import numpy as np
from scipy.spatial import KDTree, cKDTree

from crane.data import meshContainer
from crane.files import gmshInterface
from crane.files import gr3Interface
from crane.files import gmshInterface


def readAnyMeshFile(inFile, dataFile=None):
    """Reads mesh data in a meshContainer. Supported formats are GR3 (SELFE), MSH (GMSH) and meshContainer netCDF."""
    fname, ext = os.path.splitext(inFile)
    if ext == '.gr3':
        mc = gr3Interface.readGR3FileToMC(inFile)
    elif ext == '.msh':
        gmsh = gmshMesh.fromFile(inFile)
        if dataFile:
            mc = gmsh.getMeshContainer(fromFile=dataFile)
        else:
            mc = gmsh.getMeshContainer(fromFile=inFile)
    elif ext == '.nc':
        mc = meshContainer.loadFromNetCDF(inFile)
    else:
        raise Exception('Unknown file extension: ' + ext + '\n' + inFile)
    return mc


#-------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------
class verticalCoordinates:
    """
    A class for computingComputes SELFE vertical coordinates.
    """

    def __init__(
            self,
            nvrt,
            nz,
            h_s,
            h_c,
            theta_b=None,
            theta_f=None,
            ztot=None,
            sigma=None,
            cs=None,
            h0=0.01):
        self.nvrt = nvrt
        self.nz = nz
        self.h0 = h0
        self.h_s = h_s
        self.h_c = h_c
        self.theta_b = theta_b
        self.theta_f = theta_f
        self.ztot = ztot
        self.sigma = sigma
        self.cs = cs
        if cs is None and sigma is not None:
            self.cs = self.compute_cs(sigma)

    @classmethod
    def fromVGridFile(cls, filename, h0=0.01):
        """Parses vgrid file and returns a verticalCoordinates object"""
        print 'Reading ' + filename + ' ...'
        f = open(filename, 'r')
        words = f.readline().split()
        nvrt = int(words[0])
        nz = int(words[1])
        h_s = float(words[2])
        print '  nvrt=', nvrt, 'nz=', nz, 'h_s=', h_s
        f.readline()  # 'Z levels'
        ztot = np.zeros(nvrt)
        sigma = np.zeros(nvrt)
        for i in range(nz):
            words = f.readline().split()
            # print i,words
            ztot[i] = float(words[1])
        f.readline()  # 'S levels'
        words = f.readline().split()
        h_c = float(words[0])
        theta_b = float(words[1])
        theta_f = float(words[2])
        print '  h_c=', h_c, 'theta_b=', theta_b, 'theta_f=', theta_f
        for i in range(nvrt - nz):
            words = f.readline().split()
            # print i,words
            sigma[i] = float(words[1])
        return cls(nvrt, nz, h_s, h_c, theta_b, theta_f, ztot, sigma, h0=h0)

    def compute_cs(self, sigma):
        return ((1 -
                 self.theta_b) *
                np.sinh(self.theta_f *
                        sigma) /
                np.sinh(self.theta_f) +
                self.theta_b *
                (np.tanh(self.theta_f *
                         (sigma +
                          0.5)) -
                 np.tanh(self.theta_f *
                         0.5)) /
                2 /
                np.tanh(self.theta_f *
                        0.5))

    def computeDryMask(self, eta, dp):
        return eta + dp < self.h0

    def computeDryElemMask(self, eta, dp, connectivity):
        """
        Computes the true drymask where
        a) element is wet if all its nodes are wet
        b) node is wet if it touches at least one wet element
        """
        # dry nodes by depth
        dryNodes = eta + dp < self.h0
        dryElems = np.max(dryNodes[connectivity].copy(), axis=1)
        # a node is wet only if at least one elem touching it is wet
        tmp = np.tile(dryElems, (3, 1)).T
        ix = tmp == 0
        wetNodes = np.unique(connectivity[ix].ravel())
        dryNodes[:] = True
        dryNodes[wetNodes] = False
        return dryElems, dryNodes

    def computeVerticalCoordinates2(self, eta, dp, z=None, kbp=None,
                                    iwet=None):
        """
        Cython implementation.

        See also
        --------
        computeVerticalCoordinates: equivalent numpy implementation
        """
        import selfeVerticalCoords
        nNodes = len(eta)
        if z is None:
            z = np.zeros((self.nvrt, nNodes))
        if kbp is None:
            kbp = np.zeros((nNodes,), dtype=int)
        if iwet is None:
            iwet = np.zeros((nNodes,), dtype=int)
        selfeVerticalCoords.vertCoords(
            eta,
            dp,
            z,
            kbp,
            iwet,
            self.ztot,
            self.cs,
            self.sigma,
            self.nz,
            self.h0,
            self.h_c,
            self.h_s)
        z = np.ma.masked_invalid(z)
        return z, kbp, iwet

    def computeVerticalCoordinates(self, eta, dp, ztmp=None, kbp2=None):
        """
        Parameters
        ----------
        eta : array_like (nPoints,)
              water elevation at certain points
        dp  : array_like (nPoints,)
              bathymetry at same points

        Returns
        -------
        Z : array_like (nVert,nPoints)
            z coordinates
        kbp2 : array_like (nPoint,)
            bottom level indices
        iwet : array_like (nPoint,)
            wet mask (per node, isolated wet nodes are not removed like in SELFE)
        """

        # Compute z coordinates (only for 3d)
        #sys.stdout.write('Generating vertical coordinates...')
        # sys.stdout.flush()
        nNodes = len(eta)
        if ztmp is None:
            ztmp = np.zeros((self.nvrt, nNodes))
        if kbp2 is None:
            kbp2 = np.zeros((nNodes,), dtype=int)
        ztmp[:, :] = np.nan

        iwet = eta + dp >= self.h0

        nz = max(self.nz, 1)  # workaround for cases with no z levels
        hc = self.h_c  # surf layer depth
        hs = self.h_s  # s-z transition depth
        hmod = np.minimum(dp, hs)
        ztmp[nz - 1, iwet] = -hmod[iwet]
        ztmp[self.nvrt - 1, iwet] = eta[iwet]

        # check validity of v.grid (positive surf. layer depth)
        nSigma = self.nvrt - self.nz
        #etaMin = -hc-(hmod-hc)*self.cs[nSigma-1]/self.sigma[nSigma-1]
        #idx = np.logical_and( eta <= etaMin, iwet )
        # if idx.any() :
        ## raise error
        # errStr = ' '.join( ['Error in v.coords: choose a larger h_c :',
        # str(hc), str(eta[idx][0]),str(etaMin[idx][0])] )
        #raise Exception( errStr )

        # case 1: dp <= hc, shallower than fine surface layer
        idx1 = np.logical_and(hmod <= hc, iwet)
        H1 = (hmod[idx1] + eta[idx1])
        eta1 = eta[idx1]
        # case 2: deeper than surf layer (normal case)
        idx2 = np.logical_and(hmod > hc, iwet)
        eta2 = eta[idx2]
        hmod2 = hmod[idx2] - hc
        for k in xrange(nz, self.nvrt - 1):
            kin = k - nz + 1
            sigma = self.sigma[kin]
            ztmp[k, idx1] = sigma * H1 + eta1  # NOTE slow
            # NOTE slow
            ztmp[
                k, idx2] = eta2 * (1 + sigma) + hc * sigma + hmod2 * self.cs[kin]
        # find bottom indices
        # pure sigma part
        idx = dp <= hs
        kbp2[idx] = nz - 1
        # find z-coordinates
        idx_zpart = (dp > hs)
        for k in xrange(nz - 1):
            idx = idx_zpart & (-dp >= self.ztot[k]) & (-dp < self.ztot[k + 1])
            kbp2[idx] = k
            ztmp[k, idx] = -dp[idx]
            idx = idx_zpart & (-dp < self.ztot[k])
            ztmp[k, idx] = self.ztot[k]
        # extend z coordinates for shaved cells
        #hGapShaved = 0.5
        # for k in range(nz,-1,-1) :
            #idx = kbp2 == k
            # if idx.any() :
            # for kk in range(k-1,-1,-1) :
            #ztmp[kk,idx] = ztmp[kk+1,idx] - hGapShaved
        ztmp = np.ma.masked_invalid(ztmp)
        # print '  done.'
        return ztmp, kbp2, iwet

#-------------------------------------------------------------------------
# Functions (NOTE: these are obsolete, use gridUtils.py instead)
#-------------------------------------------------------------------------


def constructKDTree(triX, triY):
    """Constructs KDTree for the given mesh node coordinates.

    KDTree facilitates quick search for nearest nodes in the grid:
    nix = tree.query( np.array([[x,y]]) )[1][0]

    Args:
      triX,triY -- (ndarray (nX,)) mesh node coordinates

    Returns:
      tree -- (scipy cKDTree) KDtree object built for the mesh node coordinates
    """
    return cKDTree(np.concatenate((triX[:, None], triY[:, None]), axis=1))


def constructNodeToElemTable(conn):
    """Constructs node to element table for the given connectivity array

    Args:
      conn -- (ndarray (nElems,3)) mesh connectivity array
    Returns:
      node2elem -- (dict of lists) a dictionary that maps node index to list of elements that share that node
   """
    node2elem = {}
    #t0 = timeMod.clock()
    nElems = conn.shape[0]
    for iE in range(nElems):
        for iN in conn[iE, :]:
            node2elem.setdefault(iN, []).append(iE)
    # print 'node2elem duration', timeMod.clock() - t0, 's'
    return node2elem


def findTriangle(tree, node2elem, triX, triY, conn, x, y, verbose=True):
    """Finds an element in the triangular mesh that contains the given (x ,y) point.

    Args:
      tree -- (scipy cKDTree) KDtree object built for the mesh node coordinates
      node2elem -- (dict of lists) a dictionary that maps node index to list of elements that share that node
      triX,triY -- (ndarray (nX,)) mesh node coordinates
      conn -- (ndarray (nElems,3)) mesh connectivity array
      x,y -- (float) query point coordinates

    Returns:
      nix -- (ndarray (1,3)) node indices of the vertices of the triangle
      iTri -- (integer) index of the triangle
      u -- (ndarray (1,3)) barycentric coordinates of (x,y) point in the triangle
    """
    # find nearest node
    nix = tree.query(np.array([[x, y]]))[1][0]
    # get all elements touching the node
    eix = node2elem[nix]
    # expand stencil: get all nodes associated with these elements
    nix = np.unique(conn[eix, :].flatten())
    # get all possible elements
    eix = []
    for i in nix:
        eix.extend(node2elem[i])
    eix = list(np.unique(np.array(eix, dtype=int)))
    # compute barycentric coords for those elements
    u = barycentricCoords(triX[conn[eix, :]], triY[conn[eix, :]], x, y)
    # check which triangle has the point
    it = np.nonzero(isInsideTri(x, y, u=u))[0]
    if len(it) == 0:
        if verbose:
            print 'parent triangle not found'
        return None, None, None
    if len(it) > 1:
        # multiple parent triangles found, most likely sharing a node, take 1st
        it = it[0][None]
    # get global elem ix
    iTri = eix[it[0]]
    return conn[iTri, :], iTri, u[it]


def findStations(
        tree,
        node2elem,
        triX,
        triY,
        conn,
        staX,
        staY,
        stationNames=None,
        verbose=True):
    """Finds given station locations in the mesh.

    If any station is out of the grid, a warning message is printed and the station is ignored.

    Args:
      tree -- (scipy cKDTree) KDtree object built for the mesh node coordinates
      node2elem -- (dict of lists) a dictionary that maps node index to list of elements that share that node
      triX,triY -- (ndarray (nX,)) mesh node coordinates
      conn -- (ndarray (nElems,3)) mesh connectivity array
      staX,staY -- (ndarray (nSta,)) station coordinates
      stationNames -- (list) string name of each station (optional)

    Returns:
      iTri -- (ndarray (nOKSta,)) Indices of parent triangles for each good station
      u -- (ndarray (nOKSta,3)) Barycentric coordnates for each good station
      n -- (ndarray (nOKSta,3)) Triangle vertex indices for each good station
      goodIx -- (ndarray (nOKSta,)) Indices of good stations vrt original staX array
    """
    nSta = len(staX)
    iTri = np.zeros((nSta,), dtype=int)
    u = np.zeros((nSta, 3))
    n = np.zeros((nSta, 3), dtype=int)
    goodStaIx = np.zeros((nSta,), dtype=bool)
    for i in range(nSta):
        nn, it, uu = findTriangle(
            tree, node2elem, triX, triY, conn, staX[i],
            staY[i],
            verbose)
        if it is not None:
            iTri[i] = it
            u[i, :] = uu
            n[i, :] = nn
            goodStaIx[i] = True
        else:
            sta = stationNames[i] if stationNames else ''
            if verbose:
                print 'station out of grid, skipping:', sta, staX[i], staY[i]
    goodStaIx = np.nonzero(goodStaIx)[0]
    iTri = iTri[goodStaIx]
    u = u[goodStaIx, :]
    n = n[goodStaIx, :]
    return iTri, u, n, goodStaIx


def evalOnTriangulation(tree, node2elem, triX, triY, conn, nodalVals, x, y):
    """Evaluates a field in the given location (x,y) over the triangulation.

    Args:
      triX,triY -- (ndarray (nX,)) mesh node coordinates
      conn -- (ndarray (nElems,3)) mesh connectivity array
      nodalVals -- (ndarray (nX,)) Nodal value of the field to evaluate
      x,y -- (float) query point coordinates
      tree -- (scipy cKDTree) KDtree object built for the mesh node coordinates
      node2elem -- (dict of lists) a dictionary that maps node index to list of elements that share that node

    Returns:
      val -- (float) value of the field at location (x,y)
    """

    nix, iTri, u = findTriangle(tree, node2elem, triX, triY, conn, x, y)
    val = interpolateTri(u, nodalVals[conn[iTri, :]])

    return val


def interpolateTri(u, nodalValues):
    """Interpolates nodal values in a location given by barycentric coordinates u
    for all the elements.
    Args:
      nodalValues -- (nTri,3)
      u           -- (1,3) or (nTri,3)
    Returns:
      values      -- (nTri,1)
    """
    nTri = nodalValues.shape[0]
    values = np.ma.ones((nTri,)) * np.nan
    nv = np.ma.masked_invalid(nodalValues)
    goodIx = np.all(np.logical_not(nv.mask), axis=1)
    if not goodIx.any():
        return np.ma.masked_invalid(values[:, None])
    nv = nv[goodIx, :]
    if u.shape[0] == 1 and nodalValues.shape[0] > 1:
        U = np.tile(u, (nv.shape[0], 1))
    else:
        U = u
    values[goodIx] = np.sum(u * nv, axis=1)
    return np.ma.masked_invalid(values[:, None])


def isInsideTri(x, y, u=None, nodeX=None, nodeY=None,):
    """ Tests whether (x,y) is in a triangle  whose vertices are nodeX,nodeY.
    Args:
      NodeX,nodeY -- array (nTri,3)
      x,y         -- scalar
    Returns:
      result      -- bool array (nTri,)
    """
    if u is None:
        u = barycentricCoords(nodeX, nodeY, x, y)
    return np.logical_and(u >= 0, u <= 1).all(axis=1)


def barycentricCoords(nodeX, nodeY, x, y):
    """Returns barycentric coordinates for (x,y) in triangle whose vertices are
    nodeX,nodeY.
    Args:
      NodeX,nodeY -- array (nTri,3)
      x,y         -- scalar
    Returns:
      u           -- array (nTri,3)
    """
    detT = ((nodeY[:, 1] - nodeY[:, 2]) * (nodeX[:, 0] - nodeX[:, 2]) -
            (nodeX[:, 1] - nodeX[:, 2]) * (nodeY[:, 0] - nodeY[:, 2]))
    u1 = ((nodeY[:, 1] - nodeY[:, 2]) * ((x - nodeX[:, 2])) -
          (nodeX[:, 1] - nodeX[:, 2]) * (y - nodeY[:, 2])) / detT
    u2 = ((nodeY[:, 2] - nodeY[:, 0]) * ((x - nodeX[:, 2])) -
          (nodeX[:, 2] - nodeX[:, 0]) * (y - nodeY[:, 2])) / detT
    u3 = 1.0 - u1 - u2
    u = np.array([u1, u2, u3]).T
    return u
