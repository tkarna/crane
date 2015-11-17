import numpy as np

from crane.data import timeArray
from crane.data import meshContainer


class gmshMesh(object):
    """Object that represents GMSH mesh"""

    def __init__(
            self,
            x,
            y,
            z,
            connectivity,
            physicalNames=None,
            lines=None,
            lineTags=None):
        """Creates a new object.
        Arguments:
        x,y,z - (nNodes,) node coordinates
        connectivity - (nElems,nNodes) connectivity array of elements
        physicalNames - (dict) dictionary of tags identifying each 1D boundary line,
                        2D surface etc. (dim,iName) = 'tag' e.g. (1,1):'open bnd pacific'
        lines - (nLines,2) connectivity table for all 1D boundary elements
        lineTags - (nLine,1) associates each line to physicalName iName
        """
        self.x = x
        self.y = y
        self.z = z
        self.connectivity = connectivity.astype(int)
        self.lines = lines
        self.lineTags = lineTags
        if physicalNames is not None:
            self.physicalNames = physicalNames
        else:
            self.physicalNames = {}
        self.invPhysicalNames = {}
        for key in self.physicalNames:
            self.invPhysicalNames[self.physicalNames[key]] = key

    @staticmethod
    def readNodeBlock(infile):
        nPoints = np.fromfile(infile, dtype=int, count=1, sep=' ')[0]
        nodes = np.zeros((nPoints, 4))
        for i in range(nPoints):
            nodes[i, :] = np.fromfile(infile, dtype=float, count=4, sep=' ')
        return nodes

    @staticmethod
    def readElemBlock(infile):
        ne = np.fromfile(infile, dtype=int, count=1, sep=' ')
        # n columns unknown, must be read line by line
        lines = []  # size unknown, use list for efficiency
        triangles = []
        for i in range(ne):
            # for i in range(10) :
            line = infile.readline().strip()
            #data = [ int(word) for word in line.split() ]
            data = np.array(line.split()).astype(int)
            if data[1] == 1:
                lines.append(data)
            elif data[1] == 2:
                triangles.append(data)
            else:
                print 'Unknown element type: ', data[1]
        lines = np.array(lines) if lines else None
        triangles = np.array(triangles) if triangles else None
        return lines, triangles

    @staticmethod
    def readPhysicalNamesBlock(infile):
        physicalNames = dict()
        n = np.fromfile(infile, dtype=int, count=1, sep=' ')[0]
        for i in range(n):
            line = infile.readline().strip()
            words = line.split()
            if len(words) < 3:
                raise Exception('Corrupted line')
            dim = int(words[0])
            tag = int(words[1])
            name = line.split('"')[1]
            physicalNames[(dim, tag)] = name
        return physicalNames

    @classmethod
    def fromFile(cls, filename):
        """Reads mesh information from msh file."""
        print 'reading mesh from', filename, '...'
        infile = open(filename, 'r')
        nodes = lines = triangles = None
        physicalNames = None
        while True:
            line = infile.readline().strip()
            if line == '':
                break  # end of file
            if line == '$Nodes':
                nodes = cls.readNodeBlock(infile)
            if line == '$Elements':
                lines, triangles = cls.readElemBlock(infile)
            if line == '$PhysicalNames':
                physicalNames = cls.readPhysicalNamesBlock(infile)
        # convert to normal lists
        # nodes
        nodenum = nodes[:, 0]
        x = nodes[:, 1]
        y = nodes[:, 2]
        z = nodes[:, 3]
        # lines
        if lines is not None and len(lines) > 0:
            linenum = lines[:, 0]
            linetype = lines[0, 1]
            nlinetags = lines[0, 2]
            lineTags = lines[:, 3:3 + nlinetags]
            lineNodes = lines[:, -2:] - 1  # two last columns
            # keep only physical tag
            lineTags = lineTags[:, -2]
        else:
            lineNodes = lineTags = None
        # triangles
        trinum = triangles[:, 0]
        tritype = triangles[0, 1]
        ntritags = triangles[0, 2]
        tritags = triangles[3:3 + ntritags]
        triNodes = triangles[:, -3:] - 1  # three last columns

        # for better interpolation, round coordinates to 1e-4
        nDig = 4
        x = np.round(x, nDig)
        y = np.round(y, nDig)
        z = np.round(z, nDig)

        return cls(x, y, z, triNodes, physicalNames, lineNodes, lineTags)

    @classmethod
    def fromMeshContainer(cls, mc):
        """Creates an GMSH mesh object from MeshContainer."""

        # construct physicalNames: (dim,tagNb):'tagStr'
        bndTags = []
        for bnd in mc.boundaries:
            # print bnd.type,bnd.tag
            bndTags.append(bnd.type + ' ' + bnd.tag)
        physicalNames = {}
        dim = 1  # line
        for i, tag in enumerate(bndTags):
            print dim, i + 1, tag
            physicalNames[(dim, i + 1)] = tag
        # mesh plane (dim=2)
        physicalNames[(2, len(physicalNames) + 1)] = 'volume'
        # check duplicates
        #multiplicates = {}
        # for key in physicalNames :
        #t = physicalNames[key]
        # multiplicates.setdefault(t,0)
        #multiplicates[t] += 1
        # if multiplicates[t] > 1 :
        #physicalNames[key] = t+'_'+str(multiplicates[t])
        invPhysicalNames = {}
        for key in physicalNames:
            invPhysicalNames[physicalNames[key]] = key
        if len(physicalNames) != len(invPhysicalNames):
            print 'Warning: mesh physical names are not set correctly'
            print physicalNames
            print invPhysicalNames
        # check boundaries for gaps
        bnds = []
        for bnd in mc.boundaries:
            bnds.append(bnd.copy())
        for bnd in bnds:
            other_bnds = list(bnds)
            other_bnds.remove(bnd)
            for isEndPoint, p in enumerate([bnd.nodes[0], bnd.nodes[-1]]):
                # search for connection
                connected = False
                for ob in other_bnds:
                    if ob.nodes[0] == p or ob.nodes[-1] == p:
                        connected = True
                        break
                print bnd.tag, isEndPoint, connected
                if not connected:
                    # find nearest point
                    distList = []
                    candList = []
                    for ob in other_bnds:
                        for isEndPoint2, q in enumerate(
                                [ob.nodes[0], ob.nodes[-1]]):
                            p2 = bnd.nodes[0] if isEndPoint else bnd.nodes[-1]
                            if q != p2:
                                d = np.sqrt(
                                    (mc.x[p] - mc.x[q]) ** 2 +
                                    (mc.y[p] - mc.y[q]) ** 2 +
                                    (mc.z[p] - mc.z[q]) ** 2)
                                distList.append(d)
                                candList.append((ob, isEndPoint2))
                    distList = np.array(distList)
                    imin = distList.argmin()
                    ob, isEndPoint2 = candList[imin]
                    q = ob.nodes[-1] if isEndPoint2 else ob.nodes[0]
                    print 'this', p, mc.x[p], mc.y[p], mc.z[p]
                    print 'othe', q, mc.x[q], mc.y[q], mc.z[q]
                    if bnd.type != 'open':
                        # add line segment to current bnd
                        if isEndPoint:
                            bnd.nodes = np.hstack((bnd.nodes, [q]))
                        else:
                            bnd.nodes = np.hstack(([q], bnd.nodes))
                    else:
                        # add line segment to other bnd
                        if isEndPoint2:
                            ob.nodes = np.hstack((ob.nodes, [p]))
                        else:
                            ob.nodes = np.hstack(([p], ob.nodes))

        # construct boundary line information
        lines = []
        lineTags = []
        for bnd in bnds:
            nP = len(bnd.nodes)
            nodes = np.zeros((nP - 1, 2), dtype=int)
            dim, tagInt = invPhysicalNames[bnd.type + ' ' + bnd.tag]
            tags = np.ones((nP - 1,), dtype=int) * tagInt
            nodes[:, 0] = bnd.nodes[:-1]
            nodes[:, 1] = bnd.nodes[1:]
            lines.append(nodes)
            lineTags.append(tags)
        # merge arrays
        if lines:
            lines = np.vstack(tuple(lines))
            lineTags = np.hstack(tuple(lineTags))
        else:
            lines = lineTags = None

        return cls(
            mc.x,
            mc.y,
            mc.z,
            mc.connectivity,
            physicalNames,
            lines,
            lineTags)

    def startSection(self, outfile, tag, verbose):
        if verbose:
            print '  ' + tag
        outfile.write('${0:s}\n'.format(tag))

    def endSection(self, outfile, tag, verbose):
        outfile.write('$End{0:s}\n'.format(tag))

    def writeMesh(self, outfilename, verbose=False):
        prismList = []  # TODO update
        outfile = open(outfilename, 'w')
        print 'Writing mesh to', outfilename, '...'
        nNodes = len(self.x)
        if self.lines is not None:
            nLines = self.lines.shape[0]
        else:
            nLines = 0
        nTriangles = self.connectivity.shape[0]
        nPrisms = len(prismList)
        nElements = nLines + nTriangles + nPrisms

        # write header
        self.startSection(outfile, 'MeshFormat', verbose)
        outfile.write('2.2 0 8\n')
        self.endSection(outfile, 'MeshFormat', verbose)
        # physicalNames
        if self.physicalNames:
            self.startSection(outfile, 'PhysicalNames', verbose)
            outfile.write('{0:d}\n'.format(len(self.physicalNames)))
            for dim, tagNb in self.physicalNames.keys():
                phyStr = '"' + self.physicalNames[(dim, tagNb)] + '"'
                outfile.write('{0:d} {1:d} {2:s}\n'.format(dim, tagNb, phyStr))
            self.endSection(outfile, 'PhysicalNames', verbose)

        # nodes section
        self.startSection(outfile, 'Nodes', verbose)
        outfile.write('{0:d}\n'.format(nNodes))
        tmp = np.vstack((np.arange(nNodes) + 1, self.x, self.y, self.z)).T
        np.savetxt(outfile, tmp, fmt=['%d', '%.16g', '%.16g', '%.16g'])
        self.endSection(outfile, 'Nodes', verbose)

        # elements section
        self.startSection(outfile, 'Elements', verbose)
        outfile.write('{0:d}\n'.format(nElements))
        if nLines:
            elemType = 1  # 2 node line
            nTags = 2
            tmp = np.hstack((np.arange(1, nLines + 1)[:, None],
                             np.ones((nLines, 1), dtype=int) * elemType,
                             np.ones((nLines, 1), dtype=int) * nTags,
                             np.tile(self.lineTags[:, None], (1, 2)),
                             self.lines + 1))
            np.savetxt(outfile, tmp, fmt='%d')
        elemType = 2  # 3 node triangle
        nTags = 2
        elemTag = 0
        for key in self.physicalNames:
            if key[0] == 2:  # dim == 2
                elemTag = key[1]
                break
        tmp = np.hstack((np.arange(1, nTriangles + 1)[:, None] + nLines,
                         np.ones((nTriangles, 1), dtype=int) * elemType,
                         np.ones((nTriangles, 1), dtype=int) * nTags,
                         np.ones(
            (nTriangles, 1), dtype=int) * elemTag,  # physical entity
            # geometrical entity
            np.ones((nTriangles, 1), dtype=int) * elemTag,
            self.connectivity + 1))
        np.savetxt(outfile, tmp, fmt='%d')
        for i in range(nPrisms):
            # TODO update
            iElem = prismList[i][0]
            elemType = 6  # 6 node prism
            # 213 6 2 1 12 5 207 155 1 205 4
            nTags = 2
            tags = [11, 7]  # physical entity, geometrical entity
            nodes = prismList[i][1:]
            outfile.write(
                '{i:d} {type:d} {ntags:d} '.format(
                    i=iElem, type=elemType, ntags=nTags))
            for tag in tags:
                outfile.write('{0:d} '.format(tag))
            for n in nodes:
                outfile.write('{0:d} '.format(n))
            outfile.write('\n')
        self.endSection(outfile, 'Elements', verbose)
        outfile.close()

    def getBoundaries(self):
        if self.lines is not None and len(self.physicalNames) > 0:
            return findBoundarySegments(
                self.lines, self.lineTags, self.physicalNames)

    def getMeshContainer(
            self,
            nodalValues=None,
            fieldNames=None,
            fromFile=None):
        if nodalValues is None and fromFile is None:
            raise Exception('either nodalValues or fromFile must be specified')
        if nodalValues is not None and fieldNames is None:
            raise Exception('fieldNames unspecified')
        if nodalValues is None:
            try:
                nodalValues, fieldName, time, it, nScalars = self.readNodalValues(
                    fromFile)
            except Exception as e1:
                try:
                    # discontinuous values in each element
                    elemValues, fieldName, time, it, nScalars = self.readElemValues(
                        fromFile)
                    # convert to continuos by taking a node (by averaging
                    # around nodes)
                except Exception as e2:
                    print self.connectivity.shape
                    print e1
                    print e2
                    raise Exception('Could not read nodal values')
            if fieldNames is None:
                fieldNames = [fieldName]
            # except Exception as e :
                # print 'Warning: could not read nodal values, assigning zeros'
                # print e
                #nodalValues = np.zeros((len(self.x),1))
                #fieldNames = ['zero']
        if isinstance(fieldNames, str):
            fieldNames = [fieldNames]
        description = ''
        ta = timeArray.timeArray(np.array([0]), 'epoch')
        if not isinstance(nodalValues, np.ndarray):
            nodalValues = nodalValues * np.zeros((self.x.shape[0], 1))
        data = nodalValues[:, :, None]
        mc = meshContainer(
            description,
            ta,
            self.x,
            self.y,
            self.z,
            data,
            self.connectivity,
            fieldNames,
            coordSys='spcs')
        # add boundary information
        bnds = self.getBoundaries()
        if bnds is not None:
            for bnd in bnds:
                mc.addBoundary(bnd)
        return mc

    def getMeshContainerWithZeroData(self):
        fieldNames = ['depth']
        description = ''
        ta = timeArray.timeArray(np.array([0]), 'epoch')
        nodalValues = np.zeros((self.x.shape[0], 1))
        data = nodalValues[:, :, None]
        mc = meshContainer(
            description,
            ta,
            self.x,
            self.y,
            self.z,
            data,
            self.connectivity,
            fieldNames,
            coordSys='spcs')
        # add boundary information
        bnds = self.getBoundaries()
        if bnds is not None:
            for bnd in bnds:
                mc.addBoundary(bnd)
        return mc

    def writeNodalDataFromMC(self, outfilename, mc, iField=0, append=False,
                             forceVector=False, verbose=False):
        """Writes data from mesh container in GMSH format."""
        fieldName = mc.fieldNames[iField]
        nodalValues = mc.data[:, iField, 0]
        if not np.array_equal(mc.connectivity, self.connectivity):
            raise Exception(
                'Given MeshContainer not compatible with this mesh.')
        time = 0
        iteration = 0
        #values = nodalValues[mc.connectivity]
        # self.writeElemData(outfilename,values,fieldName,
        # time,iteration,append,forceVector,verbose)
        self.writeNodalData(outfilename, nodalValues[:, None], fieldName,
                            time, iteration, append, forceVector, verbose)

    def writeNodalData(self, outfilename, values, fieldName,
                       time, iteration,
                       append=False, forceVector=False, verbose=False):
        """Writes given element values to disk in GMSH field format."""
        if append:
            outfile = open(outfilename, 'a')
            print 'Appending data to', outfilename, '...'
        else:
            outfile = open(outfilename, 'w')
            print 'Writing data to', outfilename, '...'
            self.startSection(outfile, 'MeshFormat', verbose)
            outfile.write('2.2 0 8\n')
            self.endSection(outfile, 'MeshFormat', verbose)
        nNodes = values.shape[0]
        if self.lines is not None:
            nLines = self.lines.shape[0]
        else:
            nLines = 0
        # write header

        # NodeData
        if values.shape[0] > 0:
            nScalars = values.shape[1]
            valsIsMatrix = nScalars > 1
            #if forceVector : nScalars = 3
            # print 'nScalars',nScalars
            self.startSection(outfile, 'NodeData', verbose)
            outfile.write('{0:d}\n'.format(1))  # number of str tags
            outfile.write('\"{0:s}\"\n'.format(fieldName))  # field name
            outfile.write('{0:d}\n'.format(1))  # number of float tags
            outfile.write('{0:g}\n'.format(time))  # time step
            outfile.write('{0:d}\n'.format(3))  # number of int tags
            # iteration / export count
            outfile.write('{0:d}\n'.format(iteration))
            # number of fields per node
            outfile.write('{0:d}\n'.format(nScalars))
            # number of nodes with values
            outfile.write('{0:d}\n'.format(nNodes))
            num = np.arange(nNodes)[:, None] + 1
            tmp = np.hstack((num, values))
            fmt = ['%d']
            for i in range(nScalars):
                fmt.append('%.16g')
            np.savetxt(outfile, tmp, fmt=fmt)

            self.endSection(outfile, 'NodeData', verbose)
        outfile.close()

    def writeElemData(self, outfilename, values, fieldName,
                      time, iteration,
                      append=False, forceVector=False, verbose=False):
        """Writes given element values to disk in GMSH field format."""
        if append:
            outfile = open(outfilename, 'a')
            print 'Appending data to', outfilename, '...'
        else:
            outfile = open(outfilename, 'w')
            print 'Writing data to', outfilename, '...'
            self.startSection(outfile, 'MeshFormat', verbose)
            outfile.write('2.2 0 8\n')
            self.endSection(outfile, 'MeshFormat', verbose)
        nElems = values.shape[0]
        if self.lines is not None:
            nLines = self.lines.shape[0]
        else:
            nLines = 0
        # write header

        # ElementNodeData
        if values.shape[0] > 0:
            nNodes = self.connectivity.shape[1]
            nScalars = values.shape[1] / nNodes
            valsIsMatrix = nScalars > 1
            #if forceVector : nScalars = 3
            # print 'nScalars',nScalars
            self.startSection(outfile, 'ElementNodeData', verbose)
            outfile.write('{0:d}\n'.format(1))  # number of str tags
            outfile.write('\"{0:s}\"\n'.format(fieldName))  # field name
            outfile.write('{0:d}\n'.format(1))  # number of float tags
            outfile.write('{0:g}\n'.format(time))  # time step
            outfile.write('{0:d}\n'.format(3))  # number of int tags
            # iteration / export count
            outfile.write('{0:d}\n'.format(iteration))
            # number of fields per node
            outfile.write('{0:d}\n'.format(nScalars))
            # number of elements with values
            outfile.write('{0:d}\n'.format(nElems))
            num = np.arange(nElems)[:, None] + 1 + nLines
            tmp = np.hstack(
                (num, np.ones_like(num) * nNodes * nScalars, values))
            fmt = ['%d', '%d']
            for i in range(nNodes * nScalars):
                fmt.append('%.16g')
            np.savetxt(outfile, tmp, fmt=fmt)

            self.endSection(outfile, 'ElementNodeData', verbose)
        outfile.close()

    def readElemValues(self, filename):
        print 'reading element values from', filename, '...'
        infile = open(filename, 'r')
        while True:
            line = infile.readline().strip()
            if line == '':
                break  # end of file
            if line == '$ElementNodeData':
                break

        if self.lines is not None:
            nLines = self.lines.shape[0]
        else:
            nLines = 0
        nb = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb of str tags
        fieldName = infile.readline().strip().strip('\"')  # field name
        if nb > 1:
            for i in range(1, nb):
                infile.readline().strip().strip('\"')  # other tags
        nb = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb of float tags
        time = np.fromfile(
            infile,
            dtype=float,
            count=1,
            sep=' ')[0]  # time step
        nb = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb of int tags
        it = np.fromfile(infile, dtype=int, count=1, sep=' ')[0]  # iteration
        nScalars = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb scalars per node
        nTriangles = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb elems with values

        # print fieldName, time, it, nScalars, nTriangles

        vals = []  # all nodal values, each row for an element
        iElems = []  # all elements with values
        for i in range(nTriangles):
            ie, nbN = np.fromfile(infile, dtype=int, count=2, sep=' ')
            elemVals = np.fromfile(
                infile, dtype=float, count=3 * nScalars, sep=' ')
            vals.append(elemVals)
            iElems.append(ie)
        vals = np.array(vals)
        iElems = np.array(iElems, dtype=int) - 1 - nLines
        # iElems=np.arange(len(iElems)) # HACK assume all elements have values!
        nCols = vals.shape[1]
        elemValues = np.ones((self.connectivity.shape[0], nCols)) * np.nan
        elemValues[iElems, :] = vals
        return elemValues, fieldName, time, it, nScalars

    def readNodalValues(self, filename):
        print 'reading nodal values from', filename, '...'
        infile = open(filename, 'r')
        while True:
            line = infile.readline().strip()
            if line == '':  # end of file
                raise Exception(
                    '$NodeData section not found in file ' + filename)
            if line == '$NodeData':
                break

        nb = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb of str tags
        fieldName = infile.readline().strip().strip('\"')  # field name
        if nb > 1:
            for i in range(1, nb):
                infile.readline().strip().strip('\"')  # other tags
        nb = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb of float tags
        time = np.fromfile(
            infile,
            dtype=float,
            count=1,
            sep=' ')[0]  # time step
        nb = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb of int tags
        it = np.fromfile(infile, dtype=int, count=1, sep=' ')[0]  # iteration
        nScalars = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb scalars per node
        nNodes = np.fromfile(infile, dtype=int, count=1, sep=' ')[
            0]  # nb nodes with values

        # print fieldName, time, it, nScalars, nNodes

        # vals = [] # all nodal values, each row for an element
        # iElems = [] # all elements with values
        # for i in range(nTriangles) :
        #ie,nbN = np.fromfile(infile,dtype=int,count=2,sep=' ')
        #elemVals = np.fromfile(infile,dtype=float,count=3*nScalars,sep=' ')
        #vals.append( elemVals )
        #iElems.append( ie )
        # vals=np.array(vals)
        # iElems=np.array(iElems,dtype=int)-1-nLines
        #nCols = vals.shape[1]
        #elemValues = np.ones((self.connectivity.shape[0],nCols))*np.nan
        #elemValues[iElems,:] = vals
        vals = np.fromfile(infile, dtype=float, count=nNodes * 2, sep=' ')
        try:
            vals = np.reshape(vals, (nNodes, -1))
        except Exception as e:
            print e
            print vals.shape, nNodes
            raise Exception('Reading nodal values failed, too short data?')
        nodeNums = vals[:, 0].astype(int) - 1
        vals = vals[:, 1]
        maxNodes = nodeNums.max() + 1
        nodalValues = np.ones((maxNodes, 1)) * np.nan
        nodalValues[nodeNums, 0] = vals
        return nodalValues, fieldName, time, it, nScalars

    def writePosFileFromMC(self, outfilename, mc, iField=0, append=False):
        self.writePosFile(
            outfilename,
            mc.data[
                :,
                iField,
                0],
            mc.fieldNames[0],
            append)

    def writePosFile(
            self,
            outfilename,
            nodalValues,
            fieldName='field',
            append=False):
        """Writes data in meshContainer ti disk in GMSH POS format."""
        if append:
            outfile = open(outfilename, 'a')
            print 'Appending data to', outfilename, '...'
        else:
            outfile = open(outfilename, 'w')
            print 'Writing data to', outfilename, '...'

        elemVals = nodalValues[self.connectivity]
        X = self.x[self.connectivity]
        Y = self.y[self.connectivity]
        Z = self.z[self.connectivity]
        coords = np.vstack((X[:, 0], Y[:, 0], Z[:, 0],
                            X[:, 1], Y[:, 1], Z[:, 1],
                            X[:, 2], Y[:, 2], Z[:, 2])).T
        allData = np.hstack((coords, elemVals))
        # write header
        outfile.write('View \"%s\" {\n' % fieldName)  # field block
        nElems = self.connectivity.shape[0]
        # write triangles ST(x1,y1,z1,x2,y2,z2,x3,y3,z3){val1,val2,val3}
        fmt = 'ST(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n'
        for iElem in range(nElems):
            outfile.write(fmt % tuple(allData[iElem, :]))

        # write footer
        outfile.write('};\n')  # closes field block

    def writePosTensorFile(self, outfilename, M11, M12, M22,
                           fieldName='field', append=False):
        """Writes data in meshContainer ti disk in GMSH POS format."""
        if append:
            outfile = open(outfilename, 'a')
            print 'Appending data to', outfilename, '...'
        else:
            outfile = open(outfilename, 'w')
            print 'Writing data to', outfilename, '...'

        eM11 = M11[self.connectivity]
        eM12 = M12[self.connectivity]
        eM22 = M22[self.connectivity]
        tensor1 = np.vstack((eM11[:, 0], eM12[:, 0], eM12[:, 0], eM22[:, 0])).T
        tensor2 = np.vstack((eM11[:, 1], eM12[:, 1], eM12[:, 1], eM22[:, 1])).T
        tensor3 = np.vstack((eM11[:, 2], eM12[:, 2], eM12[:, 2], eM22[:, 2])).T
        X = self.x[self.connectivity]
        Y = self.y[self.connectivity]
        Z = self.z[self.connectivity]
        coords = np.vstack((X[:, 0], Y[:, 0], Z[:, 0],
                            X[:, 1], Y[:, 1], Z[:, 1],
                            X[:, 2], Y[:, 2], Z[:, 2])).T
        # write header
        outfile.write('View \"%s\" {\n' % fieldName)  # field block
        nElems = self.connectivity.shape[0]
        # write 3-by-3 tensor for each node in triangle
        # TT(x1,y1,z1,x2,y2,z2,x3,y3,z3){t1_11,t1_12,t1_13,t1_21,t1_22,...,t2_11,...}
        fmtHead = 'TT(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){'
        # write only 2-by-2 (x,y) part
        fmtNode = '%.16g,%.16g,0,%.16g,%.16g,0,0,0,1'
        fmtFoot = '};\n'
        for iElem in range(nElems):
            outfile.write(fmtHead % tuple(coords[iElem, :]))
            outfile.write(fmtNode % tuple(tensor1[iElem, :]))
            outfile.write(',')
            outfile.write(fmtNode % tuple(tensor2[iElem, :]))
            outfile.write(',')
            outfile.write(fmtNode % tuple(tensor3[iElem, :]))
            outfile.write(fmtFoot)

        # write footer
        outfile.write('};\n')  # closes field block


def findBoundarySegments(lines, lineTags, physicalNames):
    # build bnd tags for nodes
    tags = np.unique(lineTags)
    landNames = ['land']
    lineDim = 1
    landTag = [n for n in tags if physicalNames[(lineDim, n)][:4] in landNames]
    # find open bnd tags
    openTag = [
        n for n in tags if physicalNames[
            (lineDim, n)][
            :4] not in landNames]
    tags = list(openTag)
    tags.extend(landTag)

    class lineObj(object):

        def __init__(self, i, nodes):
            self.i = i
            self.nodes = nodes

    from collections import deque

    def sortLines(lines):
        segments = []  # list of deques representing continuous lines
        # sort function for lines
        for i, lineNodes in enumerate(lines):
            L = lineObj(i, list(lineNodes))
            # find if connects to any existing seqments
            parentSegment = None
            for s in segments:
                if s[0].nodes[0] in L.nodes:
                    # add in beginning
                    if s[0].nodes[0] == L.nodes[0]:
                        L.nodes = L.nodes[::-1]  # reverse
                    s.appendleft(L)
                    parentSegment = s
                    break
                elif s[-1].nodes[-1] in L.nodes:
                    # add in end
                    if s[-1].nodes[-1] == L.nodes[-1]:
                        L.nodes = L.nodes[::-1]  # reverse
                    s.append(L)
                    parentSegment = s
                    break
            if parentSegment is not None:
                # check for merged segments
                nodes = [
                    parentSegment[0].nodes[0],
                    parentSegment[-1].nodes[-1]]
                merged = False
                for s in segments:
                    if s == parentSegment:
                        continue
                    if s[0].nodes[0] in nodes:
                        # add in beginning
                        seg = deque(parentSegment)
                        if s[0].nodes[0] == nodes[0]:
                            seg.reverse()  # reverse
                        s.extendleft(seg)
                        merged = True
                        break
                    elif s[-1].nodes[-1] in nodes:
                        # add in end
                        seg = deque(parentSegment)
                        if s[-1].nodes[-1] == nodes[-1]:
                            seg.reverse()  # reverse
                        s.extend(seg)
                        merged = True
                        break
                if merged:
                    segments.remove(parentSegment)
            else:
                # add as a new segment
                segments.append(deque([L]))
        lineSegs = []
        for s in segments:
            line = np.array([L.nodes for L in s])
            lineSegs.append(line)
        return lineSegs

    boundaries = []
    #processedNodes = set()
    for tag in tags:
        taggedLines = lines[lineTags == tag, :]
        # need to find connected line segments
        segments = sortLines(taggedLines)
        tagStr = physicalNames[(1, tag)]
        bndType = 'land'
        if 'island' in tagStr:
            bndType = 'island'
        if 'open' in tagStr:
            bndType = 'open'
        # remove "bndType" from tagStr
        tagStr = tagStr.replace(bndType, '').strip()
        # print 'bnd tag: ',tag, tagStr
        for segment in segments:
            # convert list of lines to list of nodes
            tmp = np.array(list(segment))
            nodes = np.hstack((tmp[:, 0], [tmp[-1, 1]]))
            # remove all nodes already added
            # if nodes[0] in processedNodes :
            #nodes = nodes[1:]
            # if nodes[-1] in processedNodes :
            #nodes = nodes[:-1]
            boundaries.append(meshBoundary(bndType, tagStr, nodes))
            #processedNodes.update( nodes )

    return boundaries
