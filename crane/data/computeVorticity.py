"""
Computes the relative vorticity from SELFE output files.

Relative vorticity is defined as a curl of the horizontal velocity field
zeta = curl_z(u) = dv/dx - du/dy

This is typically much smaller than the Coriolis parameter f. Typical
non-dimensional form is
zeta^* = zeta/f

Absolute vorticity is defined as
zeta_a = f + zeta

Potential vorticity is a conserved quantity. In the case of barotropic flow
it becomes
Pi = (zeta + f)/H
where H is the depth of the water column. In the case of stratified
(baroclinic) flows
Pi = (zeta + f)/rho drho/dz
where rho is the water density.

Here we only compute the relative vorticity. Making use of the Stokes' theorem

zeta = 1/A Int_e (dv/dx - du/dy) dx = 1/A Int_de (u . dr)
where dr is the counter-clockwise surface vector of surface de of element e,
and A is the area of the element.

For a linear non-conforming velocity field (hvel.70) the latter integral can
be evaluated as
zeta = 1/A sum(u_i . r_i)
where u_i is the horizontal velocity vector at edges and r_i is the
counter-clockwise oriented edge vector.

Tuomas Karna 2016-03-15
"""
import numpy
from crane import *
from crane.data.ncExtract import selfeNCFile, selfeExtractBase, splitVarToFileType
from crane.data.meshContainer import computeAreas
import netCDF4 as nc

MISSINGVAL = numpy.float32(-9999.0)


def create_new_nc_file(old, replace_dict, new_path, verbose=False):
    """
    Create a netcdf file using another as a template but with new data.

    """
    new = nc.Dataset(new_path, 'w', format='NETCDF4')
    # Copy dimensions
    for dim_name, dim in old.dimensions.iteritems():
        if verbose:
            print ' copying dim', dim_name, len(dim)
        dim_len = len(dim)
        if dim.isunlimited() or dim_name in ['time']:
            dim_len = None
        new.createDimension(dim_name, dim_len)

    # Copy variables with skip/replace
    bookeeping_dict = dict(replace_dict)
    for name, var in old.variables.iteritems():
        if name not in bookeeping_dict:
            # copy as it is
            out_var = new.createVariable(name, var.dtype, var.dimensions)
            out_var[:] = var[:]
            if verbose:
                print ' copying', name, var.dimensions
        elif name in bookeeping_dict and bookeeping_dict[name] is None:
            # skip variable
            bookeeping_dict.pop(name)
            if verbose:
                print ' skipping', name
        else:
            # replace variable
            data, dtype, dims, meta = bookeeping_dict.pop(name)
            out_var = new.createVariable(name, dtype, dims)
            out_var[:] = data[:]
            if meta is not None:
                raise NotImplementedError('replacing meta data not implemented')
            if verbose:
                print ' replacing', name, var.dimensions, dims
        # Copy variable attributes
        for name in vars(var):
            # Skip private
            if not name.startswith('_'):
                out_var.setncattr(name, var.getncattr(name))
    for name in bookeeping_dict:
        data, dtype, dims, meta = bookeeping_dict[name]
        out_var = new.createVariable(name, dtype, dims)
        if data is not None:
            out_var[:] = data[:]
        if verbose:
            print ' adding', name, dims
        if meta is not None:
            for attrname in meta:
                if not attrname.startswith('_'):
                    out_var.setncattr(attrname, meta[attrname])

    new.sync()
    return new


def process_file(filename, outputdir=None, output_at_element=False,
                 verbose=False):
    """
    Opens hvel.xx velocity file and computes relative vorticity for the
    entire mesh and time steps. Output is stored either as a vort.70 or
    vort.63 file that mimics SELFE netcdf output convention.
    """
    sfile = selfeNCFile(filename, verbose=verbose)
    sfile.readHeader()

    nelem = sfile.nFaces
    nnode = sfile.nNodes
    nvert = sfile.nVert
    ntime = sfile.nTime

    assert sfile.discrType in ['edge', 'node'], 'only 64 and 70 velocity files are supported'
    uv_at_edges = sfile.discrType == 'edge'

    u_var = sfile.variables['hvel_x']
    v_var = sfile.variables['hvel_y']

    conn = sfile.faceNodes
    nodex = sfile.nodeX
    nodey = sfile.nodeY
    areas = computeAreas(conn, nodex, nodey)

    elem_x = nodex[conn]  # (nelem, 3)
    elem_y = nodey[conn]  # (nelem, 3)

    centroidx = numpy.mean(elem_x, axis=1)  # (nelem, )
    centroidy = numpy.mean(elem_y, axis=1)  # (nelem, )

    edge_vec_x = numpy.diff(numpy.hstack((elem_x, elem_x[:, [0]])))  # (nelem, 3)
    edge_vec_y = numpy.diff(numpy.hstack((elem_y, elem_y[:, [0]])))  # (nelem, 3)

    if uv_at_edges:
        # invert edge_node map
        edgenodes = sfile.edgeNodes
        nodes_to_edge = {}  # maps (node1, node2) to iedgenode
        for i in range(edgenodes.shape[0]):
            key = tuple(sorted(edgenodes[i, :]))
            nodes_to_edge[key] = i

        # maps (ielem, ielemedge) to iedgenode
        elem_edge_to_edge = numpy.zeros((nelem, 3), dtype=int)  # (nelem, 3)
        for ielem in range(nelem):
            for i in range(3):
                key = tuple(sorted([conn[ielem, i], conn[ielem, (i + 1) % 3]]))
                elem_edge_to_edge[ielem, i] = nodes_to_edge[key]

    # dict of fields to replace or skip, all other vars are copied
    # key: (data_array, dtype, dims, dict_of_metadata)
    replace = {'hvel_x': None,  # skip
               'hvel_y': None,  # skip
               }
    # get elev, depth, k_bottom at nodes
    if uv_at_edges:
        # get elev, depth and k_bottom files from elev.61 file
        # NOTE these cannot be retrieved from hvel.67 file
        elevfilename = filename.replace('hvel.67', 'elev.61')
        elevfile = selfeNCFile(elevfilename, verbose=verbose)
        elev_node = elevfile.variables['elev'][:]
        depth_node = elevfile.variables['depth'][:]
        k_bottom_node = elevfile.variables['k_bottom'][:]
    else:
        # convert values at nodes to centroids
        elev_node = sfile.variables['elev'][:]
        depth_node = sfile.variables['depth'][:]
        k_bottom_node = sfile.variables['k_bottom'][:]

    # convert to correct output arrays
    if output_at_element:
        # generate correct elev,depth and k_bottom arrays for .70 output file
        if uv_at_edges:
            # convert values at edges to centroids
            elev_out = sfile.variables['elev'][:][:, elem_edge_to_edge].mean(axis=2)
            depth_out = sfile.variables['depth'][:][elem_edge_to_edge].mean(axis=1)
            k_bottom_out = sfile.variables['k_bottom'][:][elem_edge_to_edge].min(axis=1)
        else:
            # convert values at nodes to centroids
            elev_out = sfile.variables['elev'][:][:, conn].mean(axis=2)
            depth_out = sfile.variables['depth'][:][conn].mean(axis=1)
            k_bottom_out = sfile.variables['k_bottom'][:][conn].min(axis=1)
        replace['elev'] = (elev_out, numpy.float32, ('time', 'face'), None)
        replace['depth'] = (depth_out, numpy.float32, ('face',), None)
        replace['k_bottom'] = (k_bottom_out, numpy.int32, ('face',), None)
    else:
        if uv_at_edges:
            elev_out = elev_node
            depth_out = depth_node
            k_bottom_out = k_bottom_node
            replace['elev'] = (elev_out, numpy.float32, ('time', 'node'), None)
            replace['depth'] = (depth_out, numpy.float32, ('node',), None)
            replace['k_bottom'] = (k_bottom_out, numpy.int32, ('node',), None)
        else:
            # elev, depth and k_bottom already at nodes
            pass

    # add new variable with no data
    meta = {'missing_value': MISSINGVAL,
            'standard_name': 'ocean_relative_vorticity',
            'long_name': 'ocean_relative_vorticity',
            'units': 's-1',
            }
    nodedim = 'face' if output_at_element else 'node'
    replace['vort'] = (None, numpy.float32, ('time', 'layers', nodedim), meta)

    # generate output file name
    fext = sfile.fileTypeStr
    inputdir, inputfile = os.path.split(filename)
    if outputdir is None:
        outputdir = inputdir
    outext = '70' if output_at_element else '63'
    outputfile = inputfile.replace('hvel.' + fext, 'vort.' + outext)
    outputfile = os.path.join(createDirectory(outputdir), outputfile)

    # open a new nc file by copying most fields from the input file
    print('Saving to {:}'.format(outputfile))
    out_nc = create_new_nc_file(sfile.ncfile, replace, outputfile,
                                verbose=verbose)

    # compute vorticity
    # vorticity field for one time step
    zeta = numpy.zeros((nvert, nelem), dtype=numpy.float32)
    if not output_at_element:
        zeta_node = numpy.zeros((nvert, nnode), dtype=numpy.float32)

    for itime in range(ntime):
        u = u_var[itime, :, :]  # shape (ntime, nvert, nnodes)
        v = v_var[itime, :, :]  # shape (ntime, nvert, nnodes)

        # compute dry elements
        dry_elems, dry_nodes = sfile.vCoords.computeDryElemMask(elev_node[itime, :], depth_node, sfile.faceNodes)

        # mask bad values
        u[u < -99.0] = 0.0
        v[v < -99.0] = 0.0
        if uv_at_edges:
            dry_edges = dry_nodes[edgenodes].max(axis=1)
            u[:, dry_edges] = 0.0
            v[:, dry_edges] = 0.0
        else:
            u[:, dry_nodes] = 0.0
            v[:, dry_nodes] = 0.0

        # compute vorticity with Stokes theorem
        zeta[:] = 0.0
        for i in range(3):
            # compute dot(uv, r)/area for edge i in all triangles
            sign = (-edge_vec_x[:, i]*(centroidy - elem_y[:, i]) +
                    edge_vec_y[:, i]*(centroidx - elem_x[:, i]))
            assert (sign < 0).all(), 'edge oriented in wrong direction'
            if uv_at_edges:
                j = elem_edge_to_edge[:, i]
                u_edge = u[:, j]
                v_edge = v[:, j]
            else:
                n1 = conn[:, i]
                n2 = conn[:, (i + 1) % 3]
                u_edge = 0.5*(u[:, n1] + u[:, n2])
                v_edge = 0.5*(v[:, n1] + v[:, n2])
            zeta[:, :] += (edge_vec_x[:, i]*u_edge +
                           edge_vec_y[:, i]*v_edge)/areas

        # compute mean over elements, shape (ntime, nvert, nelem)
        zeta[:-1, :] = 0.5*(zeta[:-1, :] + zeta[1:, :])
        # NOTE last value is not used, -2 is the surface for half levels
        zeta[-1, :] = zeta[-2, :]
        bad_values = np.logical_or(~numpy.isfinite(zeta),
                                   numpy.abs(zeta) > 99.0)
        zeta[bad_values] = numpy.nan
        zeta[:, dry_elems] = numpy.nan

        if output_at_element:
            # append to netcdf file
            zeta[numpy.isnan(zeta)] = MISSINGVAL
            out_nc.variables['vort'][itime, :, :] = zeta
        else:
            # convert element-wise zeta to nodes
            node2elem = sfile.meshSearch2d.node2elem
            zeta_node[:] = 0.0
            for inode in xrange(nnode):
                elemlist = node2elem[inode]
                elemlist = elemlist[numpy.logical_not(dry_elems[elemlist])]
                # average to nodes and half levels
                vals = zeta[:, elemlist]
                multiplicity = len(elemlist)
                if multiplicity == 0:
                    zeta_node[:, inode] = numpy.nan
                    continue
                zeta_node[:, inode] = numpy.nansum(vals, axis=1)/multiplicity
                # average to full levels
                zeta_node[-1, inode] = zeta_node[-2, inode]  # surface
                zeta_node[k_bottom_node[inode]+1:-2, inode] = 0.5*(zeta_node[k_bottom_node[inode]:-3, inode] +
                                                                   zeta_node[k_bottom_node[inode]+1:-2, inode])
                zeta_node[k_bottom_node[inode]-1, inode] = zeta_node[k_bottom_node[inode], inode]  # bottom
            zeta_node[numpy.isnan(zeta_node)] = MISSINGVAL
            out_nc.variables['vort'][itime, :, :] = zeta_node

    out_nc.close()


def process(datadir, outputdir, var, stacks=None, startdate=None,
            enddate=None, output_at_element=False):
    """
    Extracts vorticity for the given stacks
    """
    varstr, filetypestr = splitVarToFileType(var)
    seb = selfeExtractBase(datadir, var=varstr, fileTypeStr=filetypestr)

    if stacks is None:
        stacks = seb.dataFile.getStacks(startdate, enddate, wholeDays=True)

    for stack in stacks:
        ncfilename = seb.generateFileName(iStack=stack)
        process_file(ncfilename, outputdir, output_at_element=output_at_element)


def parse_commandline():
    """Parse command line arguments"""
    import argparse
    parser = argparse.ArgumentParser(description='Computes relative vorticity from SELFE output files. Output is a netcdf file akin to SELFE output files.')
    parser.add_argument('-d', dest='datadir', help='path to combined data directory', required=True)
    parser.add_argument('-o', dest='outputdir', help='path where output will be stored (default: combined data dir)')
    parser.add_argument('-v', dest='var', choices=['hvel.64', 'hvel.67'], help='variable to operate on (default hvel.67)', default='hvel.67')
    parser.add_argument('-S', dest='stackStr', type=str, help='range of output files to read '
                        '(e.g 1,14) if start,end not given')
    parser.add_argument('-s', dest='startdate', type=str, help='First date to process (e.g. 2012-5-1)')
    parser.add_argument('-e', dest='enddate', type=str, help='Last date to process (e.g. 2012-5-1)')
    parser.add_argument('--output-at-element', action='store_true', help='Store vorticity at each 3d element instead of nodes (default False)')

    args = parser.parse_args()

    startdate = enddate = stacks = None
    if args.stackStr is None:
        assert args.startdate is not None, 'startdate must be given'
        assert args.enddate is not None, 'enddate must be given'
        startdate = parseTimeStr(args.startdate)
        enddate = parseTimeStr(args.enddate)
    else:
        first, last = [int(p) for p in args.stackStr.split(',')]
        stacks = range(first, last + 1)

    print('Parsed options:')
    print(' - datadir: {:}'.format(args.datadir))
    if args.outputdir:
        print(' - output dir: {:}'.format(args.outputdir))
    print(' - variable: {:}'.format(args.var))
    if stacks is None:
        print(' - time range: {:} -> {:}'.format(str(startdate), str(enddate)))
    else:
        print(' - stacks: {:} -> {:}'.format(str(stacks[0]), str(stacks[-1])))
    print(' - output at element: {:}'.format(args.output_at_element))
    sys.stdout.flush()

    process(args.datadir, args.outputdir, args.var, stacks=stacks,
            startdate=startdate, enddate=enddate, output_at_element=args.output_at_element)

if __name__ == '__main__':
    parse_commandline()
