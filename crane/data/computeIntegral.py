"""
Computes spatial and temporal averages of SELFE scalar fields.

Tuomas Karna 2016-06-15
"""
from crane.utility import parseTimeStr
from crane.data.ncExtract import *
import crane.data.meshContainer
from crane.physicalVariableDefs import addTracers, fieldNameToFilename, fieldNameList
import crane.data.dirTreeManager as dtm
from crane.data.timeSeriesFilters import T_M2
from crane.utility import readAnyMeshFile


def slab_area_mean(dc, outvar, regionmc, regionname, regionid, integral=False):
    """
    Averages a slab dataContainer in an area defined in regionmc.

    Region to average over is designated by nodal value regionid in regionmc.
    region_name is a descriptive name for the region.

    If integral==True, will compute horizontal integral instead.
    """


    # compute mean over region
    areas = dc.computeAreas()
    if regionmc:
        elem_mask = np.mean(regionmc.data[regionmc.connectivity, 0, 0], axis=1) == regionid
        elem_mask = np.nonzero(elem_mask)[0]
        elem_data = np.mean(dc.data[regionmc.connectivity[elem_mask, :], :, :], axis=1)
        w = areas[elem_mask][:, np.newaxis, np.newaxis]
        name = regionname
    else:
        # integrate over the whole domain
        elem_data = np.mean(dc.data[dc.connectivity[:, :], :, :], axis=1)
        w = areas[:, np.newaxis, np.newaxis]
        name = 'domain'
    if integral:
        area_int = (np.nansum(elem_data*w, axis=0))
    else:
        area_int = (np.nansum(elem_data*w, axis=0)/np.nansum(w, axis=0))

    ta = dc.time
    if outvar is None:
        var = dc.getMetaData('variable')
    else:
        var = outvar
    fieldnames = dc.fieldNames
    meta = {}
    meta['location'] = name
    meta['variable'] = var
    meta['bracket'] = 'A'
    meta['msldepth'] = '0'
    meta['dataType'] = 'timeseries'
    meta['tag'] = dc.getMetaData('tag')
    data = area_int[np.newaxis, :, :]
    dc = dataContainer.dataContainer('', ta, 0, 0, 0, data,
                                     fieldnames, metaData=meta,
                                     coordSys='spcs')
    return dc


def slab_horizontal_integral(slab_dc, var=None, regionfile=None,
                            regionname=None, regionid=1,
                            mode='average'):
    """
    Computes a horizontal average of a slab meshContainer.

    If regionfile (.gr3 file) is defined, uses integrates over subregion defined by regionid and name regionname.
    Otherwise integrates over the whole domain.

    Returns a time series in a dataContainer.
    """
    if regionfile is not None:
        regionmc = readAnyMeshFile(regionfile)
        if regionname is None:
            raise Exception('region name is required')
    else:
        regionmc = None
    integral = mode == 'integral'
    mc = slab_area_mean(slab_dc, var, regionmc, regionname, regionid,
                        integral=integral)
    return mc

def read_3d_mesh_data(datadir, var, stack, meshsearchobj=None,
                      comp_z_coords=False):
    """
    Reads data from SELFE netcdf file to arrays.
    """
    varstr = getNCVariableName(var)
    extractor = selfeExtractBase(datadir, varstr,
                                 meshSearchObj=meshsearchobj, verbose=True)
    if meshsearchobj is None:
        meshsearchobj = extractor.dataFile.meshSearch2d
    filename = extractor.generateFileName(stack)
    if not os.path.isfile(filename):
        raise IOError('File not found: '+filename)
    ncfile = extractor.getNCFile(stack)
    if varstr in vectorVars:
        # recursion: if a vector field requested, extract components
        # separately
        varlist = [varstr + '_x', varstr + '_y']  # hvel_x, hvel_y
    else:
        varlist = [varstr]

    data = {}
    for var in varlist:
        d = ncfile.variables[var][:].astype(np.float32)
        data[var] = np.ma.masked_invalid(d, copy=False)

    time = ncfile.getTime()
    facenodes = meshsearchobj.faceNodes
    x = meshsearchobj.node_x
    y = meshsearchobj.node_y
    # compute and assign dry mask
    elev = ncfile.variables['elev'][:].astype(np.float32)
    depth = ncfile.variables['depth'][:].astype(np.float32)
    nTime, nPoints = elev.shape
    nvert = extractor.dataFile.nVert
    drynodes = np.zeros_like(elev, dtype=bool)
    vcoords = extractor.dataFile.vCoords
    for i in xrange(nTime):
        dryE, dryN = vcoords.computeDryElemMask(elev[i, :], depth, facenodes)
        drynodes[i, :] = dryN
    if comp_z_coords:
        z = np.zeros((nTime, nvert, nPoints))
        for i in xrange(nTime):
            # zi shape is (nVert, nPoints)
            zi, kbp2, iwet = vcoords.computeVerticalCoordinates(elev[i, :], depth)
            z[i, :, :] = zi
    else:
        z = None
    # mask out bad data values
    print 'dry nodes', np.sum(drynodes)
    for var in varlist:
        d = data[var]
        d[d < -900] = np.nan
        d.mask[:, :, :] = drynodes[:, np.newaxis, :]
        d[d.mask] = np.nan
    return data, time, x, y, z, facenodes, meshsearchobj


def slab_array_to_meshcontainer(data, time, x, y, faceNodes, varstr, fieldnames, runtag):
    """
    Store 2D mesh arrays in a meshContainer
    """
    data = data.T[:, np.newaxis, :]
    ta = timeArray.timeArray(time, 'epoch')
    z = np.zeros_like(x)
    meta = {}
    meta['dataType'] = 'slab'
    meta['location'] = 'slab'
    meta['instrument'] = 'model'
    meta['variable'] = varstr
    meta['slevel'] = '0'
    meta['tag'] = runtag
    fieldnamelist = fieldnames
    mc = meshContainer.meshContainer('', ta, x, y, z, data, faceNodes, fieldnamelist,
                                     coordSys='spcs', metaData=meta,)
    return mc


def get_stacks_for_time_range(datadir, var, starttime, endtime):
        varstr = getNCVariableName(var)
        extractor = selfeExtractBase(datadir, varstr,
                                     verbose=False)
        stacks = extractor.dataFile.getStacks(starttime, endtime)
        print ' - stacks:', stacks
        return stacks


def compute_depth_average(dz, arr, depth):
    """Computes depth average of an array whose shape is (:, nVert ,:)"""
    b = 0.5*(arr[:, 1:, :] + arr[:, :-1, :])
    b = b*dz
    b = np.nansum(b, axis=1)/depth
    # copy nan mask from the surface layer
    mask = ~np.isfinite(arr[:, -1, :])
    b[mask] = np.nan
    return b


def get_depth_averaged_mc(data_dict, var, time, x, y, depth, dz,
                          facenodes, runtag, depth_integral=False):
    """
    Computes depth averages from data array and returns a meshContainer.

    Handles vector valued fields correctly.
    """
    mc = None
    depth_scalar = 1.0 if depth_integral else depth
    varprefix = 'dint' if depth_integral else 'dav'
    for comp in sorted(data_dict.keys()):
        davarr = compute_depth_average(dz, data_dict[comp], depth_scalar)
        comp_mc = slab_array_to_meshcontainer(davarr, time, x, y, facenodes,
                                     varprefix+var, [varprefix + comp], runtag)
        comp_mc.setMetaData('original-variable', var)
        if mc is None:
            mc = comp_mc
        else:
            mc.mergeFields(comp_mc)
    return mc


def extract_depth_average(datadir, runtag, varlist, stacks=None,
                          starttime=None, endtime=None, depth_integral=False):
    """
    Extracts a depth averaged (or integrated) field from SELFE netcdf outputs
    """
    # NOTE move this guy to ncExtract module
    if stacks is None:
        tmp_var = varlist[0]
        stacks = get_stacks_for_time_range(datadir, tmp_var, starttime, endtime)

    meshsearchobj = None

    def process_stack(stack, meshsearchobj):
        # array for z coordinates
        dz = None

        # stores all meshContainers
        dcs = []

        for var in varlist:
            ddict, time, x, y, z, facenodes, meshsearchobj = \
                            read_3d_mesh_data(datadir, var, stack,
                                              meshsearchobj=meshsearchobj,
                                              comp_z_coords=(dz is None))
            if dz is None:
                # compute z coords
                dz = np.diff(z, axis=1)
                depth = np.nanmax(z, axis=1) - np.nanmin(z, axis=1)
            # compute depth average
            mc = get_depth_averaged_mc(ddict, var, time, x, y, depth, dz, facenodes, runtag,
                                       depth_integral=depth_integral)
            dcs.append(mc)

        return dcs

    concatenated_dcs = []
    for stack in stacks:
        try:
            dcs = process_stack(stack, meshsearchobj)
            if len(concatenated_dcs) == 0:
                concatenated_dcs.extend(dcs)
            else:
                # merge times
                for i in range(len(dcs)):
                    concatenated_dcs[i].mergeTemporal(dcs[i])
        except Exception as e:
            print 'Extraction failed.'
            traceback.print_exc(file=sys.stdout)
            print e

    return concatenated_dcs


def process(runtag, datadir, varlist, starttime=None, endtime=None,
            stacks=None, save_monthly_file=False,
            vertical_integral='depth_average',
            horizontal_integral=None,
            time_average=None,
            slab_z=None, slab_k=None, z_rel_to_surf=False,
            regionfile=None, regionid=None, regionname=None):
    file_rule = 'monthlyFile' if save_monthly_file else 'singleFile'

    # --------------------------------------------------------------------------
    # get 2D data in a meshContainer for depth average/integral/slab
    # --------------------------------------------------------------------------
    if vertical_integral in ['depth_average', 'depth_integral']:
        depth_integral = vertical_integral == 'depth_integral'
        print('Extracting depth average')
        dcs = extract_depth_average(datadir, runtag, varlist,
                                    starttime=starttime, endtime=endtime,
                                    stacks=stacks,
                                    depth_integral=depth_integral)
    else:
        print('Extracting slab')
        dcs = extractSlabForLevel(datadir, varlist, starttime, endtime, 'slab',
                                  z=slab_z, k=slab_k,
                                  zRelToSurf=z_rel_to_surf,
                                  wholeDays=True, stacks=stacks,
                                  verbose=False)
        for dc in dcs:
            dc.setMetaData('tag', runtag)

    for dc in dcs:
        print dc

    dtm.saveDataContainerInTree(dcs, rule=file_rule, dtype=np.float32,
                                overwrite=True, compress=True)

    # --------------------------------------------------------------------------
    # Optionally compute horisontal integral
    # --------------------------------------------------------------------------
    if horizontal_integral is not None:
        hor_domain, hor_mode = horizontal_integral.split('_')
        rfile = regionfile if hor_domain == 'region' else None
        new_dcs = []
        for dc in dcs:
            orig_var = dc.getMetaData().get('original-variable', dc.getMetaData('variable'))
            new_dc = slab_horizontal_integral(dc, var=orig_var, regionfile=rfile,
                                             regionname=regionname,
                                             regionid=regionid,
                                             mode=hor_mode)
            new_dcs.append(new_dc)
        dcs = new_dcs

        for dc in dcs:
            print dc
        dtm.saveDataContainerInTree(dcs, rule=file_rule, dtype=np.float32,
                                    overwrite=True, compress=True)

    # --------------------------------------------------------------------------
    # Optionally compute time average
    # --------------------------------------------------------------------------
    if time_average is not None:
        if time_average == 'lowpass':
            # apply removeTides filter
            raise NotImplementedError('Temporal lowpass filter has not been implemented yet')
        elif time_average == 'mean':
            # apply box car filter that handles NaNs correctly
            if starttime is None:
                starttime = dcs[0].time.getDatetime(0)
            if endtime is None:
                endtime = dcs[0].time.getDatetime(-1)
            print('Computing time average over period:\n  {:} -> {:}'.format(starttime, endtime))
            for dc in dcs:
                # ix = dc.time.getRangeIndices(st, et)
                ix = range(len(dc.time))
                dc.time.array = np.array([np.mean(dc.time.array[ix])])
                dc.data = np.nanmean(dc.data[:, :, ix], axis=2)[:, :, np.newaxis]
                var = dc.getMetaData('variable')
                new_var = 'tav' + var
                dc.setMetaData('variable', new_var)
                dc.fieldNames = ['tav' + f for f in dc.fieldNames]

            dtm.saveDataContainerInTree(dcs, rule=file_rule, dtype=np.float32,
                                        overwrite=True, compress=True)

# -----------------------
# command line interface
# -----------------------
def parse_commandline():
    import argparse

    parser = argparse.ArgumentParser(description='Extracts a vertical integral of a scalar field. Optionally also applies horisontal integral over given region or whole domain, and/or time averaging.')
    parser.add_argument('vertical_integral', choices=['slab', 'depth_average', 'depth_integral'], default='depth_average',
                        help='Vertical integration mode. If slab, a slab is extracted with given k level or depth z')
    parser.add_argument('-r', dest='runtag', required=True,
                        help='Run tag, used as a label in post-proc.')
    parser.add_argument('-d', dest='datadir',  required=True,
                        help='directory where model outputs are stored')
    parser.add_argument('-s', dest='startstr',
                        help='Date to start processing')
    parser.add_argument('-e', dest='endstr',
                        help='Date to end processing. If omitted, start time plustidal day is used.')
    parser.add_argument('--stacks',
                        dest='stacklim', help='range of output files to read (e.g 1,14) if start,end not given')
    parser.add_argument('-v', dest='varlist', required=True,
                        help='variable(s) to extract: elev,temp,salt, ... To use a specific output file define extrension, e.g. salt.70')
    parser.add_argument('--save-in-tree', action='store_true', dest='save_monthly_file',
                        help='saves extracted data in file tree with monthly files instead of a single file', default=False)
    parser.add_argument('-z', '--zCoord', type=float, dest='z_coord',
                        help='Vertical position of slab to extract: z coordinate (z<0 below datum)')
    parser.add_argument('-k', '--kLevel', type=int,
                        dest='k_level', help='vertical position of slab to extract: S level, 1 means bottom level, -1 surface')
    parser.add_argument('-S', '--zRelToSurf', action='store_true',
                        dest='z_rel_to_surf', help='Slab z coordinate (z>0) is defined as depth below free surface (default False)')
    parser.add_argument('--horizontal-integral', help='Compute integral over horizontal plane as well',
                        dest='horizontal_integral', choices=['domain_integral', 'domain_average', 'region_integral', 'region_average'])
    parser.add_argument('--region', help='GR3 region file for defining horizontal sub-regions',
                        dest='regionfile')
    parser.add_argument('--region-id', help='Integer value in the region file that defines the desired sub-region (default 1)',
                        dest='regionid', type=int, default=1)
    parser.add_argument('--region-name', help='Human readable name for the subregion',
                        dest='regionname')
    parser.add_argument('--time-average', help='Compute temporal average as well',
                        dest='time_average', choices=['lowpass', 'mean'])


    args = parser.parse_args()

    if args.stacklim is None and args.startstr is None:
        parser.error('Either -s or --stacks must be given')
    if args.stacklim is not None:
        a, b = [int(s) for s in args.stacklim.split(',')]
        stacks = np.arange(a, b + 1)
        st = et = None
    else:
        stacks = None
        if args.startstr is not None:
            st = parseTimeStr(args.startstr)
        if args.endstr is not None:
            et = parseTimeStr(args.endstr)
        else:
            et = st + datetime.timedelta(seconds=2*T_M2)

    if args.horizontal_integral == 'region':
        if args.regionfile is None:
            parser.error('region file is required')
        if args.regionname is None:
            parser.error('region name is required')

    varlist = args.varlist.split(',')
    print('runtag        : {:}'.format(args.runtag))
    print('data dir      : {:}'.format(args.datadir))
    if st is not None:
        print('start time    : {:}'.format(st))
        print('end time      : {:}'.format(et))
    else:
        print('stacks        : {:}'.format(stacks))
    print('variables     : {:}'.format(varlist))
    print('vertical int : {:}'.format(args.vertical_integral))
    if args.horizontal_integral is not None:
        print('horisontal int : {:}'.format(args.horizontal_integral))
    if args.horizontal_integral == 'region':
        print('region file    : {:}'.format(args.regionfile))
        print('region name    : {:}'.format(args.regionname))
        print('region id      : {:}'.format(args.regionid))
    if args.horizontal_integral == 'slab':
        if args.k_level:
            print('slab k level    : {:}'.format(args.k_level))
        elif args.z_coord:
            print('slab z coord    : {:}'.format(args.k_level))
            print('z relative to surf : {:}'.format(args.z_rel_to_surf))
    if args.time_average is not None:
        print('time average : {:}'.format(args.time_average))
    print('save-in-tree  : {:}'.format(args.save_monthly_file))

    process(args.runtag, args.datadir, varlist, starttime=st, endtime=et,
            stacks=stacks,
            save_monthly_file=args.save_monthly_file,
            vertical_integral=args.vertical_integral,
            horizontal_integral=args.horizontal_integral,
            time_average=args.time_average,
            slab_z=args.z_coord, slab_k=args.k_level,
            z_rel_to_surf=args.z_rel_to_surf,
            regionfile=args.regionfile,
            regionname=args.regionname,
            regionid=args.regionid)

if __name__ == '__main__':
    parse_commandline()
