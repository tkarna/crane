"""
Integrates scalar fields over a domain that is defined by min/max values of another scalar field.

Tuomas Karna 2016-10-11
"""
from crane import *
from crane.utility import parseTimeStr
from crane.data.ncExtract import *
import crane.data.meshContainer
from crane.data.meshContainer import computeAreas
from crane.physicalVariableDefs import addTracers, fieldNameToFilename, fieldNameList
import crane.data.dirTreeManager as dtm
from crane.data.timeSeriesFilters import T_M2
from crane.utility import readAnyMeshFile

import warnings


class SelfeOutputReader(object):
    """
    Simplified wrapper to read 3D SELFE output files into numpy arrays
    """
    def __init__(self, datadir, varlist):
        self.datadir = datadir
        self.varlist = varlist

        self.extractors = {}
        self.meshsearchobj = None
        self._initialize()

    def _get_extractor(self, var):
        v, filetypestr = splitVarToFileType(var)
        varstr = getNCVariableName(v)
        extractor = selfeExtractBase(self.datadir, varstr, fileTypeStr=filetypestr,
                                     meshSearchObj=self.meshsearchobj, verbose=True)

        if self.meshsearchobj is None:
            self.meshsearchobj = extractor.dataFile.meshSearch2d

        msg = 'currently only nodal "*.63" or "*.72" files are supported'
        assert extractor.dataFile.discrType == 'node', msg
        return v, extractor

    def _initialize(self):
        """
        Reads file headers and initializes data structures
        """
        for var in self.varlist:
            v, e = self._get_extractor(var)
            self.extractors[var] = e

    def read_stack(self, stack, comp_z_coords=True):
        """Reads 3D fields from one file to memory"""
        # read fields
        data = {}
        for var in self.extractors.keys():
            filename = self.extractors[var].generateFileName(stack)
            if not os.path.isfile(filename):
                raise IOError('File not found: ' + filename)
            ncfile = self.extractors[var].getNCFile(stack)

            v, filetypestr = splitVarToFileType(var)
            varstr = getNCVariableName(v)
            if varstr in vectorVars:
                varlist = [varstr + '_x', varstr + '_y']  # hvel_x, hvel_y
            else:
                varlist = [varstr]

            for v in varlist:
                d = ncfile.variables[v][:].astype(np.float32)
                data[v] = np.ma.masked_invalid(d, copy=False)

        extractor = self.extractors.values()[0]

        time = ncfile.getTime()
        facenodes = self.meshsearchobj.faceNodes
        x = self.meshsearchobj.node_x
        y = self.meshsearchobj.node_y
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

        # mask out bad data values
        print 'dry nodes', np.sum(drynodes)
        for var in data.keys():
            d = data[var]
            d[d < -900] = np.nan
            d.mask[:, :, :] = drynodes[:, np.newaxis, :]
            d[d.mask] = np.nan
        # bundle everything up for output
        data['x'] = x
        data['y'] = y
        if comp_z_coords:
            data['z'] = z
        data['time'] = time
        data['dry_nodes'] = drynodes
        data['elev'] = elev
        data['depth'] = depth

        return data


def integrate_scalar_field_at_nodes(field_array, data, meshsearchobj):
    """
    Computes volume integral of a field defined at 3D nodes (triangle vertices and whole levels)

    field_array : ndarray (nTime, nVert, nNodes)
        scalar field to integrate
    data : dict
        dictionary from read_3d_mesh_data function with keys
        'x', 'y', 'z'.

    Returns
    -------
    total_integral : array (nTime, )
        volume integral for each time step
    """
    warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')

    x = data['x']
    y = data['y']
    z = data['z']
    facenodes = meshsearchobj.faceNodes

    # compute element height
    dz = np.diff(z, axis=1)
    depth = np.nanmax(z, axis=1) - np.nanmin(z, axis=1)

    # compute nodal volume
    # integrate over vertical edges
    elem_int = 0.5*dz*(field_array[:, 1:, :] + field_array[:, :-1, :])
    # multiply by horizontal area
    areas = computeAreas(facenodes, x, y)  # area of each triangle
    node_area = np.zeros_like(x)
    scaled_area = areas/3
    node_area[facenodes[:, 0]] += scaled_area
    node_area[facenodes[:, 1]] += scaled_area
    node_area[facenodes[:, 2]] += scaled_area

    elem_int = elem_int*node_area
    np.ma.fix_invalid(elem_int, copy=False, fill_value=0.0)
    total_integral = np.sum(elem_int, axis=1).sum(axis=1)
    return total_integral


class ScalarIntegrator(object):
    """
    Computes volume integrals of SELFE 3D output fields
    """
    def __init__(self, datadir, integrant_var, criterion_var=None,
                 integrant_lim=None, criterion_lim=None):
        self.datadir = datadir
        self.integrant_var = integrant_var
        self.integrant_lim = integrant_lim if integrant_lim is not None else [None, None]
        self.criterion_var = criterion_var
        self.criterion_lim = criterion_lim if criterion_lim is not None else [None, None]

        varlist = [integrant_var]
        if self.criterion_var:
            varlist.append(self.criterion_var)
        self.reader = SelfeOutputReader(self.datadir, varlist)

        self.output_streams = []

    def add_output_stream(self, stream):
        """Adds a stream where the result will be dumped"""
        self.output_streams.append(stream)

    def apply_limits(self, array, criterion_array):
        """
        Applies limits to array in place

        All entries out of range will be set to zero
        """
        warnings.filterwarnings('ignore', 'invalid value encountered in less')
        warnings.filterwarnings('ignore', 'invalid value encountered in greater')

        # limit integral to nodes that are within the range
        mask = np.zeros_like(array, dtype=bool)

        if self.integrant_lim[0] is not None:
            ix = array < self.integrant_lim[0]
            mask[ix] = 1
        if self.integrant_lim[1] is not None:
            ix = array > self.integrant_lim[1]
            mask[ix] = 1

        if criterion_array is not None:
            if self.criterion_lim[0] is not None:
                ix = criterion_array < self.criterion_lim[0]
                mask[ix] = 1
            if self.criterion_lim[1] is not None:
                ix = criterion_array > self.criterion_lim[1]
                mask[ix] = 1

        array[mask] = 0.0

    def process_stack(self, stack):
        """
        Computes integral of scalar field for each time step in the output file

        Returns
        -------

        time : array (nTime, )
            time stamps of the output data
        integral : array (nTime, )
            volume integral of the scalar field
        """
        data = self.reader.read_stack(stack)
        k = getNCVariableName(splitVarToFileType(self.integrant_var)[0])
        integrant_array = data[k]
        if self.criterion_var is not None:
            k = getNCVariableName(splitVarToFileType(self.criterion_var)[0])
            criterion_array = data[k]
        else:
            criterion_array = None
        self.apply_limits(integrant_array, criterion_array)

        integral = integrate_scalar_field_at_nodes(integrant_array, data,
                                                   self.reader.meshsearchobj)
        return data['time'], integral

    def process_time_range(self, starttime=None, endtime=None, stacks=None):
        """
        Processes given time range (or stacks if provided)
        """
        if stacks is None:
            stacks = self.reader.extractors.values()[0].dataFile.getStacks(starttime, endtime, wholeDays=True)

        for stack in stacks:
            t, v = self.process_stack(stack)
            for ostream in self.output_streams:
                ostream.append(t, v)
        for ostream in self.output_streams:
            ostream.close()


def makeVolumeDataContainer(times, vals, location, runtag,
                            variable):
    """
    Creates a dataContainer for time series data
    """
    ta = timeArray.timeArray(times, 'epoch')
    x = np.array([0])
    y = np.array([0])
    z = np.array([0])
    data = vals[np.newaxis, np.newaxis, :]
    meta = {}
    meta['location'] = location
    meta['instrument'] = 'model'
    meta['variable'] = variable
    meta['dataType'] = 'volintegral'
    meta['tag'] = runtag
    fieldNames = [variable]

    dc = dataContainer.dataContainer('', ta, x, y, z, data,
                                     fieldNames, coordSys='spcs', metaData=meta)
    return dc


class TimeSeriesDCStream(object):
    """
    Concatenates time series data and dumps it to disk in netcdf format
    """
    def __init__(self, variable, runtag, location, rule='singleFile'):
        self.dc = None
        self.variable =  variable
        self.runtag = runtag
        self.location = location
        self.dc_factory = lambda t, v : makeVolumeDataContainer(t, v, self.location, self.runtag, self.variable)

    def append(self, time, values):
        """Appends data to the data set"""
        new_dc = self.dc_factory(time, values)
        if self.dc is None:
            self.dc = new_dc
        else:
            self.dc.mergeTemporal(new_dc)

    def close(self):
        """
        Closes the stream

        Saves concatenated dataContainer to disk
        """
        print self.dc
        rule = 'singleFile'
        dtm.saveDataContainerInTree(self.dc, rule=rule, dtype=np.float32,
                                    overwrite=True)


# -----------------------
# command line interface
# -----------------------


def parse_commandline():
    import argparse

    parser = argparse.ArgumentParser(description='Computes volume integrals of a scalar field from SELFE outputs')
    parser.add_argument('-r', '--runtag', required=True,
                        help='Run tag, used as a label in post-proc.')
    parser.add_argument('-d', '--datadir',  required=True,
                        help='directory where model outputs are stored')
    parser.add_argument('-s', '--startstr',
                        help='Date to start processing')
    parser.add_argument('-e', '--endstr',
                        help='Date to end processing. If omitted, start time plustidal day is used.')
    parser.add_argument('--stacks',
                        dest='stacklim', help='range of output files to read (e.g 1,14) if start,end not given')
    parser.add_argument('--loc', required=True,
                        help='Human readable name of the integral domain, e.g. "etm" or "lowestuary"')
    parser.add_argument('-v', '--variable', required=True,
                        help='variable to integrate: e.g. "salt". To use a specific output file define extension, e.g. "salt.70"')
    parser.add_argument('-l', '--low-threshold', type=float,
                        help='Lower bound of integrant. Values below this value will be omitted from the integral')
    parser.add_argument('-u', '--up-threshold', type=float,
                        help='Upper bound of integrant. Values above this value will be omitted from the integral')
    parser.add_argument('-c', '--criterion-variable',
                        help='Additional criterion scalar variable to define integral domain. Used in combination with --low-cri-threshold and --up-cri-threshold')
    parser.add_argument('-L', '--low-cri-threshold', type=float,
                        help='Lower bound of criterion scalar. Areas where criterion variable is below this value will be omitted from the integral')
    parser.add_argument('-U', '--up-cri-threshold', type=float,
                        help='Upper bound of criterion scalar. Areas where criterion variable is above this value will be omitted from the integral')
    parser.add_argument('--save-in-tree', action='store_true', dest='save_monthly_file',
                        help='saves extracted data in file tree with monthly files instead of a single file', default=False)
    parser.add_argument('-T', '--tracermodel',
                        help='Enable extraction of tracers: sed, oxy, generic. Must '
                        'supply number of tracers for \'sed\' and \'generic\' '
                        'models via the -N switch. \'oxy\' model provides tracers: '
                        '\'NO3\',\'NH4\',\'phy\',\'zoo\',\'det\' and \'oxy\'.', default=None)
    parser.add_argument('-N', '--numtracers', type=int,
                        help='Number of tracers to support for \'sed\' and \'generic\' models',
                        default=None)

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

    print('runtag        : {:}'.format(args.runtag))
    print('data dir      : {:}'.format(args.datadir))
    if st is not None:
        print('start time    : {:}'.format(st))
        print('end time      : {:}'.format(et))
    else:
        print('stacks        : {:}'.format(stacks))
    print('variable      : {:}'.format(args.variable))
    if args.low_threshold is not None:
        print('  min value   : {:}'.format(args.low_threshold))
    if args.up_threshold is not None:
        print('  max value   : {:}'.format(args.up_threshold))
    print('variable      : {:}'.format(args.criterion_variable))
    if args.low_cri_threshold is not None:
        print('  min value   : {:}'.format(args.low_cri_threshold))
    if args.up_cri_threshold is not None:
        print('  max value   : {:}'.format(args.up_cri_threshold))
    print('save-in-tree  : {:}'.format(args.save_monthly_file))

    if args.tracermodel:
        if not args.numtracers and args.tracermodel.split('.')[0] in ['sed', 'generic']:
            raise Exception('numTracers must be provided if sed or generic tracer models are used.')
        addTracers(args.tracermodel, numTracers=args.numtracers)

    si = ScalarIntegrator(args.datadir, args.variable, args.criterion_variable,
                          integrant_lim=[args.low_threshold, args.up_threshold],
                          criterion_lim=[args.low_cri_threshold, args.up_cri_threshold])

    v, ext = splitVarToFileType(args.variable)
    rule = 'monthlyFile' if args.save_monthly_file else 'singleFile'
    ostream = TimeSeriesDCStream(v, args.runtag, args.loc, rule=rule)
    si.add_output_stream(ostream)

    si.process_time_range(st, et, stacks=stacks)


if __name__ == '__main__':
    parse_commandline()
