#!/usr/local/anaconda/bin/python
#
# process database outputs salinity and temperature for climatology. Extract surface, bottom, and compute stats
# output to a netcdf files based on the original
#
import netCDF4 as nc
import numpy as np
import os
import sys
import datetime
from optparse import OptionParser

# create a netcdf file based on an existing netcdf file
def makeNetcdfFile(proto, fname, var, vsurf, vbot, vmin, vmax, vavg, vstd, vstrat):

    def addVar(ncout, name, v, varin, **kwargs):
        outVar = ncout.createVariable(name, 'float32', (u'time', u'node'), zlib=True)
        outVar[:] = v
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                # print "%s=%s" % (key, value)
                setattr(outVar, key, value)
        else:
            for aname in varin.ncattrs():
                # print aname, getattr(varin, aname)
                if aname != '_FillValue':
                    setattr(outVar, aname, getattr(varin, aname))

    ncid = nc.Dataset(proto, "r")
    ncout = nc.Dataset(fname, "w")
    for dname, dim in ncid.dimensions.iteritems():
            # print dname, len(dim)
        ncout.createDimension(dname, len(dim))

# Copy variables
    for vname, varin in ncid.variables.iteritems():
        # print vname,varin.datatype, varin.dimensions
        if vname == var:
            if var == 'salt':
                addVar(ncout, '''surface_%s''' % (var,), vsurf, varin,
                       standard_name='sea_surface_salinity',
                       units="1e-3",
                       long_name="Sea surface salinity",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''bottom_%s''' % (var,), vbot, varin,
                       standard_name='sea_bottom_salinity',
                       units="1e-3",
                       long_name="Sea bottom salinity",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''min_%s''' % (var,), vmin, varin,
                       standard_name='minimum_sea_water_salinity',
                       units="1e-3",
                       long_name="Minimum sea water salinity",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''max_%s''' % (var,), vmax, varin,
                       standard_name='maximum_sea_water_salinity',
                       units="1e-3",
                       long_name="Maximum sea water salinity",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''avg_%s''' % (var,), vavg, varin,
                       standard_name='average_sea_water_salinity',
                       units="1e-3",
                       long_name="Average sea water salinity",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''std_%s''' % (var,), vstd, varin,
                       standard_name='standard_deviation_sea_water_salinity',
                       units="1e-3",
                       long_name="Standard deviation of average sea water salinity",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                if vstrat != None:
                    addVar(ncout, 'strat', vstrat, varin,
                           standard_name='normalized_bulk_salinity_stratification',
                           units="dimensionless",
                           long_name="Dimensionless normalized bulk salinity stratification",
                           missing_value=-9999.0,
                           mesh="Mesh",
                           location="node",
                           coordinates="node_lon node_lat"
                           )
            elif var == 'temp':
                addVar(ncout, '''surface_%s''' % (var,), vsurf, varin,
                       standard_name='sea_surface_temperature',
                       units="degree_Celsius",
                       long_name="Sea surface temperature",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''bottom_%s''' % (var,), vbot, varin,
                       standard_name='sea_bottom_temperature',
                       units="degree_Celsius",
                       long_name="Sea bottom temperature",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''min_%s''' % (var,), vmin, varin,
                       standard_name='minimum_sea_water_temperature',
                       units="degree_Celsius",
                       long_name="Minimum sea water temperature",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''max_%s''' % (var,), vmax, varin,
                       standard_name='maximum_sea_water_temperature',
                       units="degree_Celsius",
                       long_name="Maximum sea water temperature",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''avg_%s''' % (var,), vavg, varin,
                       standard_name='average_sea_water_temperature',
                       units="degree_Celsius",
                       long_name="Average sea water temperature",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )
                addVar(ncout, '''std_%s''' % (var,), vstd, varin,
                       standard_name='standard_deviation_sea_water_temperature',
                       units="degree_Celsius",
                       long_name="Standard deviation of average sea water temperature",
                       missing_value=-9999.0,
                       mesh="Mesh",
                       location="node",
                       coordinates="node_lon node_lat"
                       )

        else:
            outVar = ncout.createVariable(
                vname, varin.datatype, varin.dimensions)
            outVar[:] = varin[:]
            for aname in varin.ncattrs():
                # print aname, getattr(varin, aname)
                if aname != '_FillValue':
                    setattr(outVar, aname, getattr(varin, aname))

    ncid.close()
    ncout.close()

# for a variable, get stuff
def getStats(v, skip):
    nt, nv, nn = np.shape(v)
    # print nt, nv, nn
    vsurf = v[:, nv - 1, :]
    vmin = np.min(v, axis=1)
    vmax = np.max(v, axis=1)
    vavg = np.mean(v, axis=1)
    vstd = np.std(v, axis=1)
    vbot = np.zeros((nt, nn))
    for i in range(0, 96, skip):
        vtmp = v[i, :,:]
        # print 'vtmp', np.shape(vtmp)
        ind = np.ma.count_masked(vtmp, axis=0)
        # print 'ind', np.shape(ind)
        vbot[i, :] = [vtmp[ind[j], j] for j in range(nn)]
    #vbot= v[:, 18,:]
    return (vsurf, vbot, vmin, vmax, vavg, vstd)

# extract surf, bot, compute stats, and write out file
def extractStats(db, syr, sday, skip, outdir):
    basedir = '''/home/workspace/ccalmr/hindcasts/%s/''' % (db,)
    eyr = syr
    eday = sday
    for j in range(syr, eyr + 1):
        odir = '''%d''' % (j,)
        if not os.path.exists(odir):
            os.makedirs(odir)
        ndays = datetime.datetime.strptime(
            '''%4d-12-31''' % (j), '%Y-%m-%d').timetuple().tm_yday
        for k in range(sday, eday + 1):
            fnsalt = '''%s/%d/outputs/%d_salt.63.nc''' % (basedir, j, k,)
            fntemp = '''%s/%d/outputs/%d_temp.63.nc''' % (basedir, j, k,)
            ncsalt = nc.Dataset(fnsalt)
            nctemp = nc.Dataset(fntemp)
            s = ncsalt.variables['salt'][:]
            (ssurf, sbot, smin, smax, savg, sstd) = getStats(s, skip)
            sstrat = ssurf * np.nan
            # for i in range(0,96):
            stmp = savg
            stmp[stmp <= 0.0] = np.nan
            sstrat = np.divide((smax - smin), stmp)
            #    sstrat = np.divide((smax[i,savg>0.0] - smin[savg>0.0]), savg[savg>0.0])
            #print np.shape(sstrat), np.shape(ssurf)
            t = nctemp.variables['temp'][:]
            (tsurf, tbot, tmin, tmax, tavg, tstd) = getStats(t, skip)
            #    print 'getStrat',j, k
            sname = '''%s/%d_saltstats.61.nc''' % (odir, k,)
            tname = '''%s/%d_tempstats.61.nc''' % (odir, k,)
            print 'making', sname
            makeNetcdfFile(
                fnsalt, sname, 'salt', ssurf, sbot, smin, smax, savg, sstd, sstrat)
            print 'making', tname
            makeNetcdfFile(
                fntemp, tname, 'temp', tsurf, tbot, tmin, tmax, tavg, tstd, None)
            print 'Done'
            ncsalt.close()
            nctemp.close()

if __name__ == '__main__':

    usage = (
        'Usage: %prog -y [year YYYY] -j [day of year] -o [path] -d [db id as in /home/workspace/ccalmr/hindcasts/] -k [skip]\n')

    parser = OptionParser(usage=usage)
    parser.add_option('-y', '--year', action='store', type='int',
                      dest='year', help='Year to process')
    parser.add_option('-j', '--jday', action='store', type='int',
                      dest='dayOfYear', help='Day of year to process (1-365/366)')
    parser.add_option('-o', '--outputDirectory', action='store', type='string',
                      dest='outDir', help='directory where generated variable stats are stored')
    parser.add_option('-d', '--database', action='store', type='string',
                      dest='db', help='Database ID (db31, db22)')
    parser.add_option('-k', '--skip', action='store', type='int',
                      dest='skip', help='Generate every skip -th time step. (Default: %default)', default=1)

    (options, args) = parser.parse_args()
    yr = options.year
    doy = options.dayOfYear
    outDir = options.outDir
    skip = options.skip
    db = options.db

    if len(sys.argv) < 2:
        parser.print_help()

    extractStats(db, yr, doy, skip, outDir)
