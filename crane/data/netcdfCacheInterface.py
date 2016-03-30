"""
Methods for fetching data from CMOP NetCDF cache.
"""

#-------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------
import os
import time
import numpy as np
from netCDF4 import Dataset as NetCDFFile

import datetime

from crane.data import dataContainer
from crane.files import stationFile

#-------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------
DB_HOST = 'cdb02'
DB_USER = 'reader'
DB_NAME = 'cmop'

# CF standard variable names
CF = dict([['salt', 'salinity'],
           ['temp', 'temperature'],
           ['cond', 'electrical conductivity'],
           ['elev', 'elevation'],
           ['pres', 'pressure'],
           ['fluores', 'fluorescence']])

# map internal variable names to those used in the netcdf cache
varToCacheVarName = {'NO3': 'nitrate', 'oxy': 'oxygen'}
cacheVarNameToVar = dict(zip(varToCacheVarName.values(),
                             varToCacheVarName.keys()))


def getAvailableOfferings(
    startTime,
    endTime,
    variables=[
        'elev',
        'temp',
        'salt']):
    """Returns a list of available offerings for given time period.

    Returns:
      offerings -- List of observation offering dicts available
    """
    import db as ncDB
    db = ncDB.DB(user=DB_USER, host=DB_HOST, dbname=DB_NAME)
    cacheVars = [varToCacheVarName.get(v, v) for v in variables]
    varStr = '(\'' + '\',\''.join(cacheVars) + '\')'
    sql = ('''select distinct offering,variable from instrument.offeringdetails
            where variable in %s
            and variable_visibility='public'
            and instrument is not null
            and instrumenttype is not null
            and deployedon <=  '%s'
            and (retrievedon >= '%s' or retrievedon is null)
            and station not in ('kiviuq','mch01','sqx01','muk01','ccho3')'''
           % (varStr, str(endTime), str(startTime)))
    results = db.execQuery(sql)
    offerings = []
    if results:
        for row in results:
            loc, dep, bra, instr = row[0].split('.')
            var = row[1]
            offeringDict = {}
            offeringDict['location'] = loc
            offeringDict['msldepth'] = dep
            offeringDict['bracket'] = bra
            offeringDict['instrument'] = instr
            offeringDict['variable'] = cacheVarNameToVar.get(var, var)
            offerings.append(offeringDict)
    db.close()

    return offerings


def getAllOfferings(variables=['elev', 'temp', 'salt']):
    """Returns a list of all offerings for the given variables.

    Returns:
      offerings -- List of observation offering dict available
    """
    import db as ncDB
    db = ncDB.DB(user=DB_USER, host=DB_HOST, dbname=DB_NAME)
    cacheVars = [varToCacheVarName.get(v, v) for v in variables]
    varStr = '(\'' + '\',\''.join(cacheVars) + '\')'
    sql = ('''select distinct offering,variable from instrument.offeringdetails
            where variable in %s
            and variable_visibility='public'
            and instrument is not null
            and instrumenttype is not null
            and station not in ('kiviuq','mch01','sqx01','muk01','ccho3')'''
           % (varStr))
    results = db.execQuery(sql)
    offerings = []
    if results:
        for row in results:
            loc, dep, bra, instr = row[0].split('.')
            var = row[1]
            offeringDict = {}
            offeringDict['location'] = loc
            offeringDict['msldepth'] = dep
            offeringDict['bracket'] = bra
            offeringDict['instrument'] = instr
            offeringDict['variable'] = cacheVarNameToVar.get(var, var)
            offerings.append(offeringDict)
    db.close()

    return offerings


def getAUVData(missionID, variable):
    """Fetches AUV data from the database for given missionID.
    NOTE: missionID is not Craig's numbering."""
    import db as ncDB
    db = ncDB.DB(user=DB_USER, host=DB_HOST, dbname=DB_NAME)
    # short variable name to db name
    variable = CF.get(variable, variable)
    sql = """select round(x(transform(location,32026))*0.3048,1),
            round(y(transform(location,32026))*0.3048,1),
            round(depth , 2),
            round(extract('epoch' from time),2),
            round(%s, 2),
            extract('epoch' from (time at time zone 'PST'))
            from auv.auvdata_cache where
            deploymentid=%d and %s notnull order by time""" % (variable, missionID, variable)
    results = db.execQuery(sql)

    x = []
    y = []
    depth = []
    var = []
    time = []

    if results:
        for row in results:
            x.append(float(row[0]))
            y.append(float(row[1]))
            depth.append(float(row[2]))
            time.append(float(row[3]))
            var.append(float(row[4]))
    else:
        print 'No results found for AUV missionID', missionID

    x = np.array(x)
    y = np.array(y)
    depth = np.array(depth)
    var = np.array(var)
    time = np.array(time)

    return x, y, depth, time, var

#"""
# def getADPData(station, offering, starttime, endtime, variables, q='PD0') :
#  """Fetches ADP data from nc cache"""
#  t, v, u = getncdatastation(station, offering, starttime, endtime, variables, quality=q)
#  return t, v, u
#"""


class netcdfCacheReader(object):

    def __init__(self, offeringDict):

        s = offeringDict['location']
        m = offeringDict['msldepth']
        b = offeringDict['bracket']
        i = offeringDict['instrument']
        v = offeringDict['variable']
        v = varToCacheVarName.get(v, v)
        self.nc_offering = '.'.join([s, m, b, i])
        self.offerString = '.'.join([s, m, b, i, v])
        self.station = s
        self.msldepth = m
        self.z = -float(m) / 100.0
        self.bracket = b
        self.var = v

        self.validmin = -1e15
        self.validmax = 1e15
        self.sql0 = '''select deploymentid, offering, msldepth, bracket,
                          instrumenttype,schemaname, tablename, columnname,
                          units, standardname, variable_validvaluemin,
                          variable_validvaluemax, variable_displaylimitminmax,
                          platformdescription, latitude, longitude
                   from instrument.offeringdetails
                   where offering = '%s' and variable = '%s' '''

        self.sql1 = '''select extract('epoch' from time) as seconds, %s
                   from %s.%s
                   where time between '%s' and '%s'
                         and %s notNULL
                         and msldepth = %s
                         and bracket = '%s'
                         and station = '%s'
                   order by time asc'''
        self.sql2 = '''select extract('epoch' from '%s'::timestamp) as seconds,
                          extract('epoch' from '%s'::timestamp) as seconds'''
        self.sql3 = '''select
                        case when date_trunc('day', '%s'::timestamp) = '%s'::timestamp
                                then to_char('%s'::timestamp ,'YYYY-MM-DD')
                                else to_char('%s'::timestamp,'YYYY-MM-DD HH24:MI:SS')
                        end,
                        case when date_trunc('day', '%s'::timestamp) = '%s'::timestamp
                                then to_char('%s'::timestamp,'YYYY-MM-DD')
                                else to_char('%s'::timestamp,'YYYY-MM-DD HH24:MI:SS')
                        end'''
        self.sql5 = '''select extract('timezone' from now())'''
        self.db = None
        self.detailsFetched = False

    def __del__(self):
        """Clean exit: disconnect"""
        self.disconnectFromDB()

    def connectToDB(self):
        """Connects to CMOP database to enable gathering of data."""
        try:
            import db as ncDB
            self.db = ncDB.DB(user=DB_USER, host=DB_HOST, dbname=DB_NAME)
        except ImportError:
            raise Exception('Unable to import Postgres')
        except:
            print DB_NAME, DB_HOST, DB_USER
            raise Exception('Unable to connect to the database')

    def disconnectFromDB(self):
        """Disconnects from CMOP database"""
        if self.db:
            self.db.close()

    def getOfferingDetails(self):
        """Gets details about offering from database."""
        if not self.db:
            self.connectToDB()

        sql = self.sql0 % (self.nc_offering, self.var)
        self.altvar = CF.get(self.var, self.var).lower().replace('_', ' ')
        results = self.db.execQuery(sql)
        if not results:
            raise Exception('Could not retrieve metadata ' + self.offerString)
        self.depid = results[0][0]
        self.station = results[0][1].split('.')[0]
        self.msldepth = results[0][2]
        self.bracket = results[0][3]
        self.instrumenttype = results[0][4]
        self.schemaname = results[0][5]
        self.tablename = results[0][6]
        self.columnname = results[0][7]
        self.units = '$' + results[0][8] + '$'
        self.units.replace('micro', r'\mu')
        self.standardname = results[0][9]
        if results[0][10] is not None:
            self.validmin = results[0][10]
        if results[0][11] is not None:
            self.validmax = results[0][11]
        self.displayminmax = results[0][12]
        self.platformdescription = results[0][13]
        self.lat = results[0][14]
        self.lon = results[0][15]

        self.detailsFetched = True

    def getData(self, startTime, endTime, quality='best'):
        """Gets data from NetCDF cache for this offering.

        The data collected from the NetCDF cache undergoes simple check to ensure
        it is valid.  Values that are determined to be too big, too small, or NaNs
        are removed from the data (both time and y arrays).

        Args:
          startTime -- Datetime of the time to start getting data
          endTime -- Datetime of the time to stop getting data
          source -- String indicating grabbing of data from NetCDF cache.
            Also a possibility to grab it from database, but not implemented here.
        Returns:
          time -- Numpy array of the time for the collected data
          y -- Numpy array for the values of the collected data
        """
        if not self.db:
            self.connectToDB()
        if not self.detailsFetched:
            self.getOfferingDetails()
        sql = self.sql2 % (startTime, endTime)
        results = self.db.execQuery(sql)
        startEpoch, endEpoch = results[0]
        sql = self.sql3 % (startTime, startTime, startTime, startTime,
                           endTime, endTime, endTime, endTime)
        results = self.db.execQuery(sql)
        self.startstr, self.endstr = results[0]
        sql = self.sql5
        results = self.db.execQuery(sql)
        self.tzoffset = results[0]

        print 'Getting nc file: ', self.nc_offering, startTime, endTime
        (time, y, units) = getncdatastation(self.station, self.nc_offering,
                                            startEpoch, endEpoch, vars=[self.var], quality=quality)
        if len(time) == 0:
            raise Exception('Retrieving data failed ' + self.offerString)
        y = y[self.altvar]
        valid = np.isfinite(y)
        valid = np.logical_and(valid, y >= float(self.validmin))
        valid = np.logical_and(valid, y <= float(self.validmax))

        time = time[valid]
        y = y[valid]
        if len(time) == 0:
            raise Exception(
                'Retrieving data failed after quality check ' +
                self.offerString)
        if self.var == 'elev' and self.station != 'hmndb' and self.z:
            print ' ************ elev corr **************** '
            correction = self.z
            if self.station == 'saturn06':
                correction = -self.z
            print y[:20]
            y = y + correction
            print y[:20]

        self.disconnectFromDB()

        return time, y, units.get(self.var, '')


def getDataContainerFromOffering(
        offeringDict,
        startTime,
        endTime,
        quality='best'):
    nc = netcdfCacheReader(offeringDict)
    t, d, unit = nc.getData(startTime, endTime, quality)

    badIx = np.nonzero(np.diff(t) <= 0)[0]
    if len(badIx) > 0:
        print 'Time sequence not monotonically increasing, attempting to fix ...'
        print badIx
        # print t[badIx:badIx+2], np.diff(t[badIx:badIx+2])
        goodIx = np.nonzero(np.diff(t) > 0)[0]
        print goodIx.shape
        goodIx = np.hstack(([0], goodIx + 1))
        t = t[goodIx]
        d = d[goodIx]

    sta = stationFile.StationFile()
    x, y = sta.getLocation(nc.station)
    z = nc.z
    meta = {}
    meta['dataType'] = 'timeseries'
    meta['location'] = nc.station
    meta['msldepth'] = str(int(round(nc.msldepth)))
    meta['bracket'] = nc.bracket
    meta['instrument'] = nc.instrumenttype
    meta['variable'] = cacheVarNameToVar.get(nc.var, nc.var)
    if unit:
        meta['unit'] = unit
    meta['tag'] = 'obs'
    if meta['location'] == 'saturn01' and meta['msldepth'] == '0':
        meta['dataType'] = 'profiler_raw'
        meta.pop('msldepth')
    # TODO fetch and add all other useful data too, samplingrate etc
    return dataContainer.dataContainer.fromTimeSeries(
        t, d, [meta['variable']], x=x, y=y, z=z, timeFormat='epoch',
        coordSys=sta.coordSys, metaData=meta)

# copied from cmop.ncdataextract


def openncFile(filename):
    print filename
    return NetCDFFile(filename, 'r')


class DataError(Exception):
    pass


def selecttimeindex(ncfile, st, en):
    tcolumn = ncfile.variables['time']
    a = tcolumn[0]
    b = tcolumn[tcolumn.shape[0] - 1]
    c = 0
    d = tcolumn.shape[0]
    # starttime in range of file?
    if st > a:
            # index to extract extraction
        c = np.searchsorted(tcolumn[:], st)
    if en < b:
        # index where to finish extraction
        d = np.searchsorted(tcolumn[:], en)
    return (c, d)

DEF_BASEDIR = '/home/workspace/ccalmr/data/nc/'


def getncdatastation(
        station,
        offering,
        starttime,
        endtime,
        vars=['all'],
        quality='best',
        basedir=DEF_BASEDIR):
    if quality == 'best':
        home_alt0 = os.path.join(basedir, 'PD0/stations/')
        home_alt1 = os.path.join(basedir, 'PD1/stations/')
        home = os.path.join(basedir, 'PD2/stations/')
        path_alt1 = '%s%s/%s' % (home_alt1, station, offering)
        path_alt0 = '%s%s/%s' % (home_alt0, station, offering)
    elif quality in ('level_0', 'PD0', 'raw'):
        home = os.path.join(basedir, 'PD0/stations/')
    elif quality in ('level_1', 'PD1', 'preliminary'):
        home = os.path.join(basedir, 'PD1/stations/')
    elif quality in ('level_2', 'PD2', 'verified'):
        home = os.path.join(basedir, 'PD2/stations/')
    else:
        raise NameError('quality level %s is unknown.')
    path = '%s%s/%s' % (home, station, offering)
    #cmop.debug('getting file from %s' % path)
    startt = time.localtime(starttime)
    endt = time.localtime(endtime)
    startyear = time.strftime('%Y', startt)
    startmonth = time.strftime('%m', startt)
    endyear = time.strftime('%Y', endt)
    endmonth = time.strftime('%m', endt)
    agg_vars = {}
    agg_units = {}
    agg_time = np.array([])
    cf = dict(
        [["salt", "salinity"],
         ["temp", "temperature"],
         ["cond", "electrical conductivity"],
         ["elev", "elevation"],
         ["pres", "pressure"],
         ["fluores", "fluorescence"]])
    vars = [cf.get(v, v).lower() for v in vars]
    vars = [v.replace("_", " ") for v in vars]
    for year in range(int(startyear), int(endyear) + 1):
        sm = 1
        em = 12
        if year == int(startyear):
            sm = int(startmonth)
        if year == int(endyear):
            em = int(endmonth)
        #cmop.debug('%s,%s,%s,%s,%s,%s,%s' % (startyear,year,endyear,startmonth,sm,em,endmonth))
        for month in range(sm, em + 1):
            try:
                file = '%s/%s%02d.nc' % (path, year, month)
                #cmop.debug('reading data from %s' % file)
                if not os.path.exists(file):
                    #cmop.debug('file %s not found.' % file)
                    if quality == 'best':
                        file = '%s/%s%02d.nc' % (path_alt1, year, month)
                        if not os.path.exists(file):
                            #cmop.debug('PD1 file %s not found.' % file)
                            file = '%s/%s%02d.nc' % (path_alt0, year, month)
                            if not os.path.exists(file):
                                #cmop.debug('PD0 file %s not found.' % file)
                                continue
                        #cmop.debug('using file, %s' % file)
                    else:
                        continue
                ncfile = openncFile(file)
                (c, d) = selecttimeindex(ncfile, starttime, endtime)
                tmp_time = ncfile.variables['time'][c:d]
                #cmop.debug('found %s values in %s' % (len(tmp_time),file))
                #time_ind = None
                # remove times outside specified time range
                # if year == int(startyear) and sm == month and year == int(endyear) and em == month:
                #   time_ind = ((tmp_time > starttime) | (tmp_time < endtime)).nonzero()
                #   tmp_time = tmp_time[time_ind]
                #   cmop.debug('reduced to %d values by time restriction' % (len(tmp_time)))
                # elif year == int(startyear) and sm == month:
                #   time_ind = (tmp_time > starttime).nonzero()
                #   tmp_time = tmp_time[time_ind]
                #   cmop.debug('first month: reduced to %d values by time restriction' % (len(tmp_time)))
                # elif year == int(endyear) and em == month:
                #   time_ind = (tmp_time < endtime).nonzero()
                #   tmp_time = tmp_time[time_ind]
                #   cmop.debug('reduced to %d values by time restriction' % (len(tmp_time)))
                # else: cmop.debug('(%s %s) "%s" != "%s" or "%s" != "%s" and %s
                # != %s or %s != %s' %
                # ((year==int(startyear)),(sm==month),year, startyear, sm,
                # month,year, endyear, em,month))
                variables = ncfile.variables.keys()
                # First file defines variables to aggregate
                if len(agg_vars.keys()) == 0:
                    varmap = []
                    for var in variables:
                        vname = var.replace(
                            '_', ' ').replace(
                            'water', '').strip()
                        if vname in vars or vars[0] == 'all':
                            agg_vars[vname] = np.array([])
                            varmap.append(var)
                if len(varmap) != len(vars):
                    #cmop.info('data file %s is missing expected variables: %s present, %s expected' % (file,', '.join(varmap),', '.join(vars)))
                    raise DataError(
                        'data file %s is missing expected variables: %s present, %s expected' %
                        (file, ', '.join(varmap), ', '.join(vars)))
                tmp_vars = agg_vars
                # Aggregating variables loop
                for var in varmap:
                    vname = var.replace('_', ' ').replace('water', '').strip()
                    agg_var = tmp_vars[vname]
                    values = ncfile.variables[var][c:d]
                    # if time_ind is not None:
                    #   values =  values[time_ind]
                    #cmop.debug('%s: %s ---> at %s, av %s, t %s, v %s' %  (file,var,len(agg_time),len(agg_var), len(values),len(tmp_time)))
                    values = np.array(values)
                    if len(values) != len(tmp_time):
                            #cmop.info('skipping data file %s, length of variable %s does not match length of time')
                        raise DataError(
                            'skipping data file %s, length of variable %s does not match length of time' %
                            (file, var))
                    bad = np.nonzero(values == 1099511627776.)[0]
                    values[bad] = np.nan
                    agg_var = np.append(agg_var, values)
                    tmp_vars[vname] = agg_var
                    tmp_var = ncfile.variables[var]
                    agg_units[vname] = getattr(tmp_var, 'units')
                    #values = np.append(np.array([None]), values)
                    #values = np.append(values, np.array([None]))
                agg_time = np.append(agg_time, tmp_time)
                agg_vars = tmp_vars
                ncfile.close()
            except DataError as e:
                raise e

    #cmop.debug('returning %s time values' % len(agg_time))
    # for (k,v) in agg_vars.items():
        #cmop.debug('%s returning %s values' % (k,len(v)))
    return (agg_time, agg_vars, agg_units)

if __name__ == '__main__':
    # simple test
    off = 'saturn02.1600.F.CT.temp'

    sT = datetime.datetime(2011, 0o5, 14)
    eT = datetime.datetime(2011, 0o7, 20)

    nc = netcdfCacheReader(off)
    t, d = nc.getData(sT, eT)
    print t
    print d

    dc = getDataContainerFromOffering(off, sT, eT)
    dc.saveAsNetCDF('test.nc')
    print dc
