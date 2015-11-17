"""
Implementation of comma separated station file reader/writer.
Lines beginning with # are ignored.

Tuomas Karna 2013-11-05
"""
import csv
import os

from crane.data import collection


class csvStationFile(collection.tupleList):
    """Comma separated station file reader/writer

    Each station has attributes : 'x' and 'y'

    Each line in the station csv file is
    station_name, x_coord, y_coord

    Examples:
    saturn02,330029.000,284648.500
    saturn03,344208.000,287180.000
    """

    def __init__(self, source=None):
        keywords = ['location', 'x', 'y']
        super(csvStationFile, self). __init__(keywords, source)

    def readFromFile(self, filename, verbose=False):
        """Reads data from csv file
        Lines beginning with # are ignored
        """
        if filename is None:
            raise Exception('csvStationFile filename must be given.')
        if not os.path.isfile(filename):
            raise IOError('file not found: ' + filename)
        with open(filename, 'r') as f:
            reader = csv.reader(
                filter(lambda row: row[0] != '#', f),
                delimiter=',', skipinitialspace=True)
            for line in reader:
                try:
                    if len(line) == 0:
                        continue
                    if verbose:
                        print line
                    if len(line) < len(self.keywords):
                        raise Exception('Too few arguments')
                    loc = line[0]
                    x = float(line[1])
                    y = float(line[2])
                    self.addSample((loc, x, y))
                except Exception as e:
                    print e
                    raise Exception(
                        'Malformed line ' +
                        str(line) +
                        ' in file ' +
                        filename)
        return self

    def writeToDisk(self, filename):
        """Writes data to csv file"""
        fid = open(filename, 'w')
        for tup in self.getTuples():
            strList = []
            for i, v in enumerate(tup):
                if v is not None:
                    s = '%.5f' % v if i in [1, 2, 3] else v
                    strList.append(s)
            line = ', '.join(strList)
            if len(strList) < len(self.keywords):
                raise Exception('Malformed line: ' + line)
            fid.write(line + '\n')
        fid.close()

    def getStations(self):
        """Return a list of stations, may contain duplicates"""
        ix = self.keywords.index('location')
        return [t[ix] for t in self.getTuples()]

    def getX(self, location):
        """Returns the x coordinate of the given location"""
        tuples = self.getTuples(location=location)
        if len(tuples) == 0:
            raise Exception('Given station coordinates not found: ' + location)
        return tuples[0][self.keywords.index('x')]

    def getY(self, location):
        """Returns the y coordinate of the given location"""
        tuples = self.getTuples(location=location)
        if len(tuples) == 0:
            raise Exception('Given station coordinates not found: ' + location)
        return tuples[0][self.keywords.index('y')]

    def getLocation(self, location):
        """Returns x,y coordinates of a station. For backwards compatibility with
        StationFile class."""
        return (self.getX(location), self.getY(location))

    def sort(self):
        """Sorts the tuples in place."""
        self.data = sorted(self.data,
                           key=lambda o: o[self.keywords.index('location')])


class csvStationFileWithDepth(csvStationFile):
    """Comma separated station file reader/writer with depth and variable
    information.

    Each station has attributes: 'x', 'y', 'z', 'zType', 'variable'
    This object is useful when defining time series to extract.

    Each line in the station csv file is
    station_name, x_coord, y_coord, z_coord, zType, variable

    Examples:
    saturn02,330029.000,284648.500,6.000,depth,salt
    saturn03,344208.000,287180.000,-2.400,z,salt
    """

    def __init__(self, source=None):
        self.reg_keywords = ['location', 'x', 'y', 'z', 'zType']
        self.opt_keywords = ['variable']
        allkeywords = self.reg_keywords + self.opt_keywords
        super(csvStationFile, self).__init__(allkeywords, source)

    def readFromFile(self, filename):
        """Reads data from csv file"""
        with open(filename, 'r') as f:
            reader = csv.reader(
                filter(lambda row: row[0] != '#', f),
                delimiter=',', skipinitialspace=True)
            for line in reader:
                try:
                    if len(line) == 0:
                        continue
                    if len(line) < len(self.reg_keywords):
                        raise Exception('Too few arguments')
                    loc = line[0]
                    x = float(line[1])
                    y = float(line[2])
                    z = float(line[3])
                    zType = line[4]
                    if len(line) > len(self.reg_keywords):
                        var = line[5]
                    else:
                        var = None
                    self.addSample((loc, x, y, z, zType, var))
                except Exception as e:
                    print e
                    raise Exception(
                        'Malformed line ' +
                        str(line) +
                        ' in file ' +
                        filename)

    def writeToDisk(self, filename):
        """Writes data to csv file"""
        fid = open(filename, 'w')
        for tup in self.getTuples():
            strList = []
            for i, v in enumerate(tup):
                if v is not None:
                    s = '%.5f' % v if i in [1, 2, 3] else v
                    strList.append(s)
            line = ', '.join(strList)
            if len(strList) < 3:
                raise Exception('Malformed line: ' + line)
            fid.write(line + '\n')
        fid.close()

    def getZ(self, location):
        """Returns the z coordinate of the given location"""
        tuples = self.getTuples(location=location)
        if len(tuples) == 0:
            raise Exception('Given station coordinates not found: ' + location)
        return tuples[0][self.keywords.index('z')]

    def sort(self):
        """Sorts the tuples in place."""
        self.data = sorted(
            self.data, key=lambda o: o[
                self.keywords.index('z')])
        self.data = sorted(self.data,
                           key=lambda o: o[self.keywords.index('location')])
        self.data = sorted(self.data,
                           key=lambda o: o[self.keywords.index('variable')])
