"""
Read SELFE time history (.th) files to a data container.

Jesse Lopez  - 2016-04-15
"""
import datetime
import argparse

import numpy as np

from crane.data import timeArray
from crane.data import dataContainer


class thParser(object):

    def __init__(self, filename, start_time):
        self.filename = filename
        self.start_date = start_time
        self.time = None
        self.data = None

    def readFile(self):
        """Read time history file."""
        th = np.loadtxt(self.filename)
        self.time = timeArray.simulationToEpochTime(th[:, 0], self.start_date)
        self.data = th[:, 1]

    def genDataContainer(self, variable='variable', station='bvao',
                         depth='0', bracket='A', save=False):
        """Generate data container."""
        x = y = z = 0
        coordSys = ''

        meta = {}
        meta['tag'] = 'timeHistory'
        meta['variable'] = variable
        meta['location'] = station
        meta['msldepth'] = depth
        meta['bracket'] = bracket
        dc = dataContainer.dataContainer.fromTimeSeries(
            self.time, self.data, fieldNames=[variable],
            x=x, y=y, z=z, timeFormat='epoch', coordSys=coordSys,
            metaData=meta)

        if save:
            fname = './'+station+'_'+variable+'_'+'0'+'_'+self.start_date.strftime('%Y-%m-%d')+'.nc'
            print fname
            dc.saveAsNetCDF(fname)
        return dc


def parseCommandLine():
    parser = argparse.ArgumentParser(description='Read time history to dataContainer.')
    parser.add_argument('filepath', type=str, help='Path to time history file.')
    parser.add_argument('starttime', type=str, help='Start time of simulation YYYY-MM-DD')
    parser.add_argument('variable', type=str, help='Variable name (e.g. - salinity, temp, turbidity)')
    parser.add_argument('station', type=str, help='Station name (e.g. - saturn01, tpoin)')
    parser.add_argument('depth', type=str, help='Station depth (e.g. - 0.1, 4.0)')
    parser.add_argument('bracket', type=str, help='Bracket (e.g. - F, A, R)')

    args = parser.parse_args()
    st = datetime.datetime.strptime(args.starttime, '%Y-%m-%d')

    th = thParser(args.filepath, st)
    th.readFile()
    th.genDataContainer(args.variable, args.station, args.depth, args.bracket, True)

if __name__ == '__main__':
    parseCommandLine()
