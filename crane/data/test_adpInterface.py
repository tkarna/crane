import unittest
import time
import datetime
from crane.data import adpInterface
from crane.plotting import profilePlot
from crane.utility import saveFigure


class TestNetCDFInterface(unittest.TestCase):

    def test_oneMonth(self):
        Vs = 'vel_mag'
        sT = datetime.datetime(2012, 11, 1)
        eT = datetime.datetime(2012, 11, 30)
        dc = adpInterface.getADPData('saturn01.1950.A.ADP', sT, eT, Vs)
        print dc

    def test_twoMonths(self):
        Vs = 'vel_mag'
        sT = datetime.datetime(2011, 11, 1)
        eT = datetime.datetime(2011, 12, 30)
        dc = adpInterface.getADPData('saturn01.1950.A.ADP', sT, eT, Vs)
        print dc

    def test_2013(self):
        Vs = 'vel_mag'
        sT = datetime.datetime(2013, 1, 1)
        eT = datetime.datetime(2013, 1, 30)

        self.failUnlessRaises(
            Exception,
            getADPData,
            'saturn01.1950.A.ADP',
            sT,
            eT,
            Vs)

    def test_Components(self):  # FAILS
        Vs = 'vel_N'
        sT = datetime.datetime(2012, 11, 1)
        eT = datetime.datetime(2012, 11, 30)
        dc = adpInterface.getADPData('saturn01.1950.A.ADP', sT, eT, Vs)
        print dc
        Vs = 'vel_E'
        dc = adpInterface.getADPData('saturn01.1950.A.ADP', sT, eT, Vs)
        print dc
        Vs = 'alongvel'
        dc = adpInterface.getADPData('saturn01.1950.A.ADP', sT, eT, Vs)
        print dc
        Vs = 'crossvel'
        dc = adpInterface.getADPData('saturn01.1950.A.ADP', sT, eT, Vs)
        print dc
        Vs = 'vel_vert'
        dc = adpInterface.getADPData('saturn01.1950.A.ADP', sT, eT, Vs)
        print dc

    def test_plot(self):

        Vs = 'vel_mag'
        sT = datetime.datetime(2011, 11, 15)
        eT = datetime.datetime(2011, 12, 4)
        dc = adpInterface.getADPData('saturn01.1950.A.ADP', sT, eT, Vs)

        print dc

        dia = stackProfileTimeSeriesDC(clabel='Velocity mag', unit='m/s')
        dia.addSample('one', dc)
        dia.showColorBar()
        saveFigure('', 'adcp', 'png', verbose=True)

        '''
        sT = time.mktime(datetime.datetime(2012, 11, 1).timetuple())
        eT = time.mktime(datetime.datetime(2012, 11, 30).timetuple())
        Vs = ['vel_mag', 'bindepth']
        (t, v, u) = adpInterface.getADPData('saturn01', 'saturn01.1950.A.ADP', sT, eT, Vs)

        self.assertTrue(u.has_key('bindepth'))
        self.assertTrue(v['bindepth'].size>0)
        self.assertEqual(v['bindepth'].size, t.size)
        '''

    '''
     def test_ncReader(self):

        sT = datetime.datetime(2012, 11, 1)
        eT = datetime.datetime(2012, 11, 30)
        off = {'location':'saturn01', 'msldepth':'1950', 'bracket':'A', 'instrument':'ADP'}
        off['variable'] = 'bindepth'
        nc = netcdfCacheReader( off )
        t,d = nc.getData( sT, eT )

        self.assertEqual(t.size, d.size)
     '''


if __name__ == "__main__":
    unittest.main()
