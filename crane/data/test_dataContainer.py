
# test dataContainer and timeArray classes

import timeArray as t
import dataContainer as d
from scipy import *
from datetime import datetime

x = array([ 0.1, 0.2, 0.3, 0.4 ] )
y = array([ 1.0, 1.0, 1.0, 1.0 ] )
z = array([ 0.0, 0.0, 0.0, 0.0 ] )

data0 = array([[0, 1, 2],
               [3, 4, 5],
               [6, 7, 8],
               [3, 3, 3]])
data1 = dstack( (data0, 2*data0, 3*data0 ) )
data2 = dstack( (4*data0, 5*data0, 6*data0 ) )

ta = t.timeArray(2193.0+arange(3)/4.0,'corie')
tb = t.timeArray(2193.0+arange(3,6)/4.0,'corie')
tc = t.timeArray(2193.5+arange(3)/8.0,'corie')

print '\ntimeArray:\n'
print 'ta:'
print 'ta:',ta.array
print ta, '\n'
print 'tb:'
print 'tb:',tb.array
print tb, '\n'
print 'tc:'
print 'tc:',tc.array
print tc, '\n'

# test conversions
print '\nConversion:'
sd = datetime(1996,12,06,0,0,1)
print 'corie',ta[0], ta.asEpoch().asSimulation(sd).asCorie()[0]
print 'epoch',ta.asEpoch()[0], ta.asSimulation(sd).asEpoch()[0]

dc1 = d.dataContainer('container1', ta, x,y,z, data1, ['a','b','c'])
dc2 = d.dataContainer('container2', tb, x,y,z, data2, ['a','b','c'])
dc3 = dc1.copy()
dc3.description = 'container3'
print '\ndataContainer:\n'
print dc1, '\n'
print dc2, '\n'
print dc3, '\n'
# test merge
dc3.mergeTemporal(dc2)
print '\nMerge time:'
print 'dc1:', dc1.time.array
print 'dc2:', dc2.time.array
print 'dc3:', dc3.time.array
print 'dc1:', dc1.data[0,1,:]
print 'dc2:', dc2.data[0,1,:]
print 'dc3:', dc3.data[0,1,:]

# test interpolation
print '\nInterpolate:'
dc4 = dc3.interpolateInTime(tc)
dc3.data[:,:]
print 'dc3:', dc3.time.array
print 'dc3:', dc3.data[0,1,:]
print 'dc4:', dc4.time.array
print 'dc4:', dc4.data[0,1,:]

print '\nMerge fields:'
dc5 = dc1.copy()
dc5.mergeFields(dc1)
print 'dc1:', dc1.fieldNames, dc1.data.shape
print 'dc5:', dc5.fieldNames, dc5.data.shape

# test copied properties, should not change the parent
print '\nCopy dependency:'
dc4.fieldNames[0] = 'I am the eggman'
print 'dc3:', dc3.fieldNames
print 'dc4:', dc4.fieldNames

dc4.description = 'I am the walrus'
print 'dc3:', dc3.description
print 'dc4:', dc4.description

#dc4.saveToDisk('dc4','')

#dc5 = d.dataContainer.loadFromDisk('dc4.npz')
#print dc5

dc5.changeTimeFormat('simulation',sd)
dc5.description = 'test 5'
print dc5
print dc5.time.startDate
#dc5.saveToDisk('dc5.npz','')

#dc5 = d.dataContainer.loadFromDisk('dc5.npz')
#print dc5
#print dc5.time.startDate

dc6 = d.dataContainer('container1', ta, x[None,0],y[None,0],z[None,0], reshape(data1[1,1,:],(1,1,-1)), ['elev'], coordSys='lonlat')
dc6.setMetaData( 'one','traffic' )
dc6.setMetaData( 'two','little wing' )
print ' *** dc6 ***',dc6
dc6.saveAsNetCDF( 'test.nc' )
dc7 = d.dataContainer.loadFromNetCDF( 'test.nc' )
print ' *** dc7 ***',dc7

