"""
Generic methods for converting data between different spatial coordinate systems.
Uses pyproj library.

Tuomas Karna 2013-01-15
"""

import pyproj
import numpy as np
import collections

UTM_ZONE10 = pyproj.Proj(proj='utm',zone=10,datum='NAD83',units='m',errcheck=True)
SPCS_N_OR = pyproj.Proj(init='nad27:3601',errcheck=True)
LL_WGS84 = pyproj.Proj(proj='latlong',datum='WGS84',errcheck=True)
LL_WO = pyproj.Proj(proj='latlong',nadgrids='WO',errcheck=True) # HARN

def convertCoords( x, y, fromSys, toSys ) :
  """Converts x,y from fromSys to toSys."""
  if isinstance(x,np.ndarray) :
    # proj may give wrong results if nans in the arrays
    lon = np.ones_like(x)*np.nan
    lat = np.ones_like(y)*np.nan
    goodIx = np.logical_and( np.isfinite(x), np.isfinite(y) )
    lon[goodIx],lat[goodIx] = pyproj.transform( fromSys, toSys, x[goodIx], y[goodIx] )
  else :
    lon,lat = pyproj.transform( fromSys, toSys, x, y )
  return lon,lat

def spcs2lonlat( x, y ) :
  """Converts SPCS to longitude-latitude."""
  return convertCoords( x, y, SPCS_N_OR, LL_WO )

def lonlat2spcs( lon, lat ) :
  """Converts longitude-latitude to SPCS."""
  return convertCoords( lon, lat, LL_WO, SPCS_N_OR )

def spcs2utm( x, y ) :
  """Converts SPCS to longitude-latitude."""
  return convertCoords( x, y, SPCS_N_OR, UTM_ZONE10 )

def utm2spcs( lon, lat ) :
  """Converts longitude-latitude to SPCS."""
  return convertCoords( lon, lat, UTM_ZONE10, SPCS_N_OR )

def WGS842spcs( lon, lat ) :
  """Converts longitude-latitude to SPCS."""
  return convertCoords( lon, lat, LL_WGS84, SPCS_N_OR )

def rotateVector( lon, lat, u, v, csys, inverse=False, xySystem=None ) :
  """Rotates a vector defined in lonlat (LL_WO) to csys coodinate system.
  
   Args:
    lon - (scalar/singleton array) lon/x coordinate position
    lat - (scalar/singleton array) lat/y coordinate position
    u - (scalar/singleton array) lon/x velocity component
    v - (scalar/singleton array) lat/y velocity component
    csys - (pyproj.Proj) coordinate system to convert to
    inverse (bool) - if True perform inverse conversion from csys to lonlat
    xySystem (pyproj.Proj) - coordinate system where lat,lon are defined
   Output:
    u2 - (scalar/singleton array) rotated x velocity component
    v2 - (scalar/singleton array) rotated y velocity component
  """
  if isinstance( u, np.ndarray ) and isinstance( lon,np.ndarray ) :
    # lat,lon varies
    if len(lon) != len(lat) or len(lat) != len(u) or len(u) != len(v) :
      print len(lon), len(lat), len(u), len(v)
      raise Exception( 'array sizes do not match')
    u2 = np.zeros(u.shape)
    v2 = np.zeros(v.shape)
    for i in range(len(lon)) :
      u2[i],v2[i] = rotateVector( lon[i],lat[i],u[i], v[i], csys, inverse, xySystem )
    return u2,v2
  if xySystem != None :
    lon2,lat2 = convertCoords( lon, lat, xySystem, LL_WO )
  else :
    lon2 = lon
    lat2 = lat
  R,_ = getVectorRotationMatrix( lon2,lat2, csys )
  # forward rotation
  R = np.matrix(R)
  if inverse :
    R = np.matrix(np.linalg.inv(R))
  u2 = u.ravel()
  v2 = v.ravel()
  V = np.matrix( np.vstack((u2,v2)) )
  U = np.array(R*V)
  if isinstance( u, np.ndarray ) :
    return U[0,:],U[1,:]
  else :
    return float(U[0,]),float(U[1,])

def getVectorRotationMatrix( lon,lat, csys ) :
  # estimate angle of rotation from LL_WO to csys system.
  x,y = pyproj.transform( LL_WO, csys, lon, lat )
  deltaDegr = 1e-7 #1./3600.*0.01
  x2,y2 = pyproj.transform( LL_WO, csys, lon, lat+deltaDegr )
  dxdl = (x2-x)/deltaDegr
  dydl = (y2-y)/deltaDegr
  theta = np.arctan2(-dxdl,dydl)
  x2,y2 = pyproj.transform( LL_WO, csys, lon+deltaDegr, lat )
  dxdl = (x2-x)/deltaDegr
  dydl = (y2-y)/deltaDegr
  theta = np.arctan2(dydl,dxdl)
  c = np.cos( theta )
  s = np.sin( theta )
  R = np.array([[c,-s],[s,c]])
  return R, theta
