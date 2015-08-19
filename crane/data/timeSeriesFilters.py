"""
Routine for filtering time series data.
NOTE could possibly be replaced by pandas or other library.

Tuomas Karna 2013-11-07
"""
import numpy as np
from collections import deque
from crane.data import dataContainer
from crane.data import timeArray

#------------------------------------------------------------------------------
# Low level routines with numpy arrays
#------------------------------------------------------------------------------

# period of M2 cycle in seconds
T_M2 = 44714.0

def computeRunningMean( t,x,T ) :
  """Computes running mean of the time series (t,x) with window length T.
  Returns new time series (tRes,xRes), one data point for each window.
  """
  tWin = deque([t[0]]) # one window of orig data
  xWin = deque([x[0]])
  tcum = t[0] # sum of windowed values
  xcum = x[0]
  tRes = [] # results
  xRes = []
  for i in range(1,len(t)) :
    # grow window to correct length
    tWin.append( t[i] )
    xWin.append( x[i] )
    tcum += t[i]
    xcum += x[i]
    tWinSpan = tWin[-1] - tWin[0]
    if tWinSpan >= T :
      # compute range/mean/whatever
      tRes.append( tcum/len(tWin) )
      xRes.append( xcum/len(xWin) )
      # shrink window by removing first element
      tcum -= tWin[0]
      xcum -= xWin[0]
      tWin.popleft()
      xWin.popleft()
  tRes = np.array(tRes)
  xRes = np.array(xRes)
  return tRes, xRes

def computeRunningMax( t,x,T ) :
  """Computes running max of the time series (t,x) with window length T.
  Returns new time series (tRes,xRes), one data point for each window.
  """
  tWin = deque([t[0]]) # one window of orig data
  xWin = deque([x[0]])
  tcum = t[0] # sum of windowed values
  xmax = x[0] # max of windowed values
  tRes = [] # results
  xRes = []
  for i in range(1,len(t)) :
    # grow window to correct length
    tWin.append( t[i] )
    xWin.append( x[i] )
    tcum += t[i]
    xmax = max(xmax,x[i])
    tWinSpan = tWin[-1] - tWin[0]
    if tWinSpan >= T :
      # compute range/mean/whatever
      tRes.append( tcum/len(tWin) )
      xRes.append( xmax )
      # shrink window by removing first element
      tcum -= tWin[0]
      if xmax == xWin[0] :
        # TODO remove list() for newer python
        xmax = max(list(xWin)[1:]) # need to recompute max
      tWin.popleft()
      xWin.popleft()
  tRes = np.array(tRes)
  xRes = np.array(xRes)
  return tRes, xRes

def computeRunningMin( t,x,T ) :
  """Computes running min of the time series (t,x) with window length T.
  Returns new time series (tRes,xRes), one data point for each window.
  """
  tWin = deque([t[0]]) # one window of orig data
  xWin = deque([x[0]])
  tcum = t[0] # sum of windowed values
  xmin = x[0] # min of windowed values
  tRes = [] # results
  xRes = []
  for i in range(1,len(t)) :
    # grow window to correct length
    tWin.append( t[i] )
    xWin.append( x[i] )
    tcum += t[i]
    xmin = min(xmin,x[i])
    tWinSpan = tWin[-1] - tWin[0]
    if tWinSpan >= T :
      # compute range/mean/whatever
      tRes.append( tcum/len(tWin) )
      xRes.append( xmin )
      # shrink window by removing first element
      tcum -= tWin[0]
      if xmin == xWin[0] :
        # TODO remove list() for newer python
        xmin = min(list(xWin)[1:]) # need to recompute min
      tWin.popleft()
      xWin.popleft()
  tRes = np.array(tRes)
  xRes = np.array(xRes)
  return tRes, xRes

def computeRunningRange( t,x,T ) :
  """Computes running range of the time series (t,x) with window length T.
  Returns new time series (tRes,xRes), one data point for each window.

  Range is defined as x_window.max()-x_window.min() for each window.
  """
  tWin = deque([t[0]]) # one window of orig data
  xWin = deque([x[0]])
  tcum = t[0] # sum of windowed values
  xmax = x[0] # max of windowed values
  xmin = x[0] # min of windowed values
  tRes = [] # results
  xRes = []
  for i in range(1,len(t)) :
    # grow window to correct length
    tWin.append( t[i] )
    xWin.append( x[i] )
    tcum += t[i]
    xmax = max(xmax,x[i])
    xmin = min(xmin,x[i])
    tWinSpan = tWin[-1] - tWin[0]
    if tWinSpan >= T :
      # compute range/mean/whatever
      tRes.append( tcum/len(tWin) )
      xRes.append( xmax - xmin )
      # shrink window by removing first element
      tcum -= tWin[0]
      if xmax == xWin[0] :
        # TODO remove list() for newer python
        xmax = max(list(xWin)[1:]) # need to recompute max
      if xmin == xWin[0] :
        # TODO remove list() for newer python
        xmin = min(list(xWin)[1:]) # need to recompute min
      tWin.popleft()
      xWin.popleft()
  tRes = np.array(tRes)
  xRes = np.array(xRes)
  return tRes, xRes

#------------------------------------------------------------------------------
# Routines with dataContainer as input/output
#------------------------------------------------------------------------------

def runningX( dc, T=T_M2, operator=computeRunningMean, gap_dt=None) :
  """Generic routine that applies given operator function for each window.
  T is window length."""
  nPoints,nFields,nTime = dc.data.shape
  # TODO extend to more points and fields
  if nPoints>1:raise NotImplementedError('Case with nPoints>1 not implemented')
  #if nFields>1:raise NotImplementedError('Case with nFields>1 not implemented')
  if nFields>1 : # process each component recursively
    components = []
    for i in range(nFields) :
      dc_comp = dc.extractField(i)
      dc_comp_new = runningX( dc_comp, T, operator )
      components.append( dc_comp_new )
    dc_new = components[0]
    for i in range(1,nFields) :
      dc_new.mergeFields( components[i] )
    return dc_new
  gaps,ranges,t = dc.detectGaps(dt=gap_dt)
  v = dc.data[0,0,:]
  x = dc.x[0,:] if dc.xDependsOnTime else dc.x
  y = dc.y[0,:] if dc.yDependsOnTime else dc.y
  z = dc.z[0,:] if dc.zDependsOnTime else dc.z
  tRes = np.array([])
  vRes = np.array([])
  if dc.xDependsOnTime : xRes = np.array([])
  if dc.yDependsOnTime : yRes = np.array([])
  if dc.zDependsOnTime : zRes = np.array([])
  for i in range( ranges.shape[0] ) :
    twin = t[ ranges[i,0]:ranges[i,1] ]
    vwin = v[ ranges[i,0]:ranges[i,1] ]
    if dc.xDependsOnTime : xwin = x[ ranges[i,0]:ranges[i,1] ]
    if dc.yDependsOnTime : ywin = y[ ranges[i,0]:ranges[i,1] ]
    if dc.zDependsOnTime : zwin = z[ ranges[i,0]:ranges[i,1] ]
    if len(twin) == 0 :
      continue
    t_op, v_op = operator( twin,vwin,T )
    tRes = np.append( tRes, t_op )
    vRes = np.append( vRes, v_op )
    if dc.xDependsOnTime : 
      tmean, xmean = computeRunningMean( twin,xwin,T )
      xRes = np.append( xRes, xmean )
    if dc.yDependsOnTime : 
      tmean, ymean = computeRunningMean( twin,ywin,T )
      yRes = np.append( yRes, ymean )
    if dc.zDependsOnTime : 
      tmean, zmean = computeRunningMean( twin,zwin,T )
      zRes = np.append( zRes, zmean )
  if len(tRes) == 0 :
    print 'Running mean could not be computed, skipping (time series too short?)'
    return
  ta = timeArray.timeArray( tRes, 'epoch' )
  data = vRes.reshape( (1,1,-1) )
  if dc.xDependsOnTime : x = xRes[None,:]
  if dc.yDependsOnTime : y = yRes[None,:]
  if dc.zDependsOnTime : z = zRes[None,:]
  meta = dc.getMetaData()
  dc2 = dataContainer( '', ta, x,y,z, data,
                            dc.fieldNames[:1], coordSys=dc.coordSys,metaData=meta)
  return dc2

def runningMean(dc, T=T_M2, gap_dt=None):
  """Computes running mean of the time series in dataContainer dc, using
  window length T. Returns new time series in dataContainer, one data point
  per window.
  """
  return runningX(dc, T, operator=computeRunningMean)

def runningMin(dc, T=T_M2, gap_dt=None):
  """Computes running min of the time series in dataContainer dc, using
  window length T. Returns new time series in dataContainer, one data point
  per window.
  """
  return runningX(dc, T, operator=computeRunningMin, gap_dt=None)

def runningMax(dc, T=T_M2, gap_dt=None):
  """Computes running max of the time series in dataContainer dc, using
  window length T. Returns new time series in dataContainer, one data point
  per window.
  """
  return runningX(dc, T, operator=computeRunningMax, gap_dt=None)

def runningRange(dc, T=T_M2, gap_dt=None):
  """Computes running range of the time series in dataContainer dc, using
  window length T. Returns new time series in dataContainer, one data point
  per window.

  Range is defined as x_window.max()-x_window.min() for each window.
  """
  return runningX(dc, T, operator=computeRunningRange)

def removeTides( dc, dt=None, gapFactor=20, T=T_M2 ) :
  """A low-pass filter to remove tidal signal from the data"""
  from scipy import signal
  time = dc.time.array
  vals = dc.data[0,0,:]
  if dt==None :
    dt = np.diff(time).mean()
  # try to calculate exact dt by omitting large gaps
  gaps,ranges,t = dc.detectGaps(dt=dt,gapFactor=gapFactor)
  diff = []
  for i in range( ranges.shape[0] ) :
    twin = time[ ranges[i,0]:ranges[i,1] ]
    diff.append(np.diff(twin))
  diff = np.concatenate( tuple(diff), axis=0 )
  dt = diff.mean()
  # filter design, low-pass butterworth
  T0 = (2*dt) # period of Nyquist frequency
  Tpass = 8*T # period of pass frequency
  Gpass = 3.0       # max dB loss in pass band 
  Tstop = 1*T # period of stop frequency
  Gstop = 30.0     # min dB atennuation in stop band
  o,Wn = signal.buttord(T0/Tpass,T0/Tstop,Gpass,Gstop)
  if o < 0 :
    raise Exception('Cannot create tidal filter. Data sampling frequency may be too low, dt='+str(dt))
  b,a = signal.butter(o,Wn,'low')
  newvals = []
  newtime = []
  # filter each contiquous data range separately
  for i in range( ranges.shape[0] ) :
    twin = time[ ranges[i,0]:ranges[i,1] ]
    vwin = vals[ ranges[i,0]:ranges[i,1] ]
    if len(vwin) > 3*len(a) :
      try :
        # default forward-backward filter
        # filtered = signal.filtfilt(b, a, vwin, padtype='constant')
        # forward-backward filter with custom boundary conditions
        # pad with mean of 1/2 pass window lenght
        N_init = int(np.ceil(Tpass/dt/2))
        # forward filter
        x_init = vwin[:N_init]
        y_init = x_init.mean()*np.ones_like(x_init)
        z_init = signal.lfiltic(b,a,y_init,x_init)
        filtered,_ = signal.lfilter(b, a, vwin, zi=z_init)
        # backward filter
        x_init = vwin[-N_init:][::-1]
        y_init = x_init.mean()*np.ones_like(x_init)
        z_init = signal.lfiltic(b,a,y_init,x_init)
        filtered,_ = signal.lfilter(b, a, filtered[::-1], zi=z_init)
        filtered = filtered[::-1]
        newvals.append( filtered )
        newtime.append( twin )
      except Exception as e :
        print a.shape, vwin.shape
        raise e
  newvals = np.concatenate( tuple(newvals), axis=0 )
  newtime = np.concatenate( tuple(newtime), axis=0 )
  dc2 = dc.copy()
  dc2.data = newvals[None,None,:]
  dc2.time.array = newtime
  return dc2

def smooth(time, vals, dt=None, gapFactor=20, T=T_M2):
  from scipy import signal
  if dt==None :
    dt = np.diff(time).mean()
  ta = timeArray.timeArray(time, 'epoch')
  # try to calculate exact dt by omitting large gaps
  gaps,ranges,t = ta.detectGaps(dt=dt,gapFactor=gapFactor)
  diff = []
  for i in range( ranges.shape[0] ) :
    twin = time[ ranges[i,0]:ranges[i,1] ]
    diff.append(np.diff(twin))
  diff = np.concatenate( tuple(diff), axis=0 )
  dt = diff.mean()
  # filter design, low-pass butterworth
  T0 = (2*dt) # period of Nyquist frequency
  Tpass = 8*T # period of pass frequency
  Gpass = 3.0       # max dB loss in pass band 
  Tstop = 1*T # period of stop frequency
  Gstop = 30.0     # min dB atennuation in stop band
  o,Wn = signal.buttord(T0/Tpass,T0/Tstop,Gpass,Gstop)
  if o < 0 :
    raise Exception('Cannot create tidal filter. Data sampling frequency may be too low, dt='+str(dt))
  b,a = signal.butter(o,Wn,'low')
  newvals = []
  newtime = []
  # filter each contiquous data range separately
  for i in range( ranges.shape[0] ) :
    twin = time[ ranges[i,0]:ranges[i,1] ]
    vwin = vals[ ranges[i,0]:ranges[i,1] ]
    if len(vwin) > 3*len(a) :
      try :
        # default forward-backward filter
        # filtered = signal.filtfilt(b, a, vwin, padtype='constant')
        # forward-backward filter with custom boundary conditions
        # pad with mean of 1/2 pass window lenght
        N_init = int(np.ceil(Tpass/dt/2/4))
        # forward filter
        x_init = vwin[:N_init]
        y_init = x_init.mean()*np.ones_like(x_init)
        z_init = signal.lfiltic(b,a,y_init,x_init)
        filtered,_ = signal.lfilter(b, a, vwin, zi=z_init)
        # backward filter
        x_init = vwin[-N_init:][::-1]
        y_init = x_init.mean()*np.ones_like(x_init)
        z_init = signal.lfiltic(b,a,y_init,x_init)
        filtered,_ = signal.lfilter(b, a, filtered[::-1], zi=z_init)
        filtered = filtered[::-1]
        newvals.append( filtered )
        newtime.append( twin )
      except Exception as e :
        print a.shape, vwin.shape
        raise e
  newvals = np.concatenate( tuple(newvals), axis=0 )
  newtime = np.concatenate( tuple(newtime), axis=0 )
  return newtime, newvals

def timeDerivative(dc):
    """Computes time derivative of a time series"""
    gaps,ranges,t = dc.detectGaps()
    v = dc.data[0,0,:]
    new_vals = []
    new_time = []
    for i in range(ranges.shape[0]):
        twin = t[ranges[i,0]:ranges[i,1]]
        vwin = v[ranges[i,0]:ranges[i,1]]
        dval = np.diff(vwin)
        dt = np.diff(twin)
        dvaldt = dval/dt
        new_t = 0.5*(twin[1:] + twin[:-1])
        new_vals.append(dvaldt)
        new_time.append(new_t)
    new_vals = np.hstack(tuple(new_vals))
    new_time = np.hstack(tuple(new_time))
    ta = timeArray.timeArray(new_time, 'epoch')
    data = new_vals.reshape((1, 1, -1))
    x = dc.x
    y = dc.y
    z = dc.z
    meta = dc.getMetaData()
    dc2 = dataContainer('', ta, x, y, z, data,
                        dc.fieldNames[:1], coordSys=dc.coordSys, metaData=meta)
    return dc2

def detectSignChanges(dc, direction='both'):
    """Detects sign changes in time series.
    direction can be 'pos2neg', 'neg2pos', or 'both'
    """
    gaps,ranges,t = dc.detectGaps()
    v = dc.data[0,0,:]
    new_vals = []
    new_time = []
    for i in range(ranges.shape[0]):
        twin = t[ranges[i,0]:ranges[i,1]]
        vwin = v[ranges[i,0]:ranges[i,1]]
        pos = vwin > 0
        if direction == 'pos2neg':
            change_ix = (pos[:-1] & ~pos[1:]).nonzero()[0]
        elif direction == 'neg2pos':
            change_ix = (~pos[:-1] & pos[1:]).nonzero()[0]
        else: # both directions
            change_ix = ((~pos[:-1] & pos[1:]) | (pos[:-1] & ~pos[1:])).nonzero()[0]
        new_time.append(twin[change_ix])
        new_vals.append(vwin[change_ix])
    new_vals = np.hstack(tuple(new_vals))
    new_time = np.hstack(tuple(new_time))
    ta = timeArray.timeArray(new_time, 'epoch')
    data = new_vals.reshape((1, 1, -1))
    x = dc.x
    y = dc.y
    z = dc.z
    meta = dc.getMetaData()
    dc2 = dataContainer('', ta, x, y, z, data,
                        dc.fieldNames[:1], coordSys=dc.coordSys, metaData=meta)
    return dc2

def detectHighLowWater(dc):
    """
    Detects high/low water time stamps from elevation time series.

    Returns two dataContainers containing only hw/lw time stamps.
    """
    detadt = timeDerivative(dc)
    detadt = removeTides(detadt, T=T_M2/24.0)
    detadt_pos2neg = detectSignChanges(detadt, direction='pos2neg')
    detadt_neg2pos = detectSignChanges(detadt, direction='neg2pos')

    # detect all high/low water time stamps
    highWaterDC = dc.interpolateInTime(detadt_pos2neg.time)
    lowWaterDC = dc.interpolateInTime(detadt_neg2pos.time)

    return highWaterDC, lowWaterDC

def computeTidalRange(dc, T=T_M2):
    gaps,ranges,t = dc.detectGaps()
    x = np.squeeze(dc.data)
    tRes = np.array([])
    xRes = np.array([])
    for i in range( ranges.shape[0] ) :
      tt = t[ ranges[i,0]:ranges[i,1] ]
      xx = x[ ranges[i,0]:ranges[i,1] ]
      tt, xx = computeRunningRange( tt,xx,2*T )
      if len(tt) == 0 :
        continue
      tt, xx = computeRunningMean( tt,xx,3*T )
      tRes = np.append( tRes, tt )
      xRes = np.append( xRes, xx )
    if len(tRes) == 0 :
      print 'tidal data could not be computed, skipping (time series too short?)'
      return
    ta = timeArray.timeArray( tRes, 'epoch' )
    data = xRes.reshape( (1,1,-1) )
    meta = dc.getMetaData()
    meta['dataType'] = 'timeseries'
    meta['tag'] = dc.getMetaData('tag')
    meta['variable'] = 'tidal_range'
    dc2 = dataContainer( '', ta, dc.x,dc.y,dc.z, data,
                              ['tidal_range'], coordSys='',metaData=meta)
    return dc2
