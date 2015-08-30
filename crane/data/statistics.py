#!/usr/bin/python

"""
Implementation of basic statistical measures.

Tuomas Karna
"""

import numpy as np

def computeSkillStatistics(reference, sample) :
  res = {}
  res['bias'] = bias(reference, sample)
  res['ioa'] = indexOfAgreement(reference, sample)
  res['murphy'] = murphySkillScore(reference, sample)
  res['stddev'] = standardDeviation(sample, ddof = 0)
  res['crmse'] = rootMeanSquareError(reference, sample, centered=True)
  res['rmse'] = rootMeanSquareError(reference, sample, centered=False)
  res['corr'] = correlationCoefficient(reference, sample)
  res['nmse'] = normalizedMeanSquareError(reference, sample)
  res['mse'] = meanSquareError(reference, sample, centered=False)
  return res

class skillAssessmentStats(object) :
  """A class for computing statistical measures for skill assessment. """
  def __init__(self, reference, refLabel='reference') :
    """Creates a new object with the given reference signal.
    The refLabel is used to identify the reference signal."""
    self.reference = reference
    self.refLabel = refLabel
    self.stats = dict()
    self.refStats = self.addSample( reference, refLabel )
    
  def addSample( self, sample, label, reference=None ) :
    """Computes statistics of the given sample versus the reference signal
    and appends the results to the statistics. The sample is identified by the
    given keyword label (typically a string). Note that the sample sample signal
    must be compatible with the reference so that the statistics make sense.
    """
    #self.stats[(label,'bias')] = bias( self.reference, sample )
    #self.stats[(label,'ioa')] = indexOfAgreement( self.reference, sample )
    #self.stats[(label,'stddev')] = standardDeviation( sample, ddof = 0 )
    #self.stats[(label,'crmse')] = rootMeanSquareError( self.reference, sample, centered=True )
    #self.stats[(label,'corr')] = correlationCoefficient( self.reference, sample )
    if reference is None :
      reference = self.reference
    res = computeSkillStatistics( reference, sample )
    for key in res:
        self.stats[(label,key)] = res[key]
    return res
  
  # emulate basic dictionary methods
  def __getitem__( self, key ) :
    """Return a value of the stats dictionary"""
    return self.stats[key]
    
  def __setitem__(self, key, value) :
    """Assign value to stats dictionary"""
    self.stats[key] = value
    
  def __iter__(self) :
    """Return iterator for the stats dictionary"""
    return self.stats.__iter__()
    
  def keys(self) :
    return self.stats.keys()
    
  def getStatistics(self) :
    """
    Returns a dictionary with all the computed statistics.
    Keys are (label, statStr) tuples, where statStr is one of the following:
    'bias'   -- bias
    'ioa'    -- index of agreement
    'stddev' -- standard deviation 
    'crmse'   -- centered root mean square error
    'corr'   -- correlation coefficient
    """
    return self.stats
    
  def getStatisticsArray(self) :
    """
    Returns an array with all the computed statistics,
    a list of all labels and a list of the variable names.
    """
    nSamples = len(self.stats)/5
    arr = np.zeros((nSamples,5))
    labels = self.getLabels()
    names = self.getStatisticsNames()
    
    for i,l in enumerate(labels) :
      for j,n in enumerate(names) :
        arr[i,j] = self.stats[(l,n)]
    
    return arr, labels, names
    
  def getStatisticsNames(self) :
    """Returns a list of all statistical measures used"""
    return ['corr','stddev','bias','ioa','crmse']
    
  def getLabels(self) :
    """Returns a list of all the labels. First one is the reference."""
    labels = set()
    for key in self.keys() :
      if key[0] != self.refLabel :
        labels.add( key[0] )
    labels = list(labels)
    labels = [self.refLabel] + labels
    return labels
    
class skillAssessmentStatsDC(skillAssessmentStats) :
  """Statistics class that uses dataContainer objects"""
  def __init__(self, reference, refLabel='reference') :
    self.checkData( reference )
    self.referenceDC = reference
    skillAssessmentStats.__init__(self, np.squeeze(reference.data), refLabel)
    
  def addSample( self, sample, label ) :
    """Computes statistics of the given sample versus the reference signal
    and appends the results to the statistics. The sample is identified by the
    given keyword label (typically a string). If the two dataContainers are not
    compatible, an error is raised.
    """
    self.checkData( sample )
    skillAssessmentStats.addSample(self, np.squeeze(sample.data), label)
   
  @staticmethod
  def checkData( sample ) :
    if not isinstance( sample, dataContainer ) :
      raise Exception( 'sample must be a dataContainer object' )
    if sample.data.shape[0] > 1 :
      raise Exception( 'spatially varying data is not supported' )
    if sample.data.shape[1] > 1 :
      raise Exception( 'multiple fields are not supported' )
    return True

def bias( reference, sample ) :
  """Computes the bias of sample versus reference signal
  bias = mean( sample - reference )"""
  bias = np.mean( sample - reference )
  return bias
  
def indexOfAgreement( reference, sample ) :
  """Computes the index of agreement between the sample and the reference signal"""
  refMean = reference.mean()
  a = sum( ( sample - reference )**2 )
  b = sum( ( abs(sample-refMean) + abs(reference-refMean) )**2 )
  if b > 1e-20 :
    ioa = 1 - a/b
  else :
    ioa = 1
  return ioa

def murphySkillScore(reference, sample, reference_model=None) :
    """
    Computes the skill metric defined in Murpy (1988)
    Monthy Weather Review paper.
    """
    # here reference forecast is assumed to be the mean of observations
    if reference_model is None:
        reference_model = reference.mean()*np.ones_like(reference)
    mse_model = np.mean((sample - reference)**2)
    mse_refmodel = np.mean((reference_model - reference)**2)
    if mse_refmodel > 0.0:
        score = 1 - mse_model/mse_refmodel
    else:
        score = -1
    return score

def normalizedMeanSquareError(reference, sample):
    """Computes the mean square error divided by the variance of the reference
    signal."""
    # here reference forecast is assumed to be the mean of observations
    reference_model = reference.mean()
    mse_model = np.mean((sample - reference)**2)
    mse_refmodel = np.mean((reference_model - reference)**2)
    if mse_refmodel > 0.0:
        score = mse_model/mse_refmodel
    else:
        score = 2.0
    return score

def standardDeviation( sample, ddof = 0 ) :
  """Computes standard deviation of the signal.
  Use ddof = 1 to compute std with N-1 in the denominator (i.e. the unbiased estimator)
  ddof=0 by default."""
  stdev = np.std(sample,ddof=ddof)
  return stdev
  
def correlationCoefficient( reference, sample, ddof = 0 ) :
  """Computes the correlation coefficient between the sample and the reference signal.
  Use ddof = 1 to compute std with N-1 in the denominator (i.e. the unbiased estimator)
  ddof=0 by default."""
  bias = int(not ddof)
  corr = np.corrcoef(reference, sample, bias=bias)[0,1]
  #corr = np.corrcoef(reference, sample, ddof=ddof)[0,1]
  return corr

def meanSquareError( reference, sample, centered=True ) :
  """Computes mean squre error.
  rmse = mean( ( r - s )**2 )
  If centered=True (default), mean is removed from signals before computing the error.
  rmse = mean( ( (r-r.mean()) - (s-s.mean()) )**2 )
  """
  err = reference - sample
  if centered == True :
    err += -reference.mean() + sample.mean()
  mse = np.mean( (err)**2 )
  return mse

def rootMeanSquareError( reference, sample, centered=True ) :
  """Computes root mean squre error.
  rmse = sqrt( mean( ( r - s )**2 ) )
  If centered=True (default), mean is removed from signals before computing the error.
  rmse = sqrt( mean( ( (r-r.mean()) - (s-s.mean()) )**2 ) )
  """
  err = reference - sample
  if centered == True :
    err += -reference.mean() + sample.mean()
  rmse = np.sqrt( np.mean( (err)**2 ) )
  return rmse
  
if __name__=='__main__':

  # brief sanity check
  x = np.linspace(0,30,100)
  ref = np.sin(x) + 0.95*np.sin(0.95*x)
  sam = 0.69*np.sin(x) + 0.2*np.sin(0.5*x+0.05) - 0.08
  
  # rmse**2 = bias**2 + crmse**2
  rmse = rootMeanSquareError(ref,sam, centered=False)
  crmse = rootMeanSquareError(ref,sam, centered=True)
  bi = bias(ref,sam)
  check = rmse - np.sqrt( crmse**2 + bi**2 )
  print 'RMSE check: ' + str(check)
  
  # crmse**2 = std(b)**2 + std(c)**2 - 2*std(b)*std(c)*corrCoef
  stdr = standardDeviation(ref,ddof=0)
  stds = standardDeviation(sam,ddof=0)
  corr = correlationCoefficient(ref,sam,ddof=0)
  check = crmse - np.sqrt(stdr**2 + stds**2 - 2*stdr*stds*corr)
  print 'Taylor check: ' + str(check)
  
  print 'bias',bi
  print 'rmse',crmse
  print 'ioa',indexOfAgreement( ref,sam )
  print 'corr',corr
  print 'stddev',stds
  
  # check stats class
  st = skillAssessmentStats(ref, 'ref')
  st.addSample(sam, 'sam')
  print st.getStatistics()