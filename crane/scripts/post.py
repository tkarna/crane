"""
Post process a model run
"""
#-------------------------------------------------------------------------------
# imports  
#-------------------------------------------------------------------------------
import os
import time
import glob
import thread
import datetime

import Queue
from optparse import OptionParser

from files.pathChecker import checkPaths
from data.runInfo import RunInfo
from data import obsData 
from data import modelData
from data.combine import ThreadedCombine
from data.modelExtractor import ModelExtractor 
from data.offeringManager import OfferingManager
from data.stationSet import * 
from plotting.plot import Plots 
#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
vars = ['elev.61','hvel.64','pres.61','airt.61','shum.61','srad.61','flsu.61',
        'fllu.61','radu.61','radd.61','flux.61','evap.61','prcp.61','wind.62',
        'wist.62','dahv.62','vert.63','temp.63','salt.63','conc.63','tdff.63',
        'vdff.63','kine.63','mixl.63','zcor.63','qnon.63','totN.63','Hsig.61',
        'WavD.61']

#-------------------------------------------------------------------------------
# Clases and functions 
#-------------------------------------------------------------------------------
class Post(object):
  """Class to organize post-processing of a model run.

  Attributes:
    paths -- Dictionary of strings to run directories
    vars -- List of variables to combine
    runInfo -- RunInfo object of run info
    stationsFile -- String of path to stations.sta file
    cQueue -- Queue of output files to combine
    eQueue -- Queue of days to extract station data
    startTime -- Datetime of time to start processsing
    endTime -- Datetime of time to end processing (inclusive)
  """
  nThreads = 1

  def __init__(self, paths, stationsFile, startDay=None, endDay=None):
    """Prepares post-processing of a model run.
    
    Args:
      paths -- Dictionary of strings to run directories
      stationsFile -- String of path to stations.sta file
      startDay -- Intger of stack (day) to start processing
      endDay -- Integer of stack (day) to end processing (inclusive)
    """
    st = time.time()
    self.paths = checkPaths(paths)
    if not os.path.isfile(stationsFile):
      raise Exception('Station file does not exist: '+stationsFile) 
    self.stationsFile = stationsFile
    self.runInfo = RunInfo(self.paths['run'])
    if startDay != None: 
      self.startDay = startDay
      self.startTime = self.runInfo.startTime + datetime.timedelta(startDay-1)
    else:
      self.startDay = 1
    if endDay != None: 
      self.endDay = endDay
      self.endTime = self.runInfo.startTime + datetime.timedelta(endDay)
    else:
      self.endDay = self.runInfo.param['rnday']
    self.vars = vars
    self.makeCombineQueue()
    self.eQueue = Queue.Queue()
    self.pQueue = Queue.Queue()

    # ----- Combine Files -----
    for i in xrange(self.nThreads):
      combine = ThreadedCombine(self.cQueue, self.paths['outputs'], self.eQueue)
      combine.setDaemon(True)
      combine.start()

    # ----- Extract stations -----
    thread.start_new_thread(self.extractStations, ())

    # Wait until all days are extracted? or Create daily plots? Wait for now.
    # If it is daily, then file could be used for getObsData.
    self.cQueue.join()
    for i in xrange(self.nThreads):
      combine.join()
    print 'Done combining'
    self.eQueue.join()
    print 'Done extracting'

    # ----- Load data -----
    print 'Loading obs and model data'
    try:
      oData = obsData.loadByTimes(self.paths['obs'], self.startTime, 
                                  self.endTime)
    except:
      oData = obsData.loadFromDatabase(self.startTime, self.endTime,
        self.paths['process'], obsPath=self.paths['obs'])
    try:
      mData = modelData.loadByTimes(self.paths['process'], self.startTime,
                                    self.endTime)
    except:
      mData = modelData.loadFromExtractedFiles(self.startDay, self.endDay,
                                               self.paths['process'])

    print 'Plotting'
    plot = Plots(self.paths['images'], oData, mData)
    plot.makeTimeSeries()
    plot.makeTaylorDiagrams()
    plot.makeTaylorDiagramsVar()

    print 'Elapsed time with %d threads: %f' % (self.nThreads, time.time()-st)

  def makeCombineQueue(self):
    """Generates a queue and dictionary of files to combine"""
    self.varCheck()
    q = Queue.Queue()
    d = {}
    for i in xrange(self.startDay, self.endDay+1):
      for v in self.vars:
        q.put((i,v))
        d[(i,v)] = 1
    print 'Number of vars to combine: '+str(q.qsize())
    print 'Combine queue ready'
    self.cQueue = q
    self.vars = d 

  def varCheck(self):
    """Checks VARS list and adds tracers if needed."""
    if self.startDay != None:
      file = os.path.join(self.paths['outputs'], '%d_0000_elev.61' % self.startDay)
    else:
      file = os.path.join(self.paths['outputs'], '1_0000_elev.61')
    while not os.path.isfile(file):
      time.sleep(60)
  
    # ----- Only output vars -----
    files = glob.glob('%s/%s_0000_[a-z]*.6?' % 
                      (self.paths['outputs'], self.startDay))
    outVars = [x.split('_')[2] for x in files]
    self.vars = set(outVars) & set(self.vars)

    # ----- Add tracers -----
    for t in xrange(1, int(self.runInfo.param.params['ntracers'])+1):
      self.vars.append('trcr_%d.63' % t)

    print 'Vars ready'

  def extractStations(self):
    """When all files for a day are combined, extracts station data"""
    while len(self.vars) > 0:
      # Timeout of 3600 there to allow for keyboard interrupt. 
      v = self.eQueue.get(True,3600)  
      print 'From eQueue.get: '+str(v)
      del self.vars[v]
      if len(self.vars) > 0:
        for key in self.vars:
          print 'From self.vars key: '+str(key)
          print 'Vars left:'
          print self.vars
          print 'key[0] '+str(key[0])
          print 'v[0] '+str(v[0])
          if key[0] == v[0]:
            print 'Still have combining to do for day '+str(key[0])
            break
          else:
            self.startExtractThread(v[0])
      else:
        self.startExtractThread(v[0])
      self.eQueue.task_done()

  def startExtractThread(self, day):
    """Starts new thread to extract stations for one day"""
    print 'Extracting day '+str(day)
    mE = ModelExtractor(self.paths['outputs'],self.paths['process'],
                        self.stationsFile, startDay=day, endDay=day)
    self.extract = thread.start_new_thread(mE.extractModelStations, ())

#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
if __name__ == '__main__':
  """Extracts model and obs data, performs skill assessment, and makes plots."""
  # Args
  usage = ('Usage: %prog [Path to run directory] [Path to station.sta file]\n')
  parser = OptionParser(usage=usage)
  parser.add_option('-s', '--startDay', action='store', type='int',
                    dest='startDay', help='Day to start extraction, combining, '
                    ' and plotting.')
  parser.add_option('-e', '--endDay', action='store', type='int',
                    dest='endDay', help='Day to stop extraction, combining, '
                    ' and plotting. (Inclusive)')

  parser.set_defaults(startDay = None)
  parser.set_defaults(endDay = None)
  (options, args) = parser.parse_args()

  if len(args) != 2:
    parser.error(usage)
  path = args[0]
  stationsFile = args[1]
 
  # Let's post-process... 
  post = Post(path, stationsFile, options.startDay, options.endDay)
