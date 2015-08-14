"""
Combines parallel SELFE output files
"""
#-------------------------------------------------------------------------------
# Imports 
#-------------------------------------------------------------------------------
import subprocess
import time
import os

import threading
import Queue
from runInfo import RunInfo
from optparse import OptionParser 
#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
BIN_PATH = "/home/users/yinglong/SELFE/combine_output5"
VARS = ('elev.61','hvel.64','pres.61','airt.61','shum.61','srad.61','flsu.61',
        'fllu.61','radu.61','radd.61','flux.61','evap.61','prcp.61','wind.62',
        'wist.62','dahv.62','vert.63','temp.63','salt.63','conc.63','tdff.63',
        'vdff.63','kine.63','mixl.63','zcor.63','qnon.63','totN.63','Hsig.61',
        'WavD.61')

#-------------------------------------------------------------------------------
# Classes and functions 
#-------------------------------------------------------------------------------
class Combine(object):
  """Combines output files from parallel Selfe model runs.

  Attributes:
      path  -- String of path to the directory with the output files
      stack -- Integer of stack (typically day) to combine
      var   -- String of a var to combine (e.g. 'salt.63', 'elev.61')
      procs -- Integer of the number of processors used in run. Number of
               partial files to combine.
  """
  def __init__(self, path, stack, var, lock=None): 
    """Creates Combine object.  Generates input file.

    Args:
      path  -- String of path to the directory with the output files
      stack -- Integer of stack (typically day) to combine
      var   -- String of a var to combine (e.g. 'salt.63', 'elev.61')
      lock  -- threading.Lock() to ensure combine.in is only written to by 1 
    Returns:
      NOTHING
    """
    if not os.path.isdir(path):
      raise Exception('Directory does not exist: '+path)

    if not var in VARS: 
      eStr = 'Var %s is not available in\n %s' % (var, VARS)
      raise Exception(eStr)

    self.path = path
    self.var = var
    self.stack = stack
    self.procs = self.readLocalToGlobal()[3]
    self.lock = lock

    while not self.isStackReady():
      print 'Stack is not ready. stack: %s - var: %s' % (stack, var)
      time.sleep(300)
    while not self.isVarReady():
      print 'Var is not ready. stack: %s - var: %s' % (stack, var)
      time.sleep(60)


  def combineEm(self):
    """Executes combine_output5 code to combine output files"""
    # Need to be in the output directory because combine_output5
    # requires the output files and local_to_global files to be there
    if self.lock != None:
      self.lock.acquire()
    self.makeInputFile()
    old = os.getcwd()
    os.chdir(self.path)

    args = [BIN_PATH]
    proc = subprocess.Popen(args, shell=False)
    if self.lock != None:
      self.lock.release()
    proc.wait()
  
  def isStackReady(self):
    """Checks if a stack is ready for combining"""
    nextStack = self.stack+1
    file = '%s/%d_0000_elev.61' % (self.path, nextStack) 
    ready = False
    if os.path.isfile(file): ready = True

    return ready

  def isVarReady(self):
    """Checks if specific variable is ready for combining"""
    for proc in xrange(self.procs):
      file = '%s/%d_%04d_%s' % (self.path, self.stack, proc, self.var)
      while not os.path.isfile(file):
        time.sleep(180)
    return True

  def checkLocalToGlobalFiles(self):
    """Checks for required local_to_global mapping files"""
    for p in xrange(self.procs):
      file = '%s/%d_%4d_%s' % (self.path, self.stack, p, self.var)
      if not os.path.isfile(file):
        raise Exception('Must have local_to_global mapping files in output '
                        'directory to combine partial files.')

  def makeInputFile(self):
    """Makes input file for use with combine_outputs5"""
    old = os.getcwd()
    os.chdir(self.path)

    file = self.path+'/combine_output.in'
    try:
      f = open(file,'w')
    except IOError:
      raise Exception('Unable to open input file for combine_outputs5: '+file)

    f.write('%s\n' % self.var)
    f.write('%d %d\n' % (self.stack, self.stack))
    f.write('0')
    f.close()
    os.chdir(old)

  def readLocalToGlobal(self, file='local_to_global_0000'):
    """Returns data about grid required for combine_output module

    Args:
      file -- String of path to a local_to_global file
    Returns:
      stats -- List of:
        [0] -- Integer of number of elements (ne_global)
        [1] -- Integer of number of nodes (np_global)
        [2] -- Integer of number of vertical levels (nvrt)
        [3] -- Integer of number of processors (nprocs)
        [4] -- Integer of number of tracers (ntracers)
    """
    file = os.path.join(self.path, file)
    if not os.path.isfile(file):
      raise Exception('local_to_global file does not exist: '+file)
    f = open(file,'r')
    stats = [int(x) for x in f.readline().rstrip('\n').split()]
      
    return stats



class ThreadedCombine(threading.Thread):
  """Threaded combiner of parallel Selfe outputs.

  !Do not use greater than one thread.  It is not threadsafe due to
  combine.in file.

  Gets tuple comboInfo from a queue.
  comboInfo[0] -- Integer of stack
  comboInfo[1] -- String of var (e.g. 'elev.61')
  
  Useful to speed up combining as long as disk I/O is not saturated.
  """
  def __init__(self, queue, path, oQueue=None):
    """
    Args:
      queue -- Queue of the variables to combine
      path -- String of path to the directory with output files
      procs -- Integer of the number of processors (# partial files)
      oQueue -- Can be used to pass args from queue to another queue
    """
    if not os.path.isdir(path):
      raise Exception('Directory to output files does not exist: '+path)
    super(ThreadedCombine, self).__init__() 
    self.queue = queue
    self.path = path
    self.oQueue = oQueue
    self.stopRequest = threading.Event()
    self.lock = threading.Lock()
  
  def run(self):
    """Combines a stack for a single var"""
    while not self.stopRequest.isSet(): 
      # Huge timeout specified to allow keyboard interrupt
      try:
        comboInfo = self.queue.get(timeout=10)
        stack = comboInfo[0]
        var = comboInfo[1]
        c = Combine(self.path, stack, var, self.lock)
        c.combineEm()
        if self.oQueue != None:
          self.oQueue.put(comboInfo)
        self.queue.task_done()
      except Queue.Empty:
        continue

  def join(self, timeout=None):
    self.stopRequest.set()
    super(ThreadedCombine, self).join(timeout)
  
  def readLocalToGlobal(self, file='local_to_global_0000'):
    """Returns data about grid required for combine_output module

    Args:
      file -- String of path to a local_to_global file
    Returns:
      stats -- List of:
        [0] -- Integer of number of elements (ne_global)
        [1] -- Integer of number of nodes (np_global)
        [2] -- Integer of number of vertical levels (nvrt)
        [3] -- Integer of number of processors (nprocs)
        [4] -- Integer of number of tracers (ntracers)
    """
    file = os.path.join(self.path, file)
    if not os.path.isfile(file):
      raise Exception('local_to_global file does not exist: '+file)
    f = open(file,'r')
    stats = [int(x) for x in f.readline().rstrip('\n').split()]
      
    return stats

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
def test_Combine():
  """Perfroms a test on combining files"""
  start = time.time()
  TEST_DIR = '/home/workspace/users/lopezj/data/test/combine'
  VARS = ['elev.61', 'wind.62', 'salt.63', 'hvel.64','elev.61', 'wind.62', 'salt.63', 'hvel.64',
          'elev.61', 'wind.62', 'salt.63', 'hvel.64','elev.61', 'wind.62', 'salt.63', 'hvel.64']
  STACK = 1
  PROCS = 24
  FILES = [('1_elev.61', 'test_elev.61'),
           ('1_wind.62', 'test_wind.62'), 
           ('1_salt.63', 'test_salt.63'),
           ('1_hvel.64', 'test_hvel.64')]
  print 'Combine files\n'
  for v in VARS:
    print 'Combining '+v
    c = Combine(TEST_DIR, STACK, v)
    c.combineEm()

  print 'Test created files against reference\n'
  good = True
  for f in FILES:
    print 'Comparing %s and %s' % (f[0], f[1])
    file1 = os.path.join(TEST_DIR, f[0])
    file2 = os.path.join(TEST_DIR, f[1])
    if os.path.getsize(file1) != os.path.getsize(file2):
      print 'Failed test, file sizes differ.'
      good = False

  if good:
     print 'Passed file combine test'
  else:
    print 'Failed file combine test'
  print 'Elapsed time: %s' % (time.time() - start)

def test_ThreadedCombine(nThreads=None):
  start = time.time()
  """Perfroms a test on combining files"""
  TEST_DIR = '/home/workspace/users/lopezj/data/test/combine'
  VARS = ['elev.61', 'wind.62', 'salt.63', 'hvel.64','elev.61', 'wind.62', 'salt.63', 'hvel.64',
          'elev.61', 'wind.62', 'salt.63', 'hvel.64','elev.61', 'wind.62', 'salt.63', 'hvel.64']
  STACK = 1
  PROCS = 24
  FILES = [('1_elev.61', 'test_elev.61'),
           ('1_wind.62', 'test_wind.62'), 
           ('1_salt.63', 'test_salt.63'),
           ('1_hvel.64', 'test_hvel.64')]
  if nThreads == None: 
    nThreads = len(VARS)

  queue = Queue.Queue()
  for var in VARS:
    queue.put((STACK,var))

  print '\nCombine files using %d threads\n' % nThreads
  for i in xrange(nThreads):
    t = ThreadedCombine(queue, TEST_DIR)
    t.setDaemon(True)
    t.start()
  queue.join() 
  for i in xrange(nThreads):
    t.join()

  print 'Test created files against reference\n'
  good = True
  for f in FILES:
    print 'Comparing %s and %s' % (f[0], f[1])
    file1 = os.path.join(TEST_DIR, f[0])
    file2 = os.path.join(TEST_DIR, f[1])
    if os.path.getsize(file1) != os.path.getsize(file2):
      print 'Failed test, file sizes differ.'
      good = False

  if good:
     print 'Passed file combine test'
  else:
    print 'Failed file combine test'
  print 'Elapsed time: %s' % (time.time() - start)

def test():
  """Performs a test on the classes.  """
  test_Combine()
  test_ThreadedCombine(4)

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
#TODO: Make thread safe!  Need to change input to the FORTRAN code! Or lock file
#TODO: Move readLocalToGlobal to a class
if __name__ == '__main__':
  test()    
