#-------------------------------------------------------------------------------
# Imports 
#-------------------------------------------------------------------------------
import os
import re
import time
import datetime

from crane.files import paramParser
#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Classes and functions 
#-------------------------------------------------------------------------------
class RunInfo(object):
  """Gathers information about a model run
 
  Attributes:
    param -- Param object that has parameters from param.in
    nprocs -- Integer of the number of processors used (need outputs dir to exist)
    system -- String describing the cluster/computer used for model run (hostname)
    name -- Name of the model run (based on SGE submit script)
    startTime -- Datetime object of the start time of the run
    endTime -- Datetime of the end time of the run
    output_dt -- Integer of model output time step (dt*nspool)
  """
  def __init__(self, path, jobSubmitScript=None):
    """Creates an instance of RunInfo

    Args:
      path -- String of path to the run directory of the run. The path with the
        the binaries and input files. 
    """
    self.path = path
    file = path+'/param.in'
    if os.path.isfile(file):
      self.param = paramParser.ParamParser(file)
    else:
      raise Exception('File does not exist: '+param)

    dir = path+'/outputs/'
    #if os.path.isdir(dir):
    #  self.nprocs = self.getNumberOfCores()
    if jobSubmitScript != None: 
      self.name = self.getRunName()
    self.system = self.getSystemDescription()
    self.startTime = self.getStartTime()
    self.endTime = self.startTime + datetime.timedelta(int(self.param['rnday']))

  def getStartTime(self):
    """Returns the start time of a model run"""
    file = '%s/bctides.in' % self.path
    try:
      f = open(file)
      return self.DateStringToDatetime(f.readline())
    except IOError, e:
      raise Exception('Unable to open file: '+file)

  def getNumberOfCores(self):
  # TODO: This is broke. Fix it.
    """Returns the number of processors used for the run
    
    Based on the number of partial files for local_to_global_??? files. 

    Returns:
      len(match) -- Number of local_to_global_??? files as proxy for number of cores      
    """
    path = self.path+'/outputs'
    if not os.path.isdir(path):
      raise Exception('Output directory does not exist: '+path)
    dirList = os.listdir(path)
    pattern = re.compile("local_to_global")
    match = filter(pattern.match, dirList)
    # if files are not created yet sleep for a bit and wait
    if len(match) == 0:
      while True:
        time.sleep(60)
        if self.getNumberOfCores() != 0:
          nprocs = self.getNumberOfCores()
          break
    else:
      nprocs = len(match)
         
    return nprocs 
     
  def getSystemDescription(self):
    """Returns string description of the system this method is run on. 

    Args:
      NONE
    Returns:
      description - Just returns the hostname. 
    """
    self.system = os.uname()[1]
 
  def getRunName(self, jobSubmitScript):
    """Returns the name of the run from the SGE submit script.

    Args:
      jobSubmitScript -- String of the name of the SGE submit script
    Returns:
      runName -- String of the run name.  If none is found 'None' is returned.
    """
    name = None
    if not os.path.isfile(jobSubmitScript):
      raise Exception('SGE submit script not found: '+jobSubmitScript)
    else:
      file = open(jobSubmitScript)
      for line in file.readlines():
        tmp = line.split(' ')
        if tmp(' ')[0] == '#$' and tmp(' ')[1] == '-N':
          name = tmp(' ')[2]

    return name

  def DateStringToDatetime(self, dateString):
    """ Takes time in form 'mm/dd/yyyy HH:MM:SS PST' and returns datetime """
    time = dateString.split()
    month, day, year = time[0].split('/')
    year = int(year); month = int(month); day = int(day)
    hour, minute, second = time[1].split(':')
    hour = int(hour); minute = int(minute); second = int(second)

    return datetime.datetime(year, month, day, hour, minute, second)
 
