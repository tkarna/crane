"""
Verifies run related paths
"""
#-------------------------------------------------------------------------------
# Imports  
#-------------------------------------------------------------------------------
import os

#-------------------------------------------------------------------------------
# Functions  
#-------------------------------------------------------------------------------
def checkPaths(paths, forecast=False):
  """Validates supplied paths. Returns paths if good.
  Args:
    paths -- Dict of paths to directories listed above or String to run directory
    forecast -- Boolean to change format of run direcotry structure
  Returns:
    paths -- Dict of paths verified to exist
  """
  if type(paths) == str:
    if not os.path.isdir(paths):
      raise Exception('Run directory path does not exist: '+path)
    newPaths = {}
    newPaths['run'] = os.path.join(paths, 'run')
    newPaths['process'] = os.path.join(paths, 'process')
    newPaths['images'] = os.path.join(paths, 'images') 
    if forecast:
      newPaths['obs'] = os.path.join(paths, 'data')
      newPaths['outputs'] = newPaths['run']
    else:
      newPaths['obs'] = os.path.join(paths, '../obs')
      newPaths['outputs'] = os.path.join(newPaths['run'], 'outputs')
    paths = newPaths
  if type(paths) == dict:
    if not os.path.isdir(paths['run']):
      raise Exception('Run directory path does not exist: '+paths['run'])
    if not os.path.isdir(paths['obs']):
      raise Exception('Obs directory path does not exist: '+paths['obs'])
    if not os.path.isdir(paths['process']):
      raise Exception('Process directory path does not exist: '+paths['process'])
    if not os.path.isdir(paths['outputs']):
      raise Exception('Outputs directory path does not exist: '+paths['outputs'])
    if not os.path.isdir(paths['images']):
      raise Exception('Images directory path does not exist: '+paths['images'])
  else:
    raise Exception('Must provide path to run directory or specify directories')

  return paths

