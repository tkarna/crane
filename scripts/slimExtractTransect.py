#!/usr/bin/python
"""
Extract transect data using the efficient SLIM extract_mod python module.

Examples:

# extract vdff, -d data dir, -t defines transect bp file, -n transect name string, -o output dir, -s -e time range
python slimExtractTransect.py -d ~pturner/db29/run29/outputs/ -v vdff -t nchannel_fine.bp -n nchannel -o run29/transect -s 2012-5-1 -e 2012-5-17

Tuomas Karna 2012-11-16
"""

import numpy as np
import sys
import datetime

from crane.physicalVariableDefs import addTracers
import crane.data.slimNcExtract as SNE

#-------------------------------------------------------------------------------
# Main: Commandline interface
#-------------------------------------------------------------------------------
def parseCommandLine() :
  from optparse import OptionParser

  parser = OptionParser()
  parser.add_option('-r', '--runTag', action='store', type='string',
                      dest='runTag', help='Run tag, used as a label in post-proc.')
  parser.add_option('-d', '--dataDirectory', action='store', type='string',
                      dest='dataDir', help='directory where model outputs are stored')
  parser.add_option('-s', '--start', action='store', type='string',
                      dest='startStr', help='Date to start processing')
  parser.add_option('-e', '--end', action='store', type='string',
                      dest='endStr', help='Date to end processing')
  parser.add_option('', '--stacks', action='store', type='string',
                    dest='stackStr', help='range of output files to read '
                    '(e.g 1,14) if start,end not given')
  parser.add_option('-v', '--variable', action='store', type='string',
                      dest='varList', help='variable(s) to extract: elev,temp,salt, ...\nTo use specific output file define extrension, e.g. salt.70')
  parser.add_option('-t', '--buildPointFile', action='store', type='string',
                      dest='bpFile', help='text file (*.bp) containing the transect (x,y) coordinates',
                             default=None)
  parser.add_option('-n', '--name', action='store', type='string',
                      dest='name', help='name of the transect (Default: first line in *.bp file)')
  parser.add_option('-c', '--modelCoordSys', action='store', type='string',
                      dest='modelCoordSys', default='spcs',
                      help='horizontal coordinate system used in model: '
                           'spcs or utm (Default: %default)')
  parser.add_option('-o', '--outDirectory', action='store', type='string',
                      dest='outDir', help='base directory for netCDF file tree '
                                                '(optional)')
  parser.add_option('-T', '--tracerModel', action='store', type='string', dest='tracerModel',
                    help='Enable extraction of tracers: sed, oxy, generic. Must '
                         'supply number of tracers for \'sed\' and \'generic\' '
                         'models via the -N switch. \'oxy\' model provides tracers: '
                         '\'NO3\',\'NH4\',\'phy\',\'zoo\',\'det\' and \'oxy\'.', default=None)
  parser.add_option('-N', '--numTracers', action='store', type='int', dest='numTracers',
                    help='Number of tracers to support for \'sed\' and \'generic\' models',
                    default=None)
  parser.add_option('', '--decimals', action='store', type='int', dest='digits',
                    help='Round extracted data to given decimal precision to save disk space', default=None)

  (options, args) = parser.parse_args()

  dataDir       = options.dataDir
  varList       = options.varList.split(',') if options.varList else None
  bpFile        = options.bpFile
  name          = options.name
  modelCoordSys = options.modelCoordSys
  outDir        = options.outDir
  runTag        = options.runTag
  startStr      = options.startStr
  endStr        = options.endStr
  stackStr      = options.stackStr
  runTag        = options.runTag
  tracerModel   = options.tracerModel
  numTracers    = options.numTracers
  digits        = options.digits

  if tracerModel :
    if not numTracers and tracerModel.split('.')[0] in ['sed','generic']:
      parser.print_help()
      parser.error('numTracers must be provided if sed or generic tracer models are used.')
    addTracers( tracerModel, varList, numTracers=numTracers)

  if not dataDir :
    parser.print_help()
    parser.error('dataDir  undefined')
  if not varList and tracerModel is None :
    parser.print_help()
    parser.error('variable undefined')
  #if not outDir :
    #parser.print_help()
    #parser.error('outDir   undefined')
  if startStr is None and stackStr is None:
    parser.print_help()
    parser.error('startStr undefined')
  if endStr is None and stackStr is None:
    parser.print_help()
    parser.error('endStr   undefined')
  if not name :
    parser.print_help()
    parser.error('name  undefined')
  if not runTag :
    parser.print_help()
    parser.error('runTag  undefined')
  if not bpFile :
    parser.print_help()
    parser.error('buildPointFile  undefined')

  if stackStr is not None:
    limits = [int(v) for v in stackStr.split(',')]
    stacks = np.arange(limits[0], limits[1]+1)
  else:
    stacks = None

  if startStr is not None:
    startTime = datetime.datetime.strptime( startStr ,'%Y-%m-%d')
    endTime = datetime.datetime.strptime( endStr ,'%Y-%m-%d')
  else:
    startTime = endTime = None

  print 'Parsed options:'
  if name :
    print ' - name ', name
  if stackStr is None:
    print ' - time range:',str(startTime),'->', str(endTime)
  else:
    print ' - stacks:',str(stacks[0]),'->', str(stacks[-1])
  print ' - dataDir',dataDir
  print ' - runTag',runTag
  print ' - transect file',bpFile
  if outDir :
    print ' - output dir',outDir
  print ' - variables ', varList
  print ' - model coord system', modelCoordSys
  sys.stdout.flush()

  se = SNE.slimExtract(dataDir, varList[0], verbose=True)
  for i in range(len(varList)): 
    dcs = se.extractTransectForBPFile(startTime, endTime, varList[i], bpFile, name)
    dcs.setMetaData( 'tag', runTag )
    import crane.data.dirTreeManager as dtm
    rule = 'singleFile'
    dtm.saveDataContainerInTree( dcs, rootPath=outDir, dtype=np.float32, overwrite=True, compress=True,
                                 digits=digits, rule=rule )

if __name__=='__main__' :
  parseCommandLine()
