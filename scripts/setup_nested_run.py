"""
One-way nested run setup

Sets up and creates files required to launch a nested run.  Essentially acts
as a wrapper around Fortran codes that actually do all the work.

lopezj - 08/06/2012
"""

#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import os
import shutil
from subprocess import call
import numpy as np
from optparse import OptionParser
from crane.files import grid
from crane.files import paramParser

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
BIN_PATH = '/home/workspace/users/lopezj/bin'

#-------------------------------------------------------------------------------
# Classes and functions
#-------------------------------------------------------------------------------
def genBoundary(sPath, ignore):
    """ Generates the boundary for creating *3D.th forcing files for nested run """
    oldWD = os.getcwd()
    # Go to directory with smaller domain hgrid.gr3
    os.chdir(sPath)

    # Create required gen_fg.in file
    makeFgIn(sPath, ignore)

    # Call Fortran program to actually create the fg.bp file
    gen_fg = "%s/gen_fg" % BIN_PATH
    call([gen_fg])      
 
    # Back to directory where you started
    os.chdir(oldWD)
    
def makeFgIn(sPath, ignore):
    """ Creates gen_fg.in file required as input for gen_fg program """
    file = os.path.join(sPath, 'hgrid.gr3')
    gridObj = grid.readHGrid(file)
    # If open boundaries are to be ignored, get them from list and create
    # a list of boundaries to be included for creation of boundary files.
    # Otherwise, all open boundaries will be forced.
    if ignore != []:
        nOB = gridObj.obNodes.nB - len(ignore)
        bndList = []
        for bnd in range(1,gridObj.obNodes.nB+1):
            if bndList.count(bnd):
                continue
            bndList.append(bnd)
    else:
        nOB = gridObj.obNodes.nB
        bndList = range(1,nOB+1)

    try:
        f = open("gen_fg.in","w")
        f.write("%d ! Total number of open boundary segments to be included in fg.bp\n" % nOB)
        for bnd in bndList:
            f.write("%d " % bnd)
        f.write(" ! List of open boundary IDs\n")
        f.close()
    except IOError:
        raise

def create3DFiles(lPath, sPath, outputsPath, nDays):
    """ Wrapper for interpolate_variables_selfe4.f90 which creates *3D.th files """
    oldWD = os.getcwd()
    # Must have permissions to write in this directory or there is troubletown
    os.chdir(outputsPath)

    # Copy required files needed as input 
    # Complaining if the files are the same, will not just overwrite...
    try:
      shutil.copy("%s/hgrid.gr3" % lPath, "bg.gr3")
    except:
      pass
    try:
      shutil.copy("%s/fg.bp" % sPath, "fg.bp")
    except:
      pass
    try:
      shutil.copy("%s/vgrid.in" % sPath, "vgrid.fg")
    except:
      pass

    # Create interpolate_variables_selfe.in and run program
    # Elevation = 1
    print 'Creating elev3D.th'
    makeInterpIn(nDays, 1)
    gen_fg = "%s/interpolate_variables_selfe4" % BIN_PATH
    call([gen_fg])  
    shutil.move('elev3D.th', 'elev900.th')

    # Salt and temp = 2
    print 'Creating salt3D.th and elev3D.th'
    makeInterpIn(nDays, 2)
    gen_fg = "%s/interpolate_variables_selfe4" % BIN_PATH
    call([gen_fg])
    shutil.move('temp3D.th', 'temp900.th') 
    shutil.move('salt3D.th', 'salt900.th') 

    # UV = 3
    #makeInterpIn(nDays, 3)
    #gen_fg = "%s/interpolate_variables_selfe4" % BIN_PATH
    #call([gen_fg])

    # Back to where we came from
    os.chdir(oldWD)

def makeTimeIn(sPath, inFile, nDays, timestep, ignore=[]):
    """ Creates input file for timeint.in code to interp *3D.th files """
    file = os.path.join(sPath,'hgrid.gr3')
    gridObj = grid.readHGrid(file)

    # If open boundaries are to be ignored, get them from list and create
    # a list of boundaries to be included for creation of boundary files.
    # Otherwise, all open boundaries will be forced.
    if ignore != []:
        nOB = gridObj.obNodes.nB - len(ignore)
        bndList = []
        for bnd in range(1,gridObj.obNodes.nB+1):
            if ignore.count(bnd):
                continue
            bndList.append(bnd)
    else:
        nOB = gridObj.obNodes.nB
        bndList = range(1,nOB+1)

    # Get number of nodes from open boundaries that will be used for forcings
    # nNodesOB = np.zeros((nOB))
    # for bnd in range(len(bndList)): 
    #nNodesOB[bnd] = bndList[bnd]

    # Need to symbolic link to original *3D.th file
    if os.path.lexists("th.old"):
        os.unlink("th.old")
    os.symlink(inFile, "th.old")

    # Need to determine number of vertical levels
    # Assumes files are named /blah/blah/????3D.th
    # Also assumes that vgrid.in exists in sPath directory
    if inFile[-9:-5] == 'elev':
        nLevels = 1
    else:
        vf = open("%s/vgrid.in" % sPath)
        nLevels = int(vf.readline().rstrip("\n").split()[0])
        vf.close()
        if inFile[-9:-5] == 'hvel':
            nLevels = nLevels * 2

    # Actually create the file
    f = open("timeint.in","w")
    f.write("%s %s\n" % (nDays, nLevels))
    f.write("%d " % nOB)
    for bnd in range(nOB):
        nodes = eval('len(grid.obNodes.ob_%d_nodes)' % (bnd+1))
        print 'grid.obsNodes.ob_%d_nodes' % (bnd+1)
        f.write(" %d " % nodes) 
        print bnd, nodes 
    f.write("\n")
    f.write("%d\n" % timestep)
    f.close()

def interpTime3D(sPath, outputsPath, nDays, timestep, ignore):
    """ Wrapper for interpolating *.th files to a timestep """
    old_cd = os.getcwd()
    os.chdir(outputsPath)
    time_int = "%s/timeint_3Dth2" % BIN_PATH

    files = ['elev900.th', 'salt900.th', 'temp900.th']
    for f in files:
      inFile = os.path.join(outputsPath, f)
      makeTimeIn(sPath, inFile, nDays, timestep, ignore)
      call([time_int])
      # Move to nested run path
      shutil.move('th.new', '%s/%s3D.th' % (sPath, f[:4]))

    os.chdir(outputsPath)
    # Move newly interpolated file to an apprpriate name

def makeInterpIn(nDays, int):
    """ Makes input file for interpolation.
  
    Args:
      nDays -- Number of days to interpolate
      int -- 1. elev3D.th 2. salt3D.th and temp3D.th 3. uv3D.th
    """
    f = open('interpolate_variables_selfe.in','w')
    print int, nDays
    f.write('%s  %s\n' % (int, nDays))
    f.close()

def setupNestedRun(lPath, sPath, timestep, outputsPath, nDays, ignore):
    genBoundary(sPath, ignore)
    create3DFiles(lPath, sPath, outputsPath, nDays)
    # If timestep of new run doesn't match old run, interpolate
    param = paramParser.ParamParser("%s/param.in" % lPath)
    interpTime3D(sPath, outputsPath, nDays, timestep, ignore)

#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
if __name__ == '__main__':
# Create a parser for command line args and options
    import sys
    usage = ("Usage: %prog largePath smallPath outputPath nDays [options]")

    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--time", action="store", type="int",
                      dest="timestep", help="Time step of smaller run")
    parser.add_option("-i", "--ignore", action="store", type="string",
                      dest="th_exist", help="Do not create forcings on these boundary numbers in the small domain")

    parser.set_defaults(timestep = 90)
    parser.set_defaults(ignore = [])

    (options, args) = parser.parse_args()

    timestep = options.timestep
    ignore = options.ignore

# Grab paths for larger domain (lPath) and smaller domain (sPath) 
    if len(args) != 4:
        parser.error("Incorrect number of arguments")
    else:
        lPath = os.path.abspath(sys.argv[1])
        sPath = os.path.abspath(sys.argv[2])
        oPath = os.path.abspath(sys.argv[3])
        nDays = int(sys.argv[4])
    setupNestedRun(lPath, sPath, timestep, oPath, nDays, ignore) 
