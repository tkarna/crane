"""
Cleans up and organizes a completed model run.

Written to be used as a stand alone script or Clean class can be imported to
be used as you please.

lopezj - 08/09/2012
"""

#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import os
import shutil
import glob
import tarfile
from optparse import OptionParser

#-------------------------------------------------------------------------------
# Classes and functions
#-------------------------------------------------------------------------------
class Clean(object):
	def __init__(self, path):
		self.path = path

	def filesToLog(self, path):
		""" Moves non-data files (not ?_[a-z]*.6?) from model to log directory """
		logDir = "%s/log" % path
		runDir = "%s/run" % path 
		try:
			os.mkdirs(path)
		except OSError, e:
			if e.errno != errno.EEXIST:
				raise

		# Files to move listed, then path added on
		filesToMove = ['centers.bp', 'coriolis.out', 'flux.out', 'mirror.out',
					   'sample_Z.out', 'sidecenters.bp', 'total.dat', 'total_ST.dat']
		for files in filesToMove:
			files = "%s/%s" % (runDir, files)

		filesToSearch = ["*.o[1-9]*","*.po[1-9]*"]
		for part in filesToSearch:
			for files in glob.glob("%s/%s" % (runDir, part)):
				filesToMove.append(files)

		# Loop over the files listed above and move to log dir
		for each in filesToMove:
			try:
				shutil.move(each, logDir)
			except:
				print "Unable to move %s to %s as part of post run cleanup" % (each, logDir)


	def makeRunArchive(self, path):
		""" Moves all input files for the run into 'input' directory and creates archive """
		inputDir = "%s/run/inputs" % path
		runDir = "%s/run" % path
		logDir = "%s/log" % path
		try:
			os.mkdirs(path)
		except IOError, e:
			if e.errno != errno.EEXIST:
				raise

		# ID all files to move in a list
		filesToMove = []
		filesToSearch = ["*.gr3","*.in","*.ic","*.th","pelfe*","sflux/sflux_inputs.txt","*.ll"]
		for part in filesToSearch:
			for files in glob.glob("%s/%s" % (runDir, part)):
				filesToMove.append(files)

		# Move files to input dir 
		for files in filesToMove:
			try:
				shutil.move(files, inputDir)
			except:
				print "Unable to move %s to %s as part of post run cleanup" % (each,path)

		# Create archive of inputs 
		arch = tarfile.open("%s/inputs.tar.gz" % runDir, "w:gz")
		print "Creating archive of input files, this may take a few minutes..."
		arch.add(inputDir)
		arch.close()

		# Remove duplicative directory - all inputs save in inputs.tar.gz file 
		shutil.rmtree(inputDir)

	def makeParaviewFile(self, path):
		""" Makes a call to create input file for Paraview """
		# TODO: Import these classes and make call and dump in post dir

	def getRunTime(self, path):
		""" Grab run time for run and stuff into model database """
		# TODO: Figure out how the heck to do this	
	
def cleanUpRun(path):
	cleanup = Clean(path)
	cleanup.filesToLog(cleanup.path)
	cleanup.makeRunArchive(cleanup.path)

#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
if __name__ == '__main__':
# Create a parser for command line args and options
	usage = ("Usage: %prog path")
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

# Grab path for run (directory above run)
	if len(args) != 1:
		parser.error("Must provide command line argument to path above run directory")
	else:
		path = args[0]
	
	cleanUpRun(path)
