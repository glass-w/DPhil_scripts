'''Script to try and center 'no jump-ed' trajectory so that it doesn't jitter, as when using center with trjconv'''

import numpy as np
import MDAnalysis
import argparse
import progressbar
import time

def centre(gro_file, nojump_xtc_file, output_xtc, centre_type = 'CoM', skip=0):
	progressBar = progressbar.AnimatedProgressBar(end=100, width=75)
	progressBar.show_progress()

	u = MDAnalysis.Universe(gro_file, nojump_xtc_file)

	progressBar+2
	progressBar.show_progress()

	allAtoms = u.selectAtoms("all")

	# open a stream to write the frames of the XTC file to
	numAtoms = allAtoms.numberOfAtoms()
	FileWriter=MDAnalysis.coordinates.core.writer(filename=output_xtc,numatoms=numAtoms)

	progressBar+3
	progressBar.show_progress()
	progressBar.show_progress()
	progressPerFrame = 95.0/(u.trajectory.numframes)

	for ts in u.trajectory[::skip]:
	
		progressBar+progressPerFrame
		progressBar.show_progress()

		if centre_type == 'CoM':
			allAtoms.translate(-allAtoms.centerOfMass())
		elif centre_type == 'CoG':
			allAtoms.translate(-allAtoms.centerOfGeometry())
		FileWriter.write(allAtoms)
		
	FileWriter.close()
	print 'File writing finished'
	
	return
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Centres no jump trajectory according to Centre of Mass (CoM) of Centre of Geometry (CoG) to prevent drift")
	parser.add_argument("-s", dest="gro_file", type=str, help='the coordinate file [.gro]')
	parser.add_argument("-f", dest="nojump_xtc_file",type=str, help='the coordinates to fix -  a "no jump" xtc file')
	parser.add_argument("-o", dest="output_xtc", type=str, help='the name of the output file. Specified by the file extension (.gro, .pdb, .xtc, etc).')
	parser.add_argument("-c", dest="centre_type", type=str, help='Choice of centring methods - either "CoM" for centre of mass, or "CoG" for centre of geometry - default: "CoM"', default='CoM', choices = ['CoM','CoG'])
	parser.add_argument("-skip", dest="skip", type=int, help='Number of frames to skip from input - default: 1 (ie. no skipping of frames)', default=1)
	options = parser.parse_args()
	# create a simple progress bar that will be updated on the command line - thanks to Phil F
	startTime = time.time()
	centre(options.gro_file, options.nojump_xtc_file, options.output_xtc, options.centre_type, options.skip)
	endTime = time.time()
	print 'Program took {:7.1f} seconds to run'.format(endTime-startTime)
