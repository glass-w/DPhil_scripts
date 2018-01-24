from __future__ import print_function
import MDAnalysis as mda
from MDAnalysis.analysis.align import *
import numpy as np
from numpy import random
import matplotlib
import matplotlib.pyplot, sys,numpy
import scipy.ndimage as ndi
import itertools
import math
import numpy.linalg as la
from itertools import cycle
import time
import os
from scipy import ndimage

def main(coord, trajs, proteins1_nb, proteins2_nb, index_prot1, index_prot2):

	MDAnalysis.core.flags['use_pbc'] = False

	# Load trajectory:
	if len(trajs) == 1:
		U = mda.Universe(coord, trajs[0]) 
	else:
		print('no files found load trajectory from  ')
		U = load_md_trajectory()
	print('Trajectory loaded.')


	# Make protein selections
	# for info: BtuB index (prot1) = 1332, OmpF index (prot2) = 2208

	proteins1 = []
	proteins2 = []
	proteins1_name = []
	proteins2_name = []
	proteins1_sele = {}
	proteins2_sele = {}
	end_prot1 = proteins1_nb*index_prot1
	end_prot2 = proteins2_nb*index_prot2+end_prot1
	
	# access to the selection using MDAnalysis:
	protein = U.selectAtoms('protein')
	
	# store index (as a string) of the different proteins BtuB
	for p1_index in range(proteins1_nb):
		proteins1_sele[p1_index] = protein.selectAtoms('bynum ' + str(1 + p1_index*index_prot1) + ':' + str((1 + p1_index)*index_prot1))
	print('Proteins BtuB found.')

	# store index (as a string) of the different proteins ompF
	for p2_index in range(proteins2_nb):
		proteins2_sele[p2_index] = protein.selectAtoms('bynum ' + str(1 + p2_index*index_prot2 + end_prot1) + ':' + str((1 + p2_index)*index_prot2 + end_prot1))
	print('Proteins OmpF found.')



#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#
#						store data in a list and then plot the results
#
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


	#define the delta t (in nb of frames) for vector calculation:
	dt= 5
	prot_number = 0
	last_frame = 2000 - 2*dt
	

	#data structures
	#---------------

	CoG_list = np.zeros((proteins1_nb, len(range(0,last_frame,dt)),1 ,3))		

		
	#browse frames
	#=============
	nb_frame = 0
	U.trajectory.rewind()
	for frame in range(0,last_frame,dt):
			
			#debug
			print(frame)

			#get frame properties
			ts = U.trajectory[frame]

			#define container for CoG coordinates
			coord_CoG = {}					
						
						
			for p1_index in range(proteins1_nb):				
						
				#store CoG for each frame
				CoG_list[p1_index, nb_frame] = proteins1_sele[p1_index].BB.centerOfMass()

			#update frame counter
			#--------------------
			nb_frame += 1
				
	CoG_list2 = CoG_list.tolist()

	for i in range(0,len(CoG_list2[40])):
		 

		CoG = CoG_list2[40][i]
		print(i)

		counter = str(i).zfill(4)
		outputname =  'plot'+'_prot40_'+counter+'.tiff'

				 		
		make_image(CoG[0][0], CoG[0][1], -200, 1200, -200, 1200, 400, outputname) 		
		

# -----------------------------------------------------------------------------------------------
#
#
#                  Functions to be called 
#
# -----------------------------------------------------------------------------------------------



def make_image(x, y, xmin, xmax, ymin, ymax, bin_num, outputname):

	fig = matplotlib.pyplot.figure()
	fig.set_size_inches(1, 1)
	ax = matplotlib.pyplot.Axes(fig, [0, 0 , 1 ,1])
	#ax =fig.add_subplot(111)
	ax.set_xlim([xmin,xmax])		
	ax.set_ylim([ymin,ymax])
	
	ax.set_axis_off()
	fig.add_axes(ax) 

	list_x = [x]
	list_y = [y]
		
	#for i in range (-resolution, resolution):
	#	for j in range (-resolution, resolution):	
	#		list_x.append(x + i)
	#		list_y.append(y + j)			
	
			
	array_op,xedges,yedges = np.histogram2d(list_x,list_y, bins=bin_num,range=[[xmin,xmax],[ymin,ymax]])
	array_op = ndimage.gaussian_filter(array_op.astype(float),4.1)

	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
	# interesting interpolation: bicubic, lanczos see http://matplotlib.org/examples/images_contours_and_fields/interpolation_methods.html
	ax.imshow(array_op.T,extent=extent,interpolation='nearest',origin='lower',cmap='bone', aspect = 'normal')
	matplotlib.pyplot.savefig(outputname, dpi = 512)



# ----------------------------------------
#                  Parser 
# ----------------------------------------


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--traj', metavar='TRAJ_FILE', type=str, help='Trajectory files', nargs='+')
    parser.add_argument('--coord', type=str, help='System coordinate file', required=True,
                        dest='coord')
    parser.add_argument('--num-prot1', type=int, help='Number of proteins 1', default=1,
                        dest='proteins1_nb')
    parser.add_argument('--num-prot2', type=int, help='Number of proteins 2', default=1,
                        dest='proteins2_nb')
    parser.add_argument('--index_prot1', type=int, help='nb. atoms of protein1', default=1,
                        dest='index_prot1')
    parser.add_argument('--index_prot2', type=int, help='nb. atoms of protein2', default=1,
                        dest='index_prot2')
    args = parser.parse_args()
    main(args.coord, args.traj, args.proteins1_nb, args.proteins2_nb, args.index_prot1, args.index_prot2)
