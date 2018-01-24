# Code written by M. C. January 2014
# Code adapted by W. G. January 2018
#!/usr/bin/env python


from __future__ import print_function
import MDAnalysis as mda
from MDAnalysis.analysis.align import *
import numpy as np
import matplotlib
import matplotlib.pyplot, sys,numpy
import scipy.ndimage as ndi
import itertools
import math
import numpy.linalg as la
from itertools import cycle
import time
import os
import operator
import scipy.stats as stats

# --------------------------------------------
#                  main function
# --------------------------------------------

def main(coord, trajs, proteins_nb, index_prot, bb_prot, res_prot, start_frame, end_frame, delta_t, lip_sel, img_name):

	mda.core.flags['use_pbc'] = False

	# Load trajectory:
	if len(trajs) == 1:
		U = mda.Universe(coord, trajs[0])
	else:
		print('no files found load trajectory from  ')
		U = load_md_trajectory()
	print('Trajectory loaded.')


	# Make protein selections

	proteins = []
	proteins_name = []
	proteins_sele = {}
	end_prot = proteins_nb*index_prot

	# access to the selection using MDAnalysis:
	protein = U.select_atoms('protein')

	# store index (as a string) of the different proteins
	for p_index in range(proteins_nb):
		proteins_sele[p_index] = protein.select_atoms('bynum ' + str(1 + p_index*index_prot) + ':' + str((1 + p_index)*index_prot))
	print('Proteins found.')


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#
#						rotation calculation
#
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


	#define the delta t (in nb of frames) for vector calculation:
	dt = delta_t
	prot_number = 0

	#cutoff (in Angstrom) to display the density around the protein
	cutoff = 120
	cutoff1 = 120

	#data structures
	angle_list = np.zeros((proteins_nb, len(range(start_frame, end_frame, dt))))
	vector_list = np.zeros((proteins_nb, len(range(start_frame, end_frame, dt)), 1, 3))
	prot_list = np.zeros((proteins_nb, len(range(0, end_frame-start_frame, dt)), bb_prot, 3))
	frame_list = np.linspace(0, end_frame-start_frame, ((end_frame-start_frame)/dt)+1)

	neighbours_prot = {}
	neighbours_lip = {}

	for p_index in range(proteins_nb):
		neighbours_prot[p_index] = {}

	for p_index in range(proteins_nb):
		neighbours_lip[p_index] = {}

	#lipid definition
	lipids_sele_def = lip_sel  # like "name PO4"
	lipids_sele = U.select_atoms(lipids_sele_def)


	#browse frames
	#=============
	frame_num = 0
	U.trajectory.rewind()
	for frame in range(start_frame,end_frame,dt):

			#debug
			print(frame)

			#get frame properties
			ts = U.trajectory[frame]
			box_size = ts.dimensions[0:3]
			frame_in_title = str(frame_num).zfill(4)
			frame_legend1 = (frame)*2
			frame_legend2 = (frame+dt)*2

			#define container for proteins coordinates
			coord_t = {}
			coord_CoG = {}
			coord_pt = {}
			coord_dt = {}
			coord_CoG_dt = {}
			coord_pt_dt = {}

			#process proteins: frame t
			#------------------

			#get coords of all proteins (ensuring they're all within 0 and box_dim)
			coords_prot_removed_whole = np.zeros((proteins_nb, bb_prot, 3))

			#get coords of lipids (ensuring they're all within 0 and box_dim)
#			coords_lip = coords_remove_whole(lipids_sele.coordinates(), box_size)


			for p_index in range(proteins_nb):

				coords_prot_removed_whole[p_index, :] = coords_remove_whole(proteins_sele[p_index].BB.coordinates(), box_size)

			for p_index in range(proteins_nb):

				# time t ---------------------------------------------------------------------
				#store CoG of current prot
				coord_CoG[p_index] = calculate_cog(coords_prot_removed_whole[p_index, :], box_size)

				#shift all coords to -box_dim/2 +box_dim/2 referential
				tmp_prot_centered = coords_center_in_box(coords_prot_removed_whole, coord_CoG[p_index], box_size)

				#store neighbouring lipids coords in the same referential
#				tmp_lipid_centered = coords_center_in_box(coords_lip, coord_CoG[p_index], box_size)

				#store centered coords of current prot
				coord_t[p_index] = tmp_prot_centered[p_index]

				#store neighbouring proteins coords in the same referential
				not_p_index = list(range(proteins_nb))
				not_p_index.remove(p_index)
				tmp_prot_centered = tmp_prot_centered[not_p_index, :]
				tmp_prot_neighbour = np.ndarray.flatten(tmp_prot_centered)
				tmp_prot_neighbour.shape = ((proteins_nb-1)*bb_prot, 3)


				#coordinates of selected residue at t
				coord_pt = [coord_t[p_index][res_prot][0], coord_t[p_index][res_prot][1], 0]
				rot = rotation_matrix(angle_2D_between(coord_pt, [0, 1, 0]))
				vector_list[p_index, frame_num, :] = rotate_coord(coord_pt, rot)
				#print ('coord_pt, angle_2D_between(coord_pt, [0,1,0]), rot, vector_list[p_index,frame_num,:]', coord_pt, angle_2D_between(coord_pt, [0,1,0]), rot, vector_list[p_index,frame_num,:])

				# TO CHECK THE VALUES IN HARD
				neighbours_prot[p_index][frame_num] = rotate_coord(tmp_prot_neighbour, rot)
#				neighbours_prot[p_index][frame_num] = rotate_coord(tmp_lipid_centered,rot)

			# go through each protein and store reorientated coords, this is really for the plotting stage...

			for p1_index in range(proteins_nb):

				#retrieve selection
				prot = proteins_sele[p1_index]

				#put coord of proteins back in the box
				coord_dt = coords_remove_whole(prot.BB.coordinates(), box_size)

				#calculate CoG
				CoG_dt = calculate_cog(coord_dt, box_size)

				#center coord in -box_dim/2 +box_dim/2 referential
				coord_dt = coords_center_in_box(coord_dt, CoG_dt, box_size)

				# coordinates of L440 at t+dt
				coord_pt_dt = [coord_dt[440][0], coord_dt[440][1], 0]

				# print("")
				# print(coord_pt_dt)

				# angle_2D_between already returns the angle in RADIANS, don't need to convert to radians again!

				#rot = rotation_matrix(math.radians(angle_2D_between(coord_pt_dt, [0, 1, 0])))
				rot = rotation_matrix(angle_2D_between(coord_pt_dt, [0, 1, 0]))
				coord_pt_dt = rotate_coord(coord_pt_dt, rot)

				# print (coord_pt_dt)
				# print("")



				#store protein coordinates
				prot_list[p1_index, frame_num, :] = rotate_coord(coord_dt, rot)

			#update frame counter
			#--------------------
			frame_num += 1


	vector_list = vector_list.tolist()
	prot_list = prot_list.tolist()

	bins = np.arange(-cutoff, cutoff)
	array_bin = np.zeros((40, 40))
	xedges = []
	yedges = []
	list_angle = []
	#print(vector_list[5:10])
	for p_index in range(proteins_nb):
		#print('Got this far...')

		for frame in range(0,len(neighbours_prot[p_index])):
			#print('How about this...')
			#fig_vector_prot = matplotlib.pyplot.figure()
			#ax_vector_prot = fig_vector_prot.add_subplot(111)

			array_op, xedges, yedges = np.histogram2d(neighbours_prot[p_index][frame][:, 0], neighbours_prot[p_index][frame][:, 1], bins=40, range=[[-cutoff, cutoff], [-cutoff, cutoff]])
			array_bin = array_bin + array_op

			for vector in range(0, len(neighbours_prot[p_index][frame][:, 0])):

					dist = np.linalg.norm([neighbours_prot[p_index][frame][vector, 0], neighbours_prot[p_index][frame][vector, 1]])

					if dist < cutoff1:
						list_angle.append(math.degrees(angle_2D_between([neighbours_prot[p_index][frame][vector, 0], neighbours_prot[p_index][frame][vector, 1], 0], [0, 1, 0])))

	# create angle diagram
	fig_angle_prot = matplotlib.pyplot.figure()
	ax_angle_prot = fig_angle_prot.add_subplot(111)
	angle_distrib(list_angle, ax_angle_prot)
	title_angle_prot = 'plot_angle_prot_final.svg'
	fig_angle_prot.savefig(title_angle_prot, dpi=200)



	fig_vector_prot = matplotlib.pyplot.figure()
	ax_vector_prot = fig_vector_prot.add_subplot(111)
	plot_prot(prot_list[5][0], -cutoff, cutoff, -cutoff, cutoff, ax_vector_prot)
	plot_vector(vector_list[5][0], -cutoff, cutoff, -cutoff, cutoff, 'black', ax_vector_prot)

	plot_density_array(array_bin, xedges, yedges, -cutoff, cutoff, -cutoff, cutoff, 50, ax_vector_prot, 'BuPu')

	title_vector_prot = img_name + '.svg'
	fig_vector_prot.savefig(title_vector_prot, dpi=200)

# -----------------------------------------------------------------------------------------------
#
#
#                  Functions to be called
#
# -----------------------------------------------------------------------------------------------


def angle_2D_between(v1, v2):

	angle = 0
	v1_u = [v1[0], v1[1]] / la.norm([v1[0], v1[1]])
	v2_u = [v2[0], v2[1]] / la.norm([v2[0], v2[1]])
	dot =  np.dot(v1_u, v2_u)
	#print(dot)
	if (dot >= 1):
		angle = 0
	elif (dot <= -1):
		angle = 1
	else:
		angle = math.acos(dot) # acos returns angle in RADIANS

	sign = v1[0]*v2[1] - v2[0]*v1[1]

	if sign < 0:
		angle = -angle

	return angle


def abs_angle_2D_between(v1, v2):

	angle = 0
	v1_u = [v1[0],v1[1]] / la.norm([v1[0],v1[1]])
	v2_u = [v2[0],v2[1]] / la.norm([v2[0],v2[1]])
	dot =  np.dot(v1_u, v2_u)
	#print(dot)
	if (dot >= 1):
		angle = 0
	elif (dot <= -1):
		angle = 180
	else:
		angle = math.degrees(math.acos(dot))

	return angle


def calculate_cog(tmp_coords, box_dim):

	#this method allows to take pbc into account when calculcating the center of geometry
	#see: http://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions

	cog_coord = np.zeros(3)
	tmp_nb_atoms = np.shape(tmp_coords)[0]

	for n in range(0,3):
		tet = tmp_coords[:,n] * 2 * math.pi / float(box_dim[n])
		xsi = np.cos(tet)
		zet = np.sin(tet)
		tet_avg = math.atan2(-np.average(zet),-np.average(xsi)) + math.pi
		cog_coord[n] = tet_avg * box_dim[n] / float(2*math.pi)

	return cog_coord

def rotation_matrix(angle):
	# angle should be in RADIANS
	matrix = np.array([[math.cos(angle), -math.sin(angle), 0], [math.sin(angle),  math.cos(angle), 0], [0, 0, 1]])

	return matrix

def load_md_trajectory():
	# Load trajectory
	print("Loading trajectory.")
	root_dir = '/sansom/n24/chavent/Big_OmpS/model/MARTINI/6x6/'
	#universe = mda.Universe( 'start_mix6_6_M.skip10_now.gro',
    #                ['md_prod_mix6_6_M.part1-3.skip10_now.xtc',
    #                 'md_prod_mix6_6_M.part0004.skip10_now.xtc',
    #                 'md_prod_mix6_6_M.part0005.skip10_now.xtc',
    #                 'md_prod_mix6_6_M.part0006.skip10_now.xtc',
    #                 'md_prod_mix6_6_M.part0007.skip10_now.xtc',
    #                 'md_prod_mix6_6_M.part0008.skip10_now.xtc',
    #                 'md_prod_mix6_6_M.part0009.skip10_now.xtc',
    #                 'md_prod_mix6_6_M.part0010.skip10_now.xtc',
    #                 'md_prod_mix6_6_M.part0011.skip10_now.xtc'])
	universe = mda.Universe( 'start_mix6_6_M.skip10_now.gro',
	                [ 'md_prod_mix6_6_M.part0011.skip10_now.xtc'])



	print("Finished loading.")
	return universe

def coords_remove_whole(coords, box_dim):
	#this function ensures the coordinates are within 0 and box_dim

	coords[:,0] -= np.floor(coords[:,0]/float(box_dim[0])) * box_dim[0]
	coords[:,1] -= np.floor(coords[:,1]/float(box_dim[1])) * box_dim[1]

	return coords

# get the coordinates of the proteins
def coords_center_in_box(coords, cog, box_dim):

	coords_loc = np.copy(coords)

	if len(np.shape(coords_loc)) > 2:
		#center lipids coordinates on cluster (x,y) coordinates
		coords_loc[:,:,0] -= cog[0]
		coords_loc[:,:,1] -= cog[1]

		#deal with pbc
		coords_loc[:,:,0] -= (np.floor(2*coords_loc[:,:,0]/float(box_dim[0])) + (1-np.sign(coords_loc[:,:,0]))/float(2)) * box_dim[0]
		coords_loc[:,:,1] -= (np.floor(2*coords_loc[:,:,1]/float(box_dim[1])) + (1-np.sign(coords_loc[:,:,1]))/float(2)) * box_dim[1]
	else:
		#center lipids coordinates on cluster (x,y) coordinates
		coords_loc[:,0] -= cog[0]
		coords_loc[:,1] -= cog[1]

		#deal with pbc
		coords_loc[:,0] -= (np.floor(2*coords_loc[:,0]/float(box_dim[0])) + (1-np.sign(coords_loc[:,0]))/float(2)) * box_dim[0]
		coords_loc[:,1] -= (np.floor(2*coords_loc[:,1]/float(box_dim[1])) + (1-np.sign(coords_loc[:,1]))/float(2)) * box_dim[1]

	return coords_loc

def rotate_coord(coords, rotation_matrix):

	coords_loc = (np.copy(coords)).transpose()
	#print (np.shape(coords_loc))
	coords_loc = rotation_matrix.dot(coords_loc)
	coords_loc = coords_loc.transpose()
	#print (np.shape(coords_loc))
	return coords_loc

def calcul_std_lower_bound(list_std, list_avg):

	length_list = len(list_avg)
	std_lower_list = []

	for i in range(0, length_list):

		lower = list_avg[i] - list_std[i]
		std_lower_list.append(lower)

	return std_lower_list

def calcul_std_upper_bound(list_std, list_avg):

	length_list = len(list_avg)
	std_upper_list = []

	for i in range(0, length_list):

		lower = list_avg[i] + list_std[i]
		std_upper_list.append(lower)

	return std_upper_list

# -------------------------------------------------------------
#
#                  plotting functions
#
# -------------------------------------------------------------

def angle_distrib (list_values, ax):

	#print(list_values)
	# ompf
	#n, bins, patches = matplotlib.pyplot.hist(list_values, 50, histtype='stepfilled',normed=1, facecolor='orange', alpha=0.5)
	# btub
	n, bins, patches = matplotlib.pyplot.hist(list_values, 50, histtype='stepfilled', facecolor='gray', alpha=0.5)

	ax.set_xlim([-180,180])
	ax.set_xticks(numpy.arange(-180,180,20))
	#ax.set_ylim([ymin,ymax])

	#matplotlib.pyplot.savefig(str(name)+'_'+str(dt)+'.svg')
	#matplotlib.pyplot.close()

# plot the coordinates of the proteins and the angle in the corner
def plot_prot_angle(coord, xmin, xmax, ymin, ymax, angle, frame, ax):

	time = frame*2

	# label matplotlib axes
	matplotlib.pyplot.xlabel('X [A]', size=15)
	matplotlib.pyplot.ylabel('Y [A]', size=15)

	ticks_x = int((xmax-xmin)/10)
	ticks_y = int((ymax-ymin)/10)


	pos_leg_angle_x = int(xmin+(xmax-xmin)*0.7)
	pos_leg_angle_y = int(ymin+(ymax-ymin)*0.9)

	angle = round(angle, 2)
	legend_angle = "angle = "+str(angle)
	ax.text( pos_leg_angle_x, pos_leg_angle_y, legend_angle)


	pos_leg_frame_x = int(xmin+(xmax-xmin)*0.7)
	pos_leg_frame_y = int(ymin+(ymax-ymin)*0.95)

	legend_frame = "time = "+str(frame)+" ns"
	ax.text( pos_leg_frame_x, pos_leg_frame_y, legend_frame)



	# definition of axes for the graph
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.set_xticks(numpy.arange(xmin,xmax,ticks_x))
	ax.set_yticks(numpy.arange(ymin,ymax,ticks_y))
	ax.set_aspect('equal')

	#BtuB
	ax.plot(map(lambda x: x[0], coord),map(lambda x: x[1], coord),'o',markerfacecolor='#88C578',markeredgecolor='#000000',markersize=9,alpha=0.8)
	#OmpF
	#ax.plot(coord[:,0],coord[:,1],'o',markerfacecolor='#B29007',markeredgecolor='#000000',markersize=10,alpha=0.8)

	#print("coordinates plotted")


# plot the coordinates of the proteins
def plot_prot(coord, xmin, xmax, ymin, ymax, ax):

	from matplotlib.ticker import FuncFormatter, MaxNLocator

	# label matplotlib axes
	matplotlib.pyplot.xlabel('X [A]', size=15)
	matplotlib.pyplot.ylabel('Y [A]', size=15)

	ticks_x = int((xmax-xmin)/100)
	ticks_y = int((ymax-ymin)/100)

	# definition of axes for the graph
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])
	ax.set_xticks(numpy.arange(xmin, xmax, ticks_x))
	ax.set_yticks(numpy.arange(ymin, ymax, ticks_y))
	ax.set_aspect('equal')

	matplotlib.pyplot.xlabel('X / nm', size=15)
	matplotlib.pyplot.ylabel('Y / nm', size=15)

	import matplotlib.ticker as ticker

	scale = 0.1
	ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x * scale))
	ax.xaxis.set_major_formatter(ticks)
	ax.yaxis.set_major_formatter(ticks)
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))

	#BtuB
	#ax.plot(map(lambda x: x[0], coord),map(lambda x: x[1], coord),'o',markerfacecolor='#88C578',markeredgecolor='#000000',markersize=7,alpha=0.8)
	#OmpF

	regions_of_interest = [coord[89: 193], coord[395:509], coord[649: 759], coord[949: 1074]]
	col = ['#ffffcc', '#a1dab4', '#41b6c4', '#225ea8']

	# Specifiy the non-VSD domains

	not_vsd = coord

	not_vsd_elements = list(range(89, 193)) + list(range(395, 509)) + list(range(629, 759)) + list(range(949, 1074))

	# filter out the VSD's, need to reverse order so that subsequent indexes do not drop off

	for index in sorted(not_vsd_elements, reverse=True):
		del not_vsd[index]

	# ax.plot(map(lambda x: x[0], coord), map(lambda x: x[1], coord),'o', markerfacecolor='#B29007', markeredgecolor='#000000', markersize=10, alpha=0.8)

	ax.plot(map(lambda x: x[0], not_vsd), map(lambda x: x[1], not_vsd), '.', markerfacecolor='snow',
			markersize=10, markeredgecolor='#000000', alpha=0.6)

	for i in range(4):
		j = i + 1

		ax.plot(map(lambda x: x[0], regions_of_interest[i]), map(lambda x: x[1], regions_of_interest[i]), '.',
				markerfacecolor=col[i],
				markeredgecolor='#000000', markersize=10, alpha=0.8, label='VSD D' + str(j))

	ax.legend()

	#ax.plot(map(lambda x: x[0], coord),map(lambda x: x[1], coord),'o',markerfacecolor='#B29007',markeredgecolor='#000000',markersize=10,alpha=0.8)

	#print("coordinates plotted")

def plot_vector(vector, xmin, xmax, ymin, ymax, color_arrow, ax):


	# label matplotlib axes
	matplotlib.pyplot.xlabel('X [A]', size=15)
	matplotlib.pyplot.ylabel('Y [A]', size=15)

	ticks_x = int((xmax-xmin)/10)
	ticks_y = int((ymax-ymin)/10)

	# definition of axes for the graph
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])
	ax.set_xticks(numpy.arange(xmin, xmax, ticks_x))
	ax.set_yticks(numpy.arange(ymin, ymax, ticks_y))
	ax.set_aspect('equal')
	ax.annotate("", xy=(vector[0][0], vector[0][1]), xytext=(0,0),arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color=color_arrow))



# plot the density around a protein and plot the result
def plot_density(list_x, list_y, xmin, xmax, ymin, ymax, bin_num, ax, color):


	# label matplotlib axes
	matplotlib.pyplot.xlabel('X [A]', size=15)
	matplotlib.pyplot.ylabel('Y [A]', size=15)
	ticks_x = int((xmax-xmin)/10)
	ticks_y = int((ymax-ymin)/10)

	# definition of axes for the graph
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])
	ax.set_xticks(numpy.arange(xmin, xmax, ticks_x))
	ax.set_yticks(numpy.arange(ymin, ymax, ticks_y))
	ax.set_aspect('equal')


	array_op,xedges,yedges = np.histogram2d(list_x,list_y, bins=bin_num,range=[[xmin,xmax],[ymin,ymax]])

	max_array = max(array_op.flatten())
	array_avg = array_op/max_array
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
	# interesting interpolation: bicubic, lanczos see http://matplotlib.org/examples/images_contours_and_fields/interpolation_methods.html
	cax = ax.imshow(array_avg.T,extent=extent, interpolation='bicubic', origin='lower', cmap=color)
	matplotlib.pyplot.colorbar(cax)
	#ax.plot(list_x,list_y,'o',markerfacecolor='#FC3147',markeredgecolor='#FC3147',markersize=0.8,alpha=0.8)

# plot the density around a protein and plot the result
def plot_density_array(array,xedges,yedges, xmin, xmax, ymin, ymax, bin_num, ax, color):


	# label matplotlib axes
	matplotlib.pyplot.xlabel('X [A]', size=15)
	matplotlib.pyplot.ylabel('Y [A]', size=15)
	#ticks_x = int((xmax-xmin)/10)
	#ticks_y = int((ymax-ymin)/10)

	ticks_x = 11
	ticks_y = 11

	# print (ymax)
	# print (ymin)
	# print (xmax)
	# print (xmin)
	# print (ticks_x)
	# print (ticks_y)

	# definition of axes for the graph
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])
	ax.set_xticks(numpy.linspace(xmin, xmax, ticks_x))
	ax.set_yticks(numpy.linspace(ymin, ymax, ticks_y))

	# ax.set_xticks(numpy.arange(xmin, xmax, ticks_x))
	# ax.set_yticks(numpy.arange(ymin, ymax, ticks_y))



	ax.set_aspect('equal')


	max_array=max(array.flatten())
	array_avg = array/max_array
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
	# interesting interpolation: bicubic, lanczos see http://matplotlib.org/examples/images_contours_and_fields/interpolation_methods.html
	cax= ax.imshow(array_avg.T,extent=extent,interpolation='bicubic',origin='lower',cmap=color)
	matplotlib.pyplot.colorbar(cax)
	#ax.plot(list_x,list_y,'o',markerfacecolor='#FC3147',markeredgecolor='#FC3147',markersize=0.8,alpha=0.8)


# plot the disp around a protein and plot the result
def plot_disp(list_x, list_y, list_disp, xmin, xmax, ymin, ymax, bin_num, ax, color):


	# label matplotlib axes
	matplotlib.pyplot.xlabel('X [A]', size=15)
	matplotlib.pyplot.ylabel('Y [A]', size=15)
	ticks_x = int((xmax-xmin)/10)
	ticks_y = int((ymax-ymin)/10)

	# definition of axes for the graph
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.set_xticks(numpy.arange(xmin,xmax,ticks_x))
	ax.set_yticks(numpy.arange(ymin,ymax,ticks_y))
	ax.set_aspect('equal')


	array_op,xedges,yedges = np.histogram2d(list_x,list_y, weights = list_disp, bins=bin_num,range=[[xmin,xmax],[ymin,ymax]])

	max_array=max(array_op.flatten())
	array_avg = array_op/max_array
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
	cax= ax.imshow(array_avg, extent=extent,interpolation='lanczos',origin='lower',cmap=color)
	matplotlib.pyplot.colorbar(cax)
	#ax.plot(list_x,list_y,'o',markerfacecolor='#FC3147',markeredgecolor='#FC3147',markersize=0.8,alpha=0.8)


# plot the curve for the rotation of a protein in function of the time
def plot_angle_curve(angle_array, frame_array, ax):


	# label matplotlib axes
	matplotlib.pyplot.xlabel('frames', size=15)
	matplotlib.pyplot.ylabel('Angle (deg.)', size=15)


	angle_list = angle_array.tolist()
	frame_list = frame_array.tolist()

	min_frame = min(frame_list)
	max_frame = max(frame_list)

	min_angle = np.amin(angle_array) - 20
	max_angle = np.amax(angle_array) + 20

	# definition of axes for the graph
	ax.set_xlim([min_frame, max_frame])
	ax.set_ylim([min_angle,max_angle])

	for angle in angle_list:

		ax.plot(frame_list, angle, color=numpy.random.rand(3,1))


# plot the curve for the rotation of a protein in function of the time
def plot_angle_curve_single(angle_array, frame_array, max_angle, ax):


	# label matplotlib axes
	matplotlib.pyplot.xlabel('frames', size=15)
	matplotlib.pyplot.ylabel('Angle (deg.)', size=15)


	min_frame = min(frame_array)
	max_frame = max(frame_array)

	min_angle = 0

	# definition of axes for the graph
	ax.set_xlim([min_frame, max_frame])
	ax.set_ylim([min_angle,max_angle])

	ax.plot(frame_array, angle_array, color=numpy.random.rand(3,1))



def plot_angle_avg_std(avg_array, std_array, frame_array, ax):

	# label matplotlib axes
	matplotlib.pyplot.xlabel('frames', size=15)
	matplotlib.pyplot.ylabel('Angle (deg.)', size=15)

	avg_list = avg_array.tolist()
	frame_list = frame_array.tolist()
	std_list = std_array.tolist()

	std_lower_list = calcul_std_lower_bound(std_list, avg_list)
	std_upper_list = calcul_std_upper_bound(std_list, avg_list)

	min_frame = min(frame_list)
	max_frame = max(frame_list)

	min_std = min(std_lower_list)
	max_std = max(std_upper_list)

	min_angle = min(avg_list) - 20 - min_std
	max_angle = max(avg_list) + 20 + max_std

	# definition of axes for the graph
	ax.set_xlim([min_frame, max_frame])
	ax.set_ylim([min_angle,max_angle])

	ax.plot(frame_list, avg_list, color='black', ls='--')
	ax.fill_between(frame_list, std_lower_list, std_upper_list, facecolor='yellow', alpha=0.5)


# ----------------------------------------
#                  Parser
# ----------------------------------------


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--traj', metavar='TRAJ_FILE', type=str, help='Trajectory files', nargs='+')
    parser.add_argument('--coord', type=str, help='System coordinate file', required=True,
                        dest='coord')
    parser.add_argument('--num-prot', type=int, help='Number of proteins', default=1,
                        dest='proteins_nb')
    parser.add_argument('--index_prot', type=int, help='nb. atoms', default=1,
                        dest='index_prot')
    parser.add_argument('--bb_prot', type=int, help='nb. backbone atoms', default=1,
                        dest='bb_prot')
    parser.add_argument('--res_prot', type=int, help='Residue nb. of protein to align the proteins', default=1,
                        dest='res_prot')
    parser.add_argument('--sf', type=int, help='first frame number to use', default=1,
                        dest='start_frame')
    parser.add_argument('--ef', type=int, help='last frame number to use', default=1,
                        dest='end_frame')
    parser.add_argument('--dt', type=int, help='delta t', default=1,
                        dest='delta_t')
    parser.add_argument('--sel', type=str, help='lipid selection', default=1,
                        dest='lip_sel')
    parser.add_argument('--name', type=str, help='name of the image output', default=1,
                        dest='img_name')
    args = parser.parse_args()


    main(args.coord, args.traj, args.proteins_nb, args.index_prot, args.bb_prot, args.res_prot, args.start_frame, args.end_frame, args.delta_t, args.lip_sel, args.img_name)
