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
import itertools
import matplotlib.mlab as mlab
import scipy.stats as stats
from scipy.stats import norm,rayleigh, halfnorm
from scipy.stats import chisquare

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#
#						store displacement data in numpy array
#
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
def main(coord, trajs, proteins1_nb, proteins2_nb, index_prot1, index_prot2, clustfile):

	mda.core.flags['use_pbc'] = False

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
	# hack for OmpF
	#end_prot1 = 72*index_prot1	
	end_prot2 = proteins2_nb*index_prot2+end_prot1
	
	# access to the selection using MDAnalysis:
	protein = U.select_atoms('protein')
	
	# store index (as a string) of the different proteins BtuB

	if (proteins1_nb != 0):
		for p1_index in range(proteins1_nb):
			proteins1_sele[p1_index] = protein.select_atoms('bynum ' + str(1 + p1_index*index_prot1) + ':' + str((1 + p1_index)*index_prot1))
		print('Proteins BtuB found.')

	# store index (as a string) of the different proteins ompF
	if (proteins2_nb != 0):
		for p2_index in range(proteins2_nb):
			proteins2_sele[p2_index] = protein.select_atoms('bynum ' + str(1 + p2_index*index_prot2 + end_prot1) + ':' + str((1 + p2_index)*index_prot2 + end_prot1))
		print('Proteins OmpF found.')



	#define the delta t (in nb of frames) for vector calculation:
	dt=200 # 1 frame per ns (dt=1000)
	prot_number = 0
	start_frame = 0
	last_frame = U.trajectory.n_frames - 2*dt

	# load the file
	filename = clustfile


	# convert the text in numpy array
	a = np.loadtxt(filename)
        
	list_lines = np.arange(start_frame,last_frame,dt)

	line_BtuB_end = proteins1_nb +1
	line_OmpF_start = proteins1_nb +1
	line_OmpF_end = proteins1_nb + proteins2_nb + 2
	#matrix_clust = a[list_lines,1:line_BtuB_end] # BtuB mixed system
	#matrix_clust = a[list_lines,line_OmpF_start:line_OmpF_end] # OmpF mixed system	
	matrix_clust = a[list_lines,1:line_BtuB_end] # single mol. system

	shape_matrix_cluster = np.shape(matrix_clust)
	print ("cluster shape: "+str(np.shape(matrix_clust)))
	max_cluster = int(np.max(matrix_clust)+1) 



#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#
#						rotation calculation
#
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


def rot_trans_vel(dt):

	#data structures
	#---------------

	BB_number1 = 0
	BB_list_1 = []	

	BB_number2 = 0
	BB_list_2 = []	

	#data structures
	#---------------

	if (proteins1_nb!=0):
		BB_number1 = len(proteins1_sele[0].BB.coordinates())
		BB_list_1 = np.zeros((proteins1_nb, len(range(0,last_frame,dt)),BB_number1,3))			

	if (proteins2_nb!=0):
		BB_number2 = len(proteins2_sele[0].BB.coordinates())
		BB_list_2 = np.zeros((proteins2_nb, len(range(0,last_frame,dt)),BB_number2,3))		



	angle_list1 =[]
	angle_list2 =[]

	distance_list1 =[]
	distance_list2 =[]

	coord_list1_x =[]
	coord_list2_x =[]

	coord_list1_y =[]
	coord_list2_y =[]



	#browse frames
	#=============
	nb_frame = 0
	U.trajectory.rewind()
	for frame in range(start_frame,last_frame,dt):
	
			#debug
			print(frame)

			#get frame properties
			ts = U.trajectory[frame]
			box_size = ts.dimensions[0:3]


			#define container for proteins coordinates
			coord_CoG1 = {}			
			coord_CoG1_dt = {}
			coord_t1 = {}
			coord_pt1 = {}
			coord_pt1_dt = {}
				
			coord_CoG2 = {}			
			coord_CoG2_dt = {}
			coord_t2 = {}			
			coord_pt2 = {}
			coord_pt2_dt = {}

			#process proteins 1: frame t
			#------------------

			#get coords ofproteins 
			coords_prot1_removed_whole = np.zeros((proteins1_nb, BB_number1, 3))
			coords_prot2_removed_whole = np.zeros((proteins2_nb, BB_number2, 3))

			if (proteins1_nb!=0):
				for p1_index in range(proteins1_nb):
					coords_prot1_removed_whole[p1_index, :] = coords_remove_whole(proteins1_sele[p1_index].BB.coordinates(), box_size)


			if (proteins2_nb!=0):
				for p2_index in range(proteins2_nb):
					coords_prot2_removed_whole[p2_index, :] = coords_remove_whole(proteins2_sele[p2_index].BB.coordinates(), box_size)




			if (proteins1_nb!=0):
				for p1_index in range(proteins1_nb):				
				
					#store CoG of current prot
					coord_CoG1[p1_index] = calculate_cog(coords_prot1_removed_whole[p1_index, :], box_size)

					#shift all coords to -box_dim/2 +box_dim/2 referential
					tmp_prot1_centered = coords_center_in_box(coords_prot1_removed_whole, coord_CoG1[p1_index], box_size)

					#store centered coords of current prot
					coord_t1[p1_index] = tmp_prot1_centered[p1_index]
			
					#coordinates of W144 at t
					coord_pt1[p1_index] = [coord_t1[p1_index][144][0], coord_t1[p1_index][144][1], 0]				

			if (proteins2_nb!=0):
				for p2_index in range(proteins2_nb):
			
					#store CoG of current prot
					coord_CoG2[p2_index] = calculate_cog(coords_prot2_removed_whole[p2_index, :], box_size)

					#shift all coords to -box_dim/2 +box_dim/2 referential
					tmp_prot2_centered = coords_center_in_box(coords_prot2_removed_whole, coord_CoG2[p2_index], box_size)

					#store centered coords of current prot
					coord_t2[p2_index] = tmp_prot2_centered[p2_index]
			
					#coordinates of W144 at t
					coord_pt2[p2_index] = [coord_t2[p2_index][185][0], coord_t2[p2_index][185][1], 0]							


			#process proteins 1: frame t+dt
			#------------------		
			ts = U.trajectory[frame+dt]
			box_size = ts.dimensions[0:3]
	
			if (proteins1_nb!=0):			
				for p1_index in range(proteins1_nb):
			

					#retrieve selection
					prot = proteins1_sele[p1_index]
		
					#put coord of proteins back in the box
					coord_dt = coords_remove_whole(proteins1_sele[p1_index].BB.coordinates(), box_size)

					#calculate CoG
					CoG1_dt = calculate_cog(coord_dt, box_size)
		
					#center coord in -box_dim/2 +box_dim/2 referential
					coord_dt = coords_center_in_box(coord_dt, CoG1_dt, box_size)

					#store protein coordinates		
					BB_list_1[p1_index, nb_frame, :] = coord_dt


					# coordinates of W144 at t+dt
					coord_pt1_dt = [coord_dt[144][0],coord_dt[144][1],0]
				
					#calculate and store angle
					if nb_frame > 1:
						angle_list1.append(abs_angle_2D_between(coord_pt1[p1_index], coord_pt1_dt))
						distance_list1.append((mda.lib.distances.distance_array(np.float32(np.array([CoG1_dt[0],CoG1_dt[1],CoG1_dt[2]]).reshape((1,3))), np.float32(np.array([coord_CoG1[p1_index][0],coord_CoG1[p1_index][1],CoG1_dt[2]]).reshape((1,3))), ts.dimensions))[0][0])
						
						x_sign= np.float32(CoG1_dt[0]-coord_CoG1[p1_index][0]) 
						x_value = (mda.lib.distances.distance_array(np.float32(np.array([CoG1_dt[0],CoG1_dt[0],CoG1_dt[0]]).reshape((1,3))), np.float32(np.array([coord_CoG1[p1_index][0],CoG1_dt[0],CoG1_dt[0]]).reshape((1,3))), ts.dimensions))[0][0]
						if x_sign < 0.0 : x_value = -x_value

						y_sign= np.float32(CoG1_dt[1]-coord_CoG1[p1_index][1]) 
						y_value = (mda.lib.distances.distance_array(np.float32(np.array([CoG1_dt[1],CoG1_dt[1],CoG1_dt[1]]).reshape((1,3))), np.float32(np.array([coord_CoG1[p1_index][1],CoG1_dt[1],CoG1_dt[1]]).reshape((1,3))), ts.dimensions))[0][0]
						if y_sign < 0.0 : y_value = -y_value

						coord_list1_x.append(x_value)
						coord_list1_y.append(y_value)

					else:
						angle_list1.append( abs_angle_2D_between(coord_pt1[p1_index], coord_pt1_dt))
						distance_list1.append((mda.lib.distances.distance_array(np.float32(np.array([CoG1_dt[0],CoG1_dt[1],CoG1_dt[2]]).reshape((1,3))), np.float32(np.array([coord_CoG1[p1_index][0],coord_CoG1[p1_index][1],CoG1_dt[2]]).reshape((1,3))), ts.dimensions))[0][0])
						
						x_sign= np.float32(CoG1_dt[0]-coord_CoG1[p1_index][0]) 
						x_value = (mda.lib.distances.distance_array(np.float32(np.array([CoG1_dt[0],CoG1_dt[0],CoG1_dt[0]]).reshape((1,3))), np.float32(np.array([coord_CoG1[p1_index][0],CoG1_dt[0],CoG1_dt[0]]).reshape((1,3))), ts.dimensions))[0][0]
						if x_sign < 0.0 : x_value = -x_value

						y_sign= np.float32(CoG1_dt[1]-coord_CoG1[p1_index][1]) 
						y_value = (mda.lib.distances.distance_array(np.float32(np.array([CoG1_dt[1],CoG1_dt[1],CoG1_dt[1]]).reshape((1,3))), np.float32(np.array([coord_CoG1[p1_index][1],CoG1_dt[1],CoG1_dt[1]]).reshape((1,3))), ts.dimensions))[0][0]
						if y_sign < 0.0 : y_value = -y_value

						coord_list1_x.append(x_value)
						coord_list1_y.append(y_value)

		
			if (proteins2_nb!=0):			
				for p2_index in range(proteins2_nb):
			

					#retrieve selection
					prot = proteins2_sele[p2_index]
		
					#put coord of proteins back in the box
					coord_dt = coords_remove_whole(proteins2_sele[p2_index].BB.coordinates(), box_size)

					#calculate CoG
					CoG2_dt = calculate_cog(coord_dt, box_size)
		
					#center coord in -box_dim/2 +box_dim/2 referential
					coord_dt = coords_center_in_box(coord_dt, CoG2_dt, box_size)

					#store protein coordinates		
					BB_list_2[p2_index, nb_frame, :] = coord_dt
			
					# coordinates of W144 at t+dt
					coord_pt2_dt = [coord_dt[185][0],coord_dt[185][1],0]


					#calculate and store angle
					if nb_frame > 1:
						angle_list2.append(abs_angle_2D_between(coord_pt2[p2_index], coord_pt2_dt))
						distance_list2.append((mda.analysis.distances.distance_array(np.float32(np.array([CoG2_dt[0],CoG2_dt[1],CoG2_dt[2]]).reshape((1,3))), np.float32(np.array([coord_CoG2[p2_index][0],coord_CoG2[p2_index][1],CoG2_dt[2]]).reshape((1,3))), ts.dimensions))[0][0])

						x_sign= np.float32(CoG2_dt[0]-coord_CoG2[p2_index][0]) 
						x_value = (mda.analysis.distances.distance_array(np.float32(np.array([CoG2_dt[0],CoG2_dt[0],CoG2_dt[0]]).reshape((1,3))), np.float32(np.array([coord_CoG2[p2_index][0],CoG2_dt[0],CoG2_dt[0]]).reshape((1,3))), ts.dimensions))[0][0]
						if x_sign < 0.0 : x_value = -x_value

						y_sign= np.float32(CoG2_dt[1]-coord_CoG2[p2_index][1]) 
						y_value = (mda.analysis.distances.distance_array(np.float32(np.array([CoG2_dt[1],CoG2_dt[1],CoG2_dt[1]]).reshape((1,3))), np.float32(np.array([coord_CoG2[p2_index][1],CoG2_dt[1],CoG2_dt[1]]).reshape((1,3))), ts.dimensions))[0][0]
						if y_sign < 0.0 : y_value = -y_value

						coord_list2_x.append(x_value)
						coord_list2_y.append(y_value)

					else:
						angle_list2.append( abs_angle_2D_between(coord_pt2[p2_index], coord_pt2_dt))
						distance_list2.append((mda.analysis.distances.distance_array(np.float32(np.array([CoG2_dt[0],CoG2_dt[1],CoG2_dt[2]]).reshape((1,3))), np.float32(np.array([coord_CoG2[p2_index][0],coord_CoG2[p2_index][1],CoG2_dt[2]]).reshape((1,3))), ts.dimensions))[0][0])

						x_sign= np.float32(CoG2_dt[0]-coord_CoG2[p2_index][0]) 
						x_value = (mda.analysis.distances.distance_array(np.float32(np.array([CoG2_dt[0],CoG2_dt[0],CoG2_dt[0]]).reshape((1,3))), np.float32(np.array([coord_CoG2[p2_index][0],CoG2_dt[0],CoG2_dt[0]]).reshape((1,3))), ts.dimensions))[0][0]
						if x_sign < 0.0 : x_value = -x_value

						y_sign= np.float32(CoG2_dt[1]-coord_CoG2[p2_index][1]) 
						y_value = (mda.analysis.distances.distance_array(np.float32(np.array([CoG2_dt[1],CoG2_dt[1],CoG2_dt[1]]).reshape((1,3))), np.float32(np.array([coord_CoG2[p2_index][1],CoG2_dt[1],CoG2_dt[1]]).reshape((1,3))), ts.dimensions))[0][0]
						if y_sign < 0.0 : y_value = -y_value

						coord_list2_x.append(x_value)
						coord_list2_y.append(y_value)	
						
			

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#
#						Superimpose the 2 numpy arrays and create a mask
#
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------			

	angle_list = angle_list1
	distance_list = distance_list1
	coord_list_x = coord_list1_x
	coord_list_y = coord_list1_y

	#angle_list = angle_list2
	#distance_list = distance_list2
	#coord_list_x = coord_list2_x
	#coord_list_y = coord_list2_y


	avg_angle = np.average(angle_list)
	avg_distance = np.average(distance_list)


	avg_coord_x = np.average(coord_list_x)
	avg_coord_y = np.average(coord_list_y)

	max_dist = np.max(distance_list)+5
	max_angle = np.max(angle_list)+5


	min_coord_x = np.min(coord_list_x)-5
	max_coord_x = np.max(coord_list_x)+5

	min_coord_y = np.min(coord_list_y)-5
	max_coord_y = np.max(coord_list_y)+5

	color_figure = 'green'
	#_figure = 'orange'

	print("global avg. angle: "+str(avg_angle)) 
	plot_distrib_angle(angle_list,'angle_distrib',dt, max_angle, 0.4,color_figure)
	print("  ")
	print("global avg. distance: "+str(avg_distance))
	plot_distrib_distance(distance_list,'distance_distrib',dt, max_dist, 0.4,color_figure)
	print("-----------------------------")
	print("  ")
	print("global avg. coord x: "+str(avg_coord_x))
	plot_distrib_coord(coord_list_x,'coord_distrib_x',dt, min_coord_x, max_coord_x, 0.4,color_figure)
	print("global avg. coord y: "+str(avg_coord_y))
	plot_distrib_coord(coord_list_y,'coord_distrib_y',dt, min_coord_y, max_coord_y, 0.4,color_figure)
	print("-----------------------------")
	print("-----------------------------")
	print("  ")
	print("  ")


	matrix_angle = np.asarray(angle_list)
	matrix_distance = np.asarray(distance_list)
	matrix_coord_x = np.asarray(coord_list_x)
	matrix_coord_y = np.asarray(coord_list_y)

	print(len(matrix_angle))

	matrix_angle.shape = shape_matrix_cluster
	matrix_distance.shape = shape_matrix_cluster
	matrix_coord_x.shape = shape_matrix_cluster
	matrix_coord_y.shape = shape_matrix_cluster



	list_cluster_size = []

	list_distance_mean = []
	list_distance_spread = []

	list_angle_spread = []


	list_coord_x_mean = []
	list_coord_x_spread = []

	list_coord_y_mean = []
	list_coord_y_spread = []

	for i in range(1,max_cluster):

	
		# create the boolean matrix
		bool_matrix = matrix_clust == i

		# extract the good values from the angle and distance based on the boolean matrix		
		matrix_angle_cluster = matrix_angle[bool_matrix]
		matrix_distance_cluster = matrix_distance[bool_matrix]

		matrix_coord_x_cluster = matrix_coord_x[bool_matrix]
		matrix_coord_y_cluster = matrix_coord_y[bool_matrix]

		list_angle_cluster = matrix_angle_cluster.flatten().tolist()
		list_distance_cluster = matrix_distance_cluster.flatten().tolist()

		list_coord_x_cluster = matrix_coord_x_cluster.flatten().tolist()
		list_coord_y_cluster = matrix_coord_y_cluster.flatten().tolist()

		if len(list_angle_cluster) != 0 : 

			list_cluster_size.append(i)

			avg_angle_cluster = np.average(list_angle_cluster)
			avg_distance_cluster = np.average(list_distance_cluster)
			avg_coord_x_cluster = np.average(list_coord_x_cluster)
			avg_coord_y_cluster = np.average(list_coord_y_cluster)

			max_dist = np.max(list_distance_cluster)+5
			max_angle = np.max(list_angle_cluster)+5

			min_coord_x = np.min(list_coord_x_cluster)-5
			max_coord_x = np.max(list_coord_x_cluster)+5
			min_coord_y = np.min(list_coord_y_cluster)-5
			max_coord_y = np.max(list_coord_y_cluster)+5



			print("avg. angle cluster "+str(i)+" : "+str(avg_angle_cluster)) 
			name_file_angle = 'angle_distrib_cluster'+str(i)
			angle = plot_distrib_angle(list_angle_cluster,name_file_angle,dt, max_angle, 0.4,color_figure)
			list_angle_spread.append(angle)			
			print("  ")
			print("avg. distance cluster "+str(i)+" : "+str(avg_distance_cluster))
			name_file_dist = 'distance_distrib_cluster'+str(i) 
			distance = plot_distrib_distance(list_distance_cluster,name_file_dist,dt, max_dist, 0.4,color_figure)
			list_distance_mean.append(distance[0])	
			list_distance_spread.append(distance[1])		
			print("-----------------------------")
			print("  ")	

			print("global avg. coord cluster x "+str(i)+" : "+str(avg_coord_x_cluster))
			name_file_coord_x = 'distance_distrib_coord_x_cluster'+str(i) 
			coord_tmp = plot_distrib_coord(list_coord_x_cluster,name_file_coord_x,dt, min_coord_x, max_coord_x, 0.4,color_figure)
			list_coord_x_mean.append(coord_tmp[0])	
			list_coord_x_spread.append(coord_tmp[1])

			print("global avg. coord cluster y "+str(i)+" : "+str(avg_coord_x_cluster))
			name_file_coord_y = 'distance_distrib_coord_y_cluster'+str(i) 
			coord_tmp = plot_distrib_coord(list_coord_y_cluster,name_file_coord_y,dt, min_coord_y, max_coord_y, 0.4,color_figure)
			list_coord_y_mean.append(coord_tmp[0])	
			list_coord_y_spread.append(coord_tmp[1])

			print("-----------------------------")
			print("  ")
			print("  ")

	figure_distance_ms = matplotlib.pyplot.figure()
	ax_distance_ms = figure_distance_ms.add_subplot(1,1,1)
	matplotlib.pyplot.errorbar(list_cluster_size,list_distance_mean,yerr= list_distance_spread,fmt='-o')
	matplotlib.pyplot.savefig('figure_distance_ms.svg')

	figure_angle_ms = matplotlib.pyplot.figure()
	ax_angle_ms = figure_angle_ms.add_subplot(1,1,1)
	matplotlib.pyplot.plot(list_cluster_size,list_angle_spread,'b-')
	matplotlib.pyplot.savefig('figure_angle_ms.svg')

	figure_coord_ms = matplotlib.pyplot.figure()
	ax_coord_ms = figure_coord_ms.add_subplot(1,1,1)
	matplotlib.pyplot.errorbar(list_cluster_size,list_coord_x_mean,yerr= list_coord_x_spread,fmt='b-')
	matplotlib.pyplot.errorbar(list_cluster_size,list_coord_y_mean,yerr= list_coord_y_spread,fmt='g-')
	matplotlib.pyplot.savefig('figure_coord_ms.svg')


# -----------------------------------------------------------------------------------------------
#
#
#                  Functions to be called 
#
# -----------------------------------------------------------------------------------------------

def plot_distrib_coord (list_values,name,dt,xmin, xmax,ymax,color):

	#print(list_values)

	fig = matplotlib.pyplot.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xlim([xmin,xmax])
	n, bins, patches = matplotlib.pyplot.hist(list_values, 50, normed=1, facecolor=color, alpha=0.5)

	param = norm.fit(sorted(list_values))
	pdf_fitted = norm.pdf(sorted(list_values),loc=param[0],scale=param[1])

	mean_norm = norm.mean(loc=param[0],scale=param[1]) 
	std_norm = norm.std(loc=param[0],scale=param[1])/math.sqrt(float(len(list_values)))

	print('mean = '+str(mean_norm))
	print('std error = '+str(std_norm))

	#chi_value = chisquare(sorted(list_values), f_exp=pdf_fitted)
	#print('chi_square = '+str(chi_value[0]))
	#print('p_value = '+str(chi_value[1]))
	
	matplotlib.pyplot.plot(sorted(list_values),pdf_fitted,'r-')

	matplotlib.pyplot.savefig(str(name)+'_'+str(dt)+'.svg')
	#matplotlib.pyplot.show()
	matplotlib.pyplot.close()
	return (mean_norm,std_norm)



def plot_distrib_distance (list_values,name,dt,xmax,ymax,color):

	#print(list_values)

	fig = matplotlib.pyplot.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xlim([0,xmax])
	n, bins, patches = matplotlib.pyplot.hist(list_values, 50, normed=1, facecolor=color, alpha=0.5)

	param = rayleigh.fit(sorted(list_values))
	pdf_fitted = rayleigh.pdf(sorted(list_values),loc=param[0],scale=param[1])

	mean_rayleigh = rayleigh.mean(loc=param[0],scale=param[1]) 
	std_rayleigh = rayleigh.std(loc=param[0],scale=param[1])/math.sqrt(float(len(list_values)))

	print('mean = '+str(mean_rayleigh))
	print('std error = '+str(std_rayleigh))

	#chi_value = chisquare(sorted(list_values), f_exp=pdf_fitted)
	#print('chi_square = '+str(chi_value[0]))
	#print('p_value = '+str(chi_value[1]))
	
	matplotlib.pyplot.plot(sorted(list_values),pdf_fitted,'g-')

	matplotlib.pyplot.savefig(str(name)+'_'+str(dt)+'.svg')
	#matplotlib.pyplot.show()
	matplotlib.pyplot.close()
	return (mean_rayleigh,std_rayleigh)

def plot_distrib_angle (list_values,name,dt,xmax,ymax,color):

	#print(list_values)

	fig = matplotlib.pyplot.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xlim([0,xmax])
	n, bins, patches = matplotlib.pyplot.hist(list_values, 50, normed=1, facecolor=color, alpha=0.5)

	param = halfnorm.fit(sorted(list_values))
	pdf_fitted = halfnorm.pdf(sorted(list_values),loc=param[0],scale=param[1])

	mean_halfnorm = halfnorm.mean(loc=param[0],scale=param[1]) 
	std_halfnorm = halfnorm.std(loc=param[0],scale=param[1])

	print('mean = '+str(mean_halfnorm))
	print('std  = '+str(std_halfnorm))

	#chi_value = chisquare(sorted(list_values), f_exp=pdf_fitted)
	#print('chi_square = '+str(chi_value[0]))
	#print('p_value = '+str(chi_value[1]))

	
	matplotlib.pyplot.plot(sorted(list_values),pdf_fitted,'g-')

	matplotlib.pyplot.savefig(str(name)+'_'+str(dt)+'.svg')
	#matplotlib.pyplot.show()
	matplotlib.pyplot.close()
	return std_halfnorm


def bin_array(array):
	
	max_array = np.max(array)
	#print max_array
	array_out =[]	
	
	for i in range(np.shape(array)[0]):

		count = []		
		for j in range(1,int(max_array)+1):
			value = np.count_nonzero(array[i] == j)
			count.append(value/j)

		array_out.append(count)
	
	array_out = np.asarray(array_out)

	return array_out

def plot_diagram(array1, array2, nb_plots, timeline):

	max_array = np.max(array1)
	x = np.arange(0,np.shape(array2)[1])

	for i in range(0,np.shape(array2)[0],int(np.shape(array2)[0]/nb_plots)):
		
		max_line = np.max(array2[i])+5
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax.set_xlim([0,max_array])		
		ax.set_ylim([0,max_line])
		#print list(array[i])
		#print list(x)
		ax.bar(list(x),list(array2[i]),color='#9999ff') 
		plt.savefig(str(timeline[i])+'.svg')
		#plt.show()
		plt.close()




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
    parser.add_argument('--clustfile', type=str, help='path of clusterfile', default=1,
                        dest='clustfile')
    args = parser.parse_args()
    main(args.coord, args.traj, args.proteins1_nb, args.proteins2_nb, args.index_prot1, args.index_prot2, args.clustfile)


# mix
# python ../../scripts/bin_rot_disp.py  --num-prot1 72 --num-prot2 72 --index_prot1 1332 --index_prot2 2208  --coord ../pbc_fix/start_now_prot.gro --traj ../pbc_fix/now_20micros_dt1000_fixed_t1_prot.xtc --clustfile ../cluster/cluster_t1/4_clusters_sizes/4_1_plots_2D/xvg/4_1_clusterprot_2D.stat 

# OmpF
# python ../../scripts/bin_rot_disp.py  --num-prot1 0 --num-prot2 100 --index_prot1 1332 --index_prot2 2208 --coord ../pbc_fix/start_now_prot.gro --traj ../pbc_fix/now_20micros_dt1000_fixed_t1_prot.xtc --clustfile ../cluster/cluster_t1/4_clusters_sizes/4_1_plots_2D/xvg/4_1_clusterprot_2D.stat 

# BtuB



