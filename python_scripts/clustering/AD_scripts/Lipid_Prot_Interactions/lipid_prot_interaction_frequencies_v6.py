'''Scheme of the script:
1. identify residue IDs for all repeats of protein (unecessary if only 1 protein obv)
?1a.  Create trajectory consisting of only relevant beads?
2. For each protein residue (or bead of residue) of interest:
	a. Calculate distances between lipid beads and protein residue
	b. If distance less than 0.7 and  add to total for the protein

Changes to v5 from v4 - re-write structure of count_frequencies algorithm

VERSION 6
*** Updates to MDAnalysis Version 0.17.0: ***
***				W. Glass 04/18 			  ***
***			william.glass@chem.ox.ac.uk   ***

The script (version 5) was updated with:
	MDAnalysis 0.17.0
	NumPy 1.11.3
	MatPlotLib 2.2.2

The underlying structure has been kept more or less the same. The main differences come from selection of atoms using
MDAnalysis 0.17.0.

	Count Frequencies Function:
		1) The distance_array function is loaded from a (slightly) different library than in 0.8.0.
		2) u.select_atoms is used instead of u.segments for protein selection.
		3) There is a small discrepancy (detailed in the function) b/w MDAnalysis version 0.8.0 and 0.17.0 in what lipids
		are considered to be "close", however when viewing the structure in VMD the beta values are essentially
		indistinguishable - may be something to test...
		4) The way in which distances are calcualted are now for EACH ATOM, i.e. not each residue (as I think this is
		what u.segments did. The conclusion is that the distance calculation is a lot faster than before).


	Show Frequencies On Structure Function:
		1) Due to the way in which the distances are calculated and the way in which MDAnalusis 0.17.0 assigns b_factors
		it is no longer possible to assign b_factors to each residue, instead you must assign to atoms. To get around
		this a three letter dictionary is referenced and the corresponding rows in the coordiantes assigned the same
		b_factor (called tempfactor in 0.17.0).

*********************************************
'''

import sys
import re
import numpy as np
import time
import MDAnalysis as mda
#from MDAnalysis.core.parallel.distances import distance_array
from MDAnalysis.lib.distances import distance_array
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cmx
import matplotlib.colors as colors

import ast
import argparse

print "##########################\n"
print "Using the following versions: \n"
print "MDAnalysis Version = {}".format(str(mda.__version__))
print "NumPy Version = {}".format(str(np.__version__))
print "MatPlotLib Version = {}\n".format(str(matplotlib.__version__))
print "##########################\n"

from interaction_freq_core import protein_info, get_real_resnums_strlist, bondini_residuesymbol_size, martini_residuesymbol_size

def find_prot_residues(res_total, residue_list, nrepeats):

	'''res_total is the total number of residues in one repeat, residue list should be integers
	 this function will find res IDs of all the repeats in the system,
	 and write selection string for them'''

	residue_dictionary = dict( zip(residue_list, [[] for i in residue_list]) )

	for residue in residue_list:

		for repeat in range(nrepeats):

			residue_dictionary[residue].append(int(residue) + repeat*res_total - 1) # the minus 1 is because residue indexes start at 0

	return residue_dictionary

def make_split_list_single_prot(prot_seq, ff_ressize_dict=martini_residuesymbol_size):

	splitlist_single = []

	counter = 0

	for res in prot_seq:

		counter += ff_ressize_dict[res]

		splitlist_single.append(counter)

	# Changed from version 5, use [:-1] to return all of the elements of the splitlist_single list.
	return splitlist_single[:-1]

def get_protresname_array():

	pass # could add this function later, to make barplots functions a bit more generic/smoother

lipid_particles = {
	'headgroup' : {'POP2': ['CP','P1','P2','C1','C2','C3'],
			'PI3': ['PO3','PO0','PO1','PO2','RP1','RP2','RP3'],
			'PSPI' : ['PO3','RP1','RP2','RP3'],
			'CHOL' : ['ROH', 'R1', 'R2', 'R3', 'R4', 'R5', 'C1', 'C2'], # this is actually the whole CHOL},
			'DPG3': ['GM1','GM2','GM3','GM4','GM5','GM6','GM13','GM14','GM15','GM16','GM17', 'AM1', 'AM2', 'T1A'],
			'DPSM' : ['PO4','NC3'],
			'POPS' : ['PO4','CN0'], #changed from CNO
			'POPC' : ['PO4','NC3'],
			'POPE' : ['PO4','NH3']},
	'phosphate' : {'PIP2': ['PO3','PO1','PO2'],
			'PI3': ['PO3','PO0','PO1','PO2'],
			'PSPI' : ['PO3'],
			'CHOL' : ['ROH'],
			'GM3': ['B1A','B2A','B3A','INA'], # This is obviously not a phosphate group, but it's in the position of one in the memb
			'PPCS' : ['PO4'],
			'POPS' : ['PO4'],
			'POPC' : ['PO4'],
			'POPE' : ['PO4']},
	'phosphate_singlebead' :{'POP2': ['CP'],# {'PIP2': ['PO3'],
			'PI3': ['PO3'],
			'PSPI' : ['PO3'],
			'CHOL' : ['ROH'],
			'GM3': ['INA'], # This is obviously not a phosphate group, but it's in the position of one in the memb
			'PPCS' : ['PO4'],
			'POPS' : ['PO4'],
			'POPC' : ['PO4'],
			'POPE' : ['PO4']}}

lipid_particles_older_vers = {
	'headgroup' : {'PIP2': ['PO3','PO1','PO2','RP1','RP2','RP3'],
			'PI3': ['PO3','PO0','PO1','PO2','RP1','RP2','RP3'],
			'PSPI' : ['PO3','RP1','RP2','RP3'],
			'CHOL' : ['ROH', 'R1', 'R2', 'R3', 'R4', 'R5', 'C1', 'C2'], # this is actually the whole CHOL},
			'GM3': ['INA','B1A','B2A','B3A','B1B','B2B','B3B','INB','B1C','B2C','B3C','INC','B4C','B5C'],
			'PPCS' : ['PO4','NC3'],
			'DPSM' : ['PO4','NC3'],
			'POPS' : ['PO4','CNO'],
			'POPC' : ['PO4','NC3'],
			'POPE' : ['PO4','NH3']},
	'phosphate' : {'PIP2': ['PO3','PO1','PO2'],
			'PI3': ['PO3','PO0','PO1','PO2'],
			'PSPI' : ['PO3'],
			'CHOL' : ['ROH'],
			'GM3': ['B1A','B2A','B3A','INA'], # This is obviously not a phosphate group, but it's in the position of one in the memb
			'PPCS' : ['PO4'],
			'POPS' : ['PO4'],
			'POPC' : ['PO4'],
			'POPE' : ['PO4']},
	'phosphate_singlebead' :{'PIP2': ['PO3'],
			'PI3': ['PO3'],
			'PSPI' : ['PO3'],
			'CHOL' : ['ROH'],
			'GM3': ['INA'], # This is obviously not a phosphate group, but it's in the position of one in the memb
			'PPCS' : ['PO4'],
			'POPS' : ['PO4'],
			'POPC' : ['PO4'],
			'POPE' : ['PO4']}}

### main function to count interactions and produce interactions file ###
def count_frequencies(gro_file, trr_file, prot_seq, protein_residue_list, nrepeats, lipid, stride, lipid_part, protein_centre, protein_centre_cutoff, cutoff):

	print "Using a protein cut off of {} and protein-lipid cut off of {} \n".format(protein_centre_cutoff, cutoff)

	u = mda.Universe(gro_file, trr_file)

	protein_res_total = len(prot_seq)
	protein_residue_dictionary = find_prot_residues(protein_res_total, protein_residue_list, nrepeats)
	prot_split = make_split_list_single_prot(prot_seq)
	
	lipid_selection = 'resname {} and ('.format(lipid)

	for bead in lipid_particles[lipid_part][lipid]:

		lipid_selection += 'name {} or '.format(bead)

	lipid_selection = lipid_selection[:-4] + ')'
	
	lipid_rep_selection = 'resname {} and name {}'.format(lipid, lipid_particles[lipid_part][lipid][0]) # ie. just choose one bead  ## DPSM and name PO4

	lipids = u.select_atoms(lipid_selection)

	n_lipid_beads = len(lipid_particles[lipid_part][lipid])
	n_lipids = lipids.atoms.n_atoms / n_lipid_beads

	lipid_reps = u.select_atoms(lipid_rep_selection)

	#initialise protein-lipid interactions frequency list
	proteinres_lipid_interactions = np.array([0 for i in protein_residue_list])

	startTime = time.time()    
	print 'Here we go...\n'

	for ts in u.trajectory[::stride]:

		for i in range(nrepeats):

			single_prot = u.select_atoms('protein and resid ' + str(i * protein_res_total) + ':' + str((i + 1) * protein_res_total))

			print "Analysing protein {}, (residues {}:{})".format((i + 1), str(i * protein_res_total), str((i + 1) * protein_res_total))

			# find protein centroid - or pick out residue to represent protein position
			if protein_centre == 'centroid':
				single_prot_cent = np.array([single_prot.centroid()])

			#TODO: fix (i.e. update to MDA 0.17.0) elif statement for if you are interested in a particular residue.
			elif protein_centre == int:
				#pick out BB of specified residue
				single_prot_cent = u.segments[0][i*protein_res_total + protein_centre - 1][0] # -1 is because res numbers are zero-indexed

			else:
				print 'Error: protein_centre should either be an integer residue number, or "centroid"'
				return None

			# find lipids within N AA of prot centroid and pick those out from lipids selection

			close = distance_array(single_prot_cent, lipid_reps.positions, ts.dimensions) < protein_centre_cutoff

			lipids_close_indices = []

			# n.b. slight discrepancy between this version and the previous MDA 0.8.0 version, need to investigate this
			# 	however it does give basically the same result. This arises from the total No. lipids considered to be
			# 	close (uncomment the print statements below to see). This may be due to the different distance_array funcs
			# 	used b/w versions (?).

			# see https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.nonzero.html
			for index in np.nonzero(close)[1]: # np.nonzero gives a tuple of arrays - the second gives the indices of lipid residues that are 'close' to the protein

				lipids_close_indices += range(index*n_lipid_beads, (index + 1)*n_lipid_beads) # convert residue IDs to atoms IDs for 'close' lipids

			lipids_close = lipids[lipids_close_indices] # stores atom info of lipid heads

			n_lipids_close = sum(close.flatten())

			# now look at prot-lipid interaction on per residue level

			if n_lipids_close == 0:

				continue # ie. nothing is added to proteinres_lipid_interactions

			else:

				all_dists = distance_array(single_prot.atoms.positions, lipids_close.atoms.positions, ts.dimensions)

				protein_lipid_dist_perresidue_all = np.array([[x.min() for x in np.split(lip, prot_split, axis=0)] for lip in np.split(all_dists, n_lipids_close, axis=1)])

				#### Below was for testing, leaving here for now ####
				# protein_lipid_dist_perresidue_all = np.array(
				# 	[[x.min() for x in np.split(lip, len(u.select_atoms("name BB"))/args.nrepeats, axis=0)] for lip in
				# 	 np.split(all_dists, n_lipids_close, axis=1)])

				#### Uncomment for loops below to check what the above is doing: ####

				# for i, lip in enumerate(np.split(all_dists, n_lipids_close, axis=1)):
                #
				# 	print i, lip
                #
				# 	for i, x in enumerate(np.split(lip, prot_split, axis=0)):
                #
				# 		print x
				# 		print x.min()

				interactions = protein_lipid_dist_perresidue_all <= cutoff

				proteinres_lipid_interactions += np.sum(interactions, axis=0) # sum over all lipids for each res

		print 'Frame {} (of {}) took {:3f} s'.format(ts.frame, len(u.trajectory), time.time() - startTime)
		startTime = time.time()

	proteinres_lipid_interactions_dict = dict(zip(protein_residue_list, proteinres_lipid_interactions))

	return proteinres_lipid_interactions_dict

### functions to produce figures from data ###
def plot_frequencies(protein_residue_list, proteinres_lipid_interactions):

	plt.bar(protein_residue_list, [proteinres_lipid_interactions[res] for res in protein_residue_list])

	plt.show()

def plot_frequencies_fancy(protein_residue_list, proteinres_lipid_interactions_dict, prot, cutoff, cutoff_value):

	protein_residue_array = np.array(protein_residue_list)

	actual_res_namenums = np.core.defchararray.add(np.array(list(protein_info[prot]['protein_seq'])), get_real_resnums_strlist(prot))

	def percent_interactions(proteinres_lipid_interactions):

		total_interactions = sum(proteinres_lipid_interactions.values())

		interactions_list = []

		for res in protein_residue_list:

			interactions_list.append(float(proteinres_lipid_interactions[res]))

		interactions_array = np.array(interactions_list)

		percent_interactions = (interactions_array / total_interactions) * 100

		print percent_interactions

		return percent_interactions

	percent_interactions_dict = {}
	lipid_list = sorted(proteinres_lipid_interactions_dict.keys())

	for lipid in lipid_list:

		percent_interactions_dict[lipid] = percent_interactions(proteinres_lipid_interactions_dict[lipid])

	if cutoff:

		# graph will include ANY residue that has more than <<cutoff>> interactions for ANY lipid
		over_cutoff = np.zeros(shape=len(protein_residue_list), dtype=bool)

		for lipid in lipid_list:
			over_cutoff += percent_interactions_dict[lipid] >= cutoff_value*100

		for lipid in lipid_list:

			print percent_interactions_dict[lipid][over_cutoff]

			percent_interactions_dict[lipid] = percent_interactions_dict[lipid][over_cutoff]

		print percent_interactions_dict

		actual_res_namenums = actual_res_namenums[over_cutoff]

	fig, ax = plt.subplots()
	#rects=ax.barh(range(len(percent_interactions)), percent_interactions, align='center' , color='yellow')
	#ax.invert_xaxis()
	ax.invert_yaxis()
	width = 0.7 / len(lipid_list)
	edge_colors = ['orange', 'blue', 'darkgrey']
	i = 0
	lipid_rects = {}

	print lipid_list

	# set colours for lipids
	for lipid in lipid_list:

		# if re.match(r'PI.*', lipid):
		if re.match(r'POP.*', lipid):
			#barcolor = 'yellow'
			#maxpercent = max(percent_interactions_dict[lipid])
			#lipid_cmap = plt.get_cmap('hot_r')
			#lipid_cmap = plt.get_cmap('inferno_r')
			lipid_cmap = colors.LinearSegmentedColormap.from_list('bluemap', ['white', 'darkblue'])
		elif re.match(r'CHOL.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('cyanmap', ['white', 'darkcyan'])
		elif re.match(r'POPS.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('bluemap', ['white','darkblue'])
		elif re.match(r'DPG3.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('magentamap', ['white', 'darkmagenta'])
		elif re.match(r'POPC.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('greymap', ['white','black'])
		elif re.match(r'POPE.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('greymap', ['white','darkgrey'])
		elif re.match(r'PPCS.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('greenmap', ['white','green'])

		# try:

		maxpercent = max(max(percent_interactions_dict[lipid]) for lipid in lipid_list)

		# except ValueError:
		# 	print "\n*** No interactions at your specified -plot_cut value ({}), try another value! ***\n".format(args.plot_cut)
		# 	sys.exit(1)

		cNorm = colors.Normalize(vmin=0,vmax=maxpercent)
		scalar_map = cmx.ScalarMappable(norm = cNorm, cmap = lipid_cmap)
		barcolor=scalar_map.to_rgba(percent_interactions_dict[lipid])
		lipid_rects[lipid] = ax.barh(np.arange(len(percent_interactions_dict[lipid]))+width*i, percent_interactions_dict[lipid], width, align='center' , color=barcolor, ec=edge_colors[i])

		i += 1

	#if len(proteinres_lipid_interactions_dict.keys()) > 1:
	#	ax.legend( lipid_rects.values(), lipid_rects.keys() )
	ax.set_ylabel('Residue')

	# sort out GM3 / PIP2 names for plots
	if lipid_list[0] == 'DPG3':
		ax.set_xlabel('% of total {} contacts'.format('GM3'))
	elif lipid_list[0] == 'POP2':
		ax.set_xlabel('% of total {} contacts'.format('PIP2'))

	else:

		ax.set_xlabel('% of total {} contacts'.format(lipid_list[0]))

	ax.set_yticks(range(len(actual_res_namenums)))
	ax.set_yticklabels(actual_res_namenums, rotation='horizontal')
	plt.axis('tight')
	fig.set_size_inches(3.7, 4.3)
	plt.savefig('prot_{}_int_v6.svg'.format(''.join(lipid_list)), orientation='portrait', bbox_inches='tight')
	#plt.show()

def show_frequencies_on_structure(nmonomers, interactions, gro_file, outfile, sensitive=True):

	#normalise interactions - not really necessary - can do this w/ VMD colorscales
	u = mda.Universe(gro_file)
	protein = u.select_atoms("protein")
	nres = len(protein.residues)

	nresmono = nres/nmonomers

	total_interactions = sum(interactions.values())

	protein_residue_list = sorted(interactions.keys())

	interactions_array = np.array([float(interactions[res]) for res in protein_residue_list])

	if sensitive:
		percent_interactions = interactions_array*10000 / total_interactions # to make smaller numbers visible (for when total_interactions is v large)

	else:
		percent_interactions = interactions_array*100 / total_interactions

	residue_bfactors = np.zeros((len(np.asarray(u.select_atoms("protein").atoms)), 1))

	# this is a work around, it is needed in order to assign b_factors later on
	u.add_TopologyAttr(mda.core.topologyattrs.Tempfactors(np.zeros(len(u.atoms))))

	residue_atoms_storage = []
	residue_resname_storage =[]

	residue_info = []

	for i in range(nresmono): # for each residue in a monomer

		for j in range(nmonomers): # for the number of monomers i.e. if your protein has multiple subunits

			residue_atoms = protein.residues[j*nresmono + i].atoms

			residue_resname_storage.append(protein.residues[j*nresmono + i].resname)

			residue_info.append(protein.residues[j*nresmono + i])

			for particle in range(len(residue_atoms)):

				residue_atoms_storage.append(residue_atoms[particle])

	# set dict based on No. particles in each residue, using the three letter AA codes
	martini_residuesymbol_size_three = {'ALA': 1, 'CYS': 2, 'ASP': 2, 'GLU': 2, 'PHE': 4, 'GLY': 1, 'HIS': 4, 'ILE': 2,
										'LYS': 3, 'LEU': 2, 'MET': 2, 'ASN': 2, 'PRO': 2, 'GLN': 2, 'ARG': 3, 'SER': 2,
										'THR': 2, 'VAL': 2, 'TRP': 5, 'TYR': 4}

	### Due to the way MDA 0.17.0 sets bfactors need to get number of particles per residue. ###
		# Done by referencing a dictionary when going through the residue_atoms_storage array:
		# 	- Set the rows to: current row up to No. particles in the current residue (-1 due to indexing).
		#	- Set all of these rows (particles in a residue) to the value obtained in the interaction dict.

	# set start values, 'start' is used to assign particles, 'i' is used to reference the residue info
	# i.e. residue_info & percent_interactions match in length
	start = 0
	i = 0

	for residue in percent_interactions:

		# get info on current residue
		current_resname = residue_info[i].resname

		current_resnumber = residue_info[i].resid

		num_particles_in_res = martini_residuesymbol_size_three[current_resname]

		# define the end of the current residue based on the No. particles it has (from the dict)
		stop = start + num_particles_in_res

		# assign multiple rows (i.e. particles of the same residue) the same b_factor value
		residue_bfactors[start : stop] = percent_interactions[i]

		# set the start of the next residue
		start = stop

		i += 1

	# assign b_factors
	protein.residues.atoms.tempfactors = residue_bfactors

	# write out to file
	u.select_atoms("protein").write(outfile + ".pdb", format("PDB"))

# functions for reading input - residue list and interactions file

def make_res_list(res_string, protname):
	# residue list should be input in the format 1:4,6:9,300,310:333
	# the following decodes that format into a straight list
	if res_string == 'all':
		res_list = np.arange(len(list(protein_info[protname]['protein_seq']))) + 1
	else:
		res_list_ranges = [x.split(':') for x in res_string.split(',')]
		res_list = []
		for res_range in res_list_ranges:
			res_list += range(int(res_range[0]), int(res_range[-1])+1 )
	return res_list

def make_interactions_dict(dictionary_filename):
	f = open(dictionary_filename)

	# for line in f.readlines():
	# 	print line

	# line = f.readlines()[2]
	dictionary_string = f.readlines()[2].rstrip()
	proteinres_lipid_interactions = ast.literal_eval(dictionary_string)

	#print proteinres_lipid_interactions

	return proteinres_lipid_interactions

### Function to sum data for all (or a selection of) lipids:
def sum_lipid_ints(interactions_dict_list):
	total_interactions = {}
	for interactions_dict in interactions_dict_list:
		for res in interactions_dict.keys():
			if res not in total_interactions:
				total_interactions[res] = interactions_dict[res]
			else:
				total_interactions[res] += interactions_dict[res]
	return total_interactions

if __name__ == "__main__":
	# sub-command funcs for the parser - ie. workflows available when running script:

	def get_ints(args):

		res_list = make_res_list(args.reslist, args.p)
		protein_seq = list(protein_info[args.p]['protein_seq'])

		if args.c != 'centroid':
			args.c = int(args.c) # this should return an error if the argument is not integer-like

			if args.c < len(protein_seq):
				#check that the residue number is valid
				print 'Error: residue number selected is larger than the number of residues in the protein sequence'
				return None

		proteinres_lipid_interactions = count_frequencies(args.f, args.xtc, protein_seq, res_list, args.nrepeats, args.lipid, args.stride, args.lipid_part, args.c, args.cd, args.cutoff)

		f = open('frequencies_dictionary_residues{}_lipid{}{}_v2_5'.format(args.reslist, args.lipid, args.lipid_part), 'w')
		f.write('Generated using the command:\npython {}\n'.format(' '.join(sys.argv)))
		f.write(str(proteinres_lipid_interactions) + '\n')
		f.close()

	def pff(args):

		res_list = make_res_list(args.reslist, args.p)
		proteinres_lipid_interactions_dict = {}

		for i in range(len(args.lipid)):

			if args.lipid[i] not in proteinres_lipid_interactions_dict:
				proteinres_lipid_interactions_dict[args.lipid[i]] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i])

			elif args.lipid[i]+'x' not in proteinres_lipid_interactions_dict:
				proteinres_lipid_interactions_dict[args.lipid[i]+'x'] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i]+'x')

			else: # this will do for 3 systems but need to update if comparing more...
				proteinres_lipid_interactions_dict[args.lipid[i]+'xx'] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i]+'xx')

		#print proteinres_lipid_interactions_dict
		proteinres_lipid_interactions = plot_frequencies_fancy(res_list, proteinres_lipid_interactions_dict, args.p, True, args.plot_cut)

	def structure(args):

		interactions = make_interactions_dict(args.i)
		show_frequencies_on_structure(args.m, interactions, args.f, args.o, args.s)
	
	def sum_ints(args):

		interactions_dict_list = []

		for i in range(len(args.i)):

			interactions_dict_list.append(make_interactions_dict(args.i[i]))

		total_interactions = sum_lipid_ints(interactions_dict_list)

		f = open('frequencies_dictionary_sumlipids_v2_5'.format(), 'w')
		f.write('Generated using the command:\npython {}\n'.format(' '.join(sys.argv)))
		f.write(str(total_interactions) + '\n')
		f.close()
		
		
	parser = argparse.ArgumentParser(description='Analyse prot-lipid interactions at the residue level', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers()

	parser_ints = subparsers.add_parser('get_ints', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_ints.add_argument('-f',help='input gro file', required=True)
	parser_ints.add_argument('-x', '--xtc', help='List of xtc files', nargs='+', required=True)
	parser_ints.add_argument('-reslist',default='all', help='''reslist should be input in the format:  1:4,6:9,300,310:333  (no spaces, residue ranges are inclusive) 
	Residue numbering starts at 1)
	OR ALL residues can be specified using: all''')
	parser_ints.add_argument('-p', help='Protein name', choices=protein_info.keys(), required=True)
	parser_ints.add_argument('-r', '--nrepeats', type=int, help='Number of protein MONOMERS in system', required=True)
	parser_ints.add_argument('-lipid', help='Name of lipid', choices=lipid_particles['headgroup'], required=True)
	parser_ints.add_argument('-stride', type=int, default = 1, help='Frame number intervals that xtc will be read at')
	parser_ints.add_argument('-lp', '--lipid_part', default='headgroup', help='Part of lipid to consider in interactions', choices=lipid_particles.keys())
	parser_ints.add_argument('-c', default='centroid', help='Part of protein to measure lipid distances from, for first approximation - can be either "centroid" for protein centroid, or an integer residue number (residue numbering starts from 1)')
	parser_ints.add_argument('-cd', type = float, default=80, help='Distance away from protein to select lipids that will undergo a closer inspection (a lower number speeds the script up, but too low may mean that some lipids aren\'t considered')
	parser_ints.add_argument('-d', '--cutoff', type=float, default=6.5, help='Interactions cutoff distance')
	parser_ints.set_defaults(func=get_ints)

	parser_pff = subparsers.add_parser('plot_ff', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_pff.add_argument('-reslist',default='all', help='''reslist should be input in the format:  1:4,6:9,300,310:333  (no spaces, residue ranges are inclusive - 
	Residue numbering starts at 1)
	OR ALL residues can be specified using: all''')
	parser_pff.add_argument('-i',help='Interaction file', nargs='+')
	parser_pff.add_argument('-p', help='Protein name', choices=protein_info.keys(), required=True)
	parser_pff.add_argument('-lipid', default='PIP2', help='Name of lipid', choices=lipid_particles['headgroup'], nargs='+')
	parser_pff.add_argument('-plot_cut', default=0.025, type=float, help='Cut off value for bar charts')
	parser_pff.set_defaults(func=pff)

	parser_struct = subparsers.add_parser('structure', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_struct.add_argument('-i', help='Interaction file')
	parser_struct.add_argument('-m',type=int, help='Number of monomers in one protein')
	parser_struct.add_argument('-f',help='input gro file')
	parser_struct.add_argument('-o',help='output gro filename')
	parser_struct.add_argument('-s',help='sensitive - amplifies all interactions in order to pick up the least frequent interactions', action='store_true', default=False)
	parser_struct.set_defaults(func=structure)
	
	parser_addints = subparsers.add_parser('add', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_addints.add_argument('-i', help='Interaction files', nargs = '+', required=True)
	parser_addints.set_defaults(func=sum_ints)
	
	args = parser.parse_args()
	args.func(args)