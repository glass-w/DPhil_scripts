#import seaborn
import MDAnalysis as mda
import numpy as np
#from MDAnalysis.analysis import distances
from MDAnalysis.lib.distances import distance_array # faster C lib
import matplotlib.pyplot as plt
import glob, os
from os.path import join
from datetime import datetime

print ' '
print "### Modules Used ####"
print ' '
print "NumPy Version: " + np.__version__
print "MDAnalysis Version: " + mda.__version__
print ' '
print "####################"
print ' '

u = mda.Universe("restraints_31.25.gro", "b3tri_WB_200ns_skip10_noPBCComplete.xtc")

num_protein = 3
cut_off_value = 15
protein_selection = u.select_atoms("protein")

def make_hmap(universe):

    ''' makes an empty numpy array to be populated
     universe = a single frame or whole trajectory
     resid = list of residue ids
     resnames = list of residue names
     '''

    resid_list = np.asarray(u.select_atoms("protein and not backbone").atoms.residues.resids) # convert list to array
    atom_list = np.asarray(u.select_atoms("protein").positions)

    #resname_list = np.asarray(u.select_atoms("protein and not backbone").atoms.residues.resnames)

    residue_heat_map = np.zeros((len(resid_list), len(resid_list))) # define the array size before making the heatmap

    atom_heat_map = np.zeros((len(atom_list), len(atom_list)))


    return resid_list, residue_heat_map, atom_heat_map

def populate_hmap(empty_heatmap, universe, resid):

    '''populates empty heatmap array for a single frame
    heatmap = empty numpy array of N x N
    universe = a single frame or whole trajectory
    resid = list of residue ids
    amino_acid = amino acid in the list
    '''

    heat_map_pop = empty_heatmap    # reset for the next cycle

    com_store = np.zeros((len(resid), 3))

    #for item in range(len(resid)):  # For each residue, calculate the c.o.m and store

        #com_store[item] = universe.select_atoms("resid " + str(resid[item])).center_of_mass()

    # Calculate a distance matrix

    #d = distance_array(np.array(com_store).astype(np.float32), np.array(com_store).astype(np.float32), box=u.dimensions,
                       #backend='OpenMP')

    # Define a new heat map based on a user defined cut-off and convert to 1's & 0's

    #heat_map_pop += (d < cut_off_value).astype(int)

    d = distance_array(np.array(universe.select_atoms("protein").positions).astype(np.float32), np.array(universe.select_atoms("protein").positions).astype(np.float32),\
                       box=u.dimensions, backend='OpenMP')

    heat_map_pop += (d < cut_off_value).astype(int)

    return heat_map_pop

def plot_hmap(data, residues):
    '''plot a single frame
    data = the populated hmap
    names = list of amino acid names
    resids = list of resids of amino acids
    '''

    data = (data / np.amax(data)) * 100  # normalise to make a percentage

    fig = plt.figure()

    ax = fig.add_subplot(111)

    #colors = cm.ScalarMappable(cmap="viridis").to_rgba(data_value)



    #surf = ax.plot_surface(x,y,z,rstride=1, cstride=1,linewidth=0, antialiased=True)

    plt.suptitle('Contacts < ' + str(cut_off_value) + ' $\AA^2$')
    plt.xlabel("Residue No.")
    plt.ylabel("Residue No.")

    plt.imshow(data, cmap='plasma')#, interpolation='spline16')

    resid_labels = range(1, (len(residues) / num_protein))

    # Add lines to see where protein repeats start and stop
    plt.axvline((len(residues)/num_protein), c='white', linestyle='--', lw=0.4)
    plt.axvline((len(residues)/num_protein)*(num_protein - 1), c='white', linestyle='--', lw=0.4)
    plt.axhline((len(residues) / num_protein), c='white', linestyle='--', lw=0.4)
    plt.axhline((len(residues) / num_protein) * (num_protein - 1), c='white', linestyle='--', lw=0.4)

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Interaction percentage over trajectory', rotation=270)
    cbar.ax.get_yaxis().labelpad = 20

    #plt.show()

    return fig

def assign_beta_values(universe, residue_ids): #### STILL IN TEST PHASE ####

    residues_per_protein = (len(residue_ids) / num_protein)

    atoms_per_protein = (len(universe.select_atoms("protein").positions) / num_protein)

    #interface_array = heat_map_pop[0:(length_of_protein * 2), length_of_protein:(length_of_protein * 3) ] # selects square array of interactions off-diagonal

    #prot_1_2_array = heat_map_pop[0:(residues_per_protein), (residues_per_protein):((residues_per_protein) * 2)]

    #prot_1_3_array = heat_map_pop[0:(residues_per_protein), (residues_per_protein * 2):(residues_per_protein * 3)]

    #
    # prot_2_3_array = heat_map_pop[(residues_per_protein):((residues_per_protein) * 2), (residues_per_protein * 2):(residues_per_protein)]

    prot_atoms_1_2_array = heat_map_pop[0:(atoms_per_protein), (atoms_per_protein):((atoms_per_protein) * 2)]

    prot_atoms_1_3_array = heat_map_pop[0:(atoms_per_protein), (atoms_per_protein * 2):(atoms_per_protein * 3)]

    prot_atoms_2_3_array = heat_map_pop[(atoms_per_protein):((atoms_per_protein) * 2),
                     (atoms_per_protein * 2):(atoms_per_protein * 3)]



    # Calculate the mean of each interaction matrix, then loop through these arrays and assign each element to a number
    # The number is equivalent to the residue number

    for i,b in enumerate(np.concatenate([prot_atoms_1_2_array.mean(0), prot_atoms_2_3_array.mean(0), prot_atoms_1_3_array.mean(0)])):

        protein_selection.select_atoms("bynum " + str(i)).set_bfactors(b)

    #for i in range(length_of_protein):

        #if 0 <= i <= (length_of_protein - 1):

            # Get average interfaces contacts on protein 1

        #protein_selection.select_atoms("resnum " + str(i + 1)).set_bfactors(prot_1_2_array.mean(0)[i])

        #protein_selection.select_atoms("resnum " + str((2 * length_of_protein) + i)).set_bfactors(prot_1_3_array.mean(0)[i])

        #for j in range(length_of_protein):

            #print "j = " + str(j)

            #if (0 <= i <= length_of_protein - 1) and (0 <= j <= (length_of_protein - 1)) and interface_array[i, j] >= 1: # if in top left quad

                # set average across rows of interaction?

                #protein_selection.select_atoms("resnum " + str(i + 1)).set_bfactors(interface_array.mean(1)[i])



                # colour protein 1 & 2

                # p1
                #protein_selection.select_atoms("resnum " + str(i + 1)).set_bfactors(interface_array[i, j])
                # p2
                #protein_selection.select_atoms("resnum " + str((length_of_protein) + j)).set_bfactors(interface_array[i, j])

            #if (0 <= i <= length_of_protein - 1) and (length_of_protein <= j <= ((length_of_protein - 1) * 2)) and interface_array[i, j] >= 1: #if in top right quad

                # colour protein 3 & 1

                # p3
                #protein_selection.select_atoms("resnum " + str((length_of_protein * 2) + i)).set_bfactors(interface_array[i, j])
                # p1
                #protein_selection.select_atoms("resnum " + str(j + 1)).set_bfactors(interface_array[i, j])

            #if (length_of_protein <= i <= (length_of_protein - 1) * 2) and (length_of_protein <= j <= ((length_of_protein - 1) * 2)) and interface_array[i, j] >= 1: #if in bottom right quad

                # colour protein 2

                #protein_selection.select_atoms("resnum " + str(length_of_protein + i)).set_bfactors(interface_array[i, j])
                #protein_selection.select_atoms("resnum " + str(length_of_protein + j)).set_bfactors(interface_array[i, j])

    protein_selection.write("test_bfact_all.pdb", format="PDB")

if __name__ == "__main__":

    startTime = datetime.now()

    residue_ids, initial_resiude_heat_map, initial_atom_heat_map = make_hmap(universe=u)

    for frame in u.trajectory:

        if u.trajectory.frame % 10 == 0:
            print "Frame = " + str(u.trajectory.frame)

        u.trajectory[frame.frame]

        heat_map_pop = populate_hmap(empty_heatmap=initial_atom_heat_map, universe=u, resid=residue_ids)

        #if u.trajectory.frame == 10:
            #break

    u.trajectory[-1] # reset the trajectory prior to writing the beta values
    assign_beta_values(universe=u, residue_ids=residue_ids)


    identify_residues = heat_map_pop.astype(str)

    fig = plot_hmap(data=heat_map_pop, residues=residue_ids)

    fig.savefig('figure.svg', format='svg')
    fig.clf()

    print "Time taken = " + str(datetime.now() - startTime)
