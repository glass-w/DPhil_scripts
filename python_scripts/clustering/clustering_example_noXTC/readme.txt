Files inside this folder:

1. cluster_prot.py
The script.

2. example_data
Folder containing some example input data and the output obtained from running the script on this data.

#######################
Running the script
#######################

To get help and information about the script, run:
> python clustering_prot.py  -h

The help output gives a lot of information about what all the options are.

For the most basic run, the required inputs are:
-f 	The structure file [.gro] 
-x 	The trajectory file [.xtc]
-o 	Name of output folder

--nx_cutoff 
	Cutoff distance for proteins to be considered 'interacting', measured from protein centroid to protein centroid.  This can be estimated by loading the simulation into VMD and measuring distances between proteins known to be interacting. Units of measurement are angstroms.  If in doubt, overestimate this number, especially if running the script mostly to get information about which residues are interacting.

--species
	Filename of a file that gives the name and sequence of the protein (ie. 'species') in your simulation.  This is necessary to give output about which residues are interacting.  The format of the file is specified in note 6. of the help but, in brief, for a single protein species which is a monomer the file should contain on line with the format: '<protein_name>,1,<protein_sequence>'  The protein sequence should the same as for the protein actually used in your simulation; protein_name is anything you want.  If there is more than one type of protein in the simulation this file would contain as many lines as different types of protein.

ie. to run the script with these inputs use:
> python clustering_prot.py -f structure.gro -x processed_trajectory.xtc -o output_folder_name --species species_info_filename --nx_cutoff cutoff_in_angstroms 

#######################
In example_data are:
#######################
INPUTS:
#######################
Required input files:

input gro file:	md_prod12x12_noW_justphosphates+prot_firstframe_pbccenter.gro
	A structure file containing the same number of atoms as in the trajectory (xtc file).  
To get a .gro file containing only protein molecules:
1. create an index file:
> make_ndx -f <original_gro_file>.gro -o protein.ndx
-> enter 'q' when prompted, to quit  - the protein group should be auto detected so no other input is needed
2. Use editconf to make the .gro file:
> editconf -f <original_gro_file>.gro -o protein_only.gro -n protein.ndx
-> enter '1' when prompted - this should be the index number of the 'Protein' group

input xtc file:	md_prod12x12_noW_justphosphates+prot.pbcmol.xtc
	Trajectory file.  This particular file is of a 10 microsecond trajectory, with frames every 50 ns and was created with the -pbc mol option in trjconv, to keep molecules whole.  The script will run faster if the xtc is filtered such that there are fewer frames (sampling at 50 ns, as here, is probably minimal, and may not be enough if the simulation is shorter than 10 us) and contains fewer atoms.  Here, I have proteins and lipid phosphate atoms, but in fact only protein atoms are necessary.
To get a suitable xtc file from the original xtc file use trjconv:
> trjconv -s <tpr_used_for_simulation>.tpr -f md_prod12x12.xtc -dt 50000 -pbc mol -o md_prod12x12.50ns.onlyprot.pbcmol.xtc
The '-dt' flag picks frames every 50 ns (i.e. 50000 ps) and the '-pbc mol' option keeps molecules whole across periodic boundary conditions. When trjconv runs, you should get the option to select which molecules are included in the output.  It's fine to select group '1' ('Protein') to get fewer atoms (protein only).

information about protein: Kir_info
	This give the sequence and oligomeric state of the protein, and is needed to label the residues in the interacting-residues plot.  In this file, for instance, the name of the protein is given (this can be anything), KIR is made up of four monomers, hence the '4', and then the sequence for one monomer is given.

Optional input:

More information about protein: Kir_COG_info
	The script works out clustering information based on distances between protein centroids.  In the case of KIR, I wanted to specify which protein residues were used to calculate the protein centroid, and that information can be input using such a file.  More details about format are in note 8 of the help.


#######################
OUTPUT
#######################

In the output folder:
Log file: clustering_prot.log
	This file saves the command used to generate the data

There are also a few folders, each containing different types of information:
1_snapshots	
Some information about the system - not very useful

2_proteins_interactions 	
Folder containing all information about which proteins and residues interact with other proteins.  On a residue-to-residue basis this is shown in 2_interactions_residues_A-A_2D.xvg.  2_proteins_interactions.svg and 2_proteins_neighbours show which proteins are interacting with one another and how many neighbours each protein has. The pdb file has gives the data in terms of structure.  To view, load the pdb into VMD and colour by 'Beta'.

3_cluster_compositions
Useful for a system with more than one protein type

4_cluster_sizes
Plots of how proteins cluster - useful for simulations with high numbers of proteins

6_VMD
Information here can be used to map protein cluster size colouring into a video of the trajectory in VMD - see 'VMD visualisation' in the help output.

In the 'example_output' folder I have left two examples of log files from two simulations, mainly to show what commands I used: clustering_prot.log and clustering_prot_6x6.log
