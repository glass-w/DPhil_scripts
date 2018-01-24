'''edited 24/07/2014, to hopefully sort out the problem Heidi had
ie. so that grep-ing for molecule names (around lines 70-80) 
doesn't pick out INA when grep-ing for NA+ - the problem that Heidi had 
was that an atom type with name INA was picked out when grep-ing for 
the molecule type NA+

further edits that would be good:
check for egrep 'special characters' in molecule names, 
and backslash them out

15 Nov 2016
Edits with the help of Sarah-Beth Amos:
- in 'make_new_top_file': add possibility of white space when searching for [ molecules ] line in .top file
- in 'reorder_genconf': add residue renumbering (via call to editconf)
'''
import sys
import re
import os

from resave_old_copies import store_old_copy_gromacs_style
from renumber_gro_file import renumber_atomnumbers_only, renumber_residues

def make_new_top_file(top_filename, n_repeats, genconf_filename):
    reordered_top_file = 'reordered_' + genconf_filename[:-4] + '.top'
    store_old_copy_gromacs_style(reordered_top_file)
    
    f = open(top_filename)
    g = open(reordered_top_file, 'w')
    read = 0
    molecules = []
    nmols = []
    print 'creating large system .top file...'
    for line in f:
        if read == 1 and infoline.match(line):
            molinfo = infoline.match(line)
            molecule = molinfo.group(1).rstrip()
            nmol = int(molinfo.group(2))
            if molecule not in molecules:
                molecules.append(molecule)
                nmols.append(nmol)
            else:
    		   nmols[molecules.index(molecule)] += nmol
        else:
            g.write(line)
        if re.search(r'^\s*\[ molecules \]',line):
            read = 1
            print '\tstarting to read molecule numbers'
            infoline = re.compile(r'(.+)\s+(\d+)')
    f.close()
    for molecule in molecules:
        g.write( '{} {}\n'.format(molecule, (nmols[molecules.index(molecule)])*n_repeats) ) 
    g.close()
    print 'Molecules: ', molecules
    
    new_top_filename = genconf_filename[:-4] + '.top'
    store_old_copy_gromacs_style(new_top_filename)
    os.rename(reordered_top_file, new_top_filename)
    
    return molecules, nmols

#Protein_list = ['OmpF', 'BtuB']
cheat_dictionary = {'OmpF':2208, 'BtuB':1332, 'OmpF_BONDINI':2049, 'BtuB_BONDINI':1148, 'Protein_nav':2993, 'Protein_beta3':381}

def reorder_genconf_file(genconf_filename, molecules, nmols, n_repeats, force_field='martini'):
    #define a new filename for the reordered genconf file:
    reordered_filename = 'reordered_'+genconf_filename
    #check that 'reordered.gro' doesn't already exist; if it does, save the old copy, gromacs style
    store_old_copy_gromacs_style(reordered_filename)

    #sort out molecules list: if there are any proteins in the list 
    # (defined using Protein list)
    # then add the number of molecules together label these as a single 'Protein' class of molecules
    Protein_list = raw_input('Enter the names of protein molecules, as in the .top file, separated by spaces: ').split(' ')
    if Protein_list != ['']:
    #find number of atoms in each protein - this is a cheat for now
    #find number of each protein
    #create variables that can be used to group proteins together
	if len(Protein_list) >1:
	    natoms_in_prot1repeat = {}
	    for prot in Protein_list:
		if force_field == 'martini':
		    natoms_in_1prot = int(cheat_dictionary[prot])
		elif force_field == 'bondini':
		    natoms_in_1prot = int(cheat_dictionary[prot+'_BONDINI'])
		else:
		    print 'Force field should be specified as either bondini or martini (default)'
		    return None
		nprot = int(nmols[molecules.index(prot)])
		natoms_in_prot1repeat[prot] = natoms_in_1prot * nprot
	    total_prot_atoms_1repeat = sum(natoms_in_prot1repeat.values())
	    print 'total_prot_atoms_1repeat: ', total_prot_atoms_1repeat 
	    
    #then change molecules list so that all proteins are named 'Protein'
        for prot in Protein_list:
            if prot not in molecules:
                print '{} was not found in the .top file'.format(prot)
                return None
            molecules.remove(prot)
        molecules.insert(0, 'Protein')

    #create the reordered file, using grep and the 'molecules' list, to get the order the same as the top file
    print 'reordering genconf output file...'
    #script assumes that the atom information starts on the third line - This is a safe assumption for gro files
    os.system('head -2 {} > {}'.format(genconf_filename, reordered_filename))

    for mol in molecules:
	if mol == 'NA+':
	    mol = 'NA\+'
	elif  mol == 'CL-':
	    mol = 'CL\-'
		
	if re.match(r'Protein', mol, flags = re.IGNORECASE):
	    if len(Protein_list) > 1:
		os.system("egrep '^[0-9 ][0-9 ][0-9 ][0-9 ][0-9] *(ALA|CYS|ASP|GLU|PHE|GLY|HIS|ILE|LYS|LEU|MET|ASN|PRO|GLN|ARG|SER|THR|VAL|TRP|TYR)' {} >> temp_protein_file.gro".format(genconf_filename)) # more precise but probs not necessary
		# write protein to temporary file so that if there is >1 prot type, they can be ordered correctly
		# ie. all proteins of the same type grouped
		protfile = open('temp_protein_file.gro')
		protfile2 = open('temp_protein_file2.gro','w')
		all_lines = protfile.readlines()
		n_previous_prot_atoms = 0
		for prot in Protein_list:
		    prot_lines = []
		    for i in range(n_repeats):
			start = total_prot_atoms_1repeat*i + n_previous_prot_atoms
			end  = start + natoms_in_prot1repeat[prot]
			prot_lines += all_lines[ start : end ]
		    n_previous_prot_atoms += natoms_in_prot1repeat[prot]
		    for line in prot_lines:
			protfile2.write(line)
		protfile.close()
		protfile2.close()
		os.system("cat temp_protein_file2.gro >> {}".format(reordered_filename))
		os.remove('temp_protein_file.gro')
		os.remove('temp_protein_file2.gro')
	    else:
		os.system("egrep '^[0-9 ][0-9 ][0-9 ][0-9 ][0-9] *(ALA|CYS|ASP|GLU|PHE|GLY|HIS|ILE|LYS|LEU|MET|ASN|PRO|GLN|ARG|SER|THR|VAL|TRP|TYR)' {} >> {}".format(genconf_filename, reordered_filename)) # more precise but probs not necessary
		#os.system("egrep '(ALA|CYS|ASP|GLU|PHE|GLY|HIS|ILE|LYS|LEU|MET|ASN|PRO|GLN|ARG|SER|THR|VAL|TRP|TYR)' {} >> {}".format(genconf_filename, reordered_filename))
	else:
	    #os.system("egrep '^[0-9 ][0-9 ][0-9 ][0-9 ][0-9] *{}' {} >> {}".format(mol, genconf_filename, reordered_filename))
	    os.system("egrep '{}' {} >> {}".format(mol, genconf_filename, reordered_filename))

    os.system('tail -1 {} >> {}'.format(genconf_filename, reordered_filename))

    # renumber the reordered file so that atom and residue numbers are numbered correctly
    # this uses the script 'renumber_gro_file'
    print 're-numbering the newly-created reordered.gro file...'
    renumber_atomnumbers_only(reordered_filename)
    
    store_old_copy_gromacs_style(genconf_filename)
    os.remove(reordered_filename)
    os.rename('new_'+reordered_filename, genconf_filename) # renumbered residues outputs filename as 'new_'+filename

    #check the number of residues expected from the toplogy file
    #protein_itp = sys.argv[4]
    #add this later - would be easy to just add up residue numbers for all lipids and solvent, but protein has more than one residue per molecule, and  making sense of that requires knowledge of n_residues per molecule from the itp file

def reorder_genconf(top_filename, n_repeats, genconf_file, force_field='martini'):
    molecules, nmols = make_new_top_file(top_filename, n_repeats, genconf_file)
    print molecules, nmols
    reorder_genconf_file(genconf_file, molecules, nmols, n_repeats, force_field)
    #renumber residues, starting from 1
    ## not yet tested
    os.system('gmx_sse editconf -f {} -resnr 1 -o {}'.format(genconf_file, genconf_file))

    

if __name__ == "__main__":
    if len(sys.argv) == 4:
        reorder_genconf(sys.argv[1], int(sys.argv[2]), sys.argv[3])
    #    make_new_top_file(sys.argv[1], sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 5:
	reorder_genconf(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4])
    else:
        print 'Arguments: top_filename n_repeats genconf_filename OPTIONAL: force_field=<martini - default - or bondini>'
        print 'NB. n_repeats = n_Xrepeats * n_Yrepreats * n_Zrepeats'
        print 'eg. if -box input to genconf was 3 3 1, n_repeats = 9 (3*3*1)'
    print '''***Warning***: assumptions made: 1. genconf_filename is a gro file
                    2. protein molecules are all listed together at the start of the top file 
		    3. NA+ and CL- are the only names expected with '+' or '-' characters  - edit lines c. 80-85 to change this'''        

