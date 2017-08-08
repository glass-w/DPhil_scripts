'''Renumber moecules in the order expected for a gromacs file, 
ie. atom numbers increase by 1 with line number, residue numbers increase by 1 or 0 with line number,
 and when 99999 is reached in either case, numbering starts again at 0 (atom and residue numbers are of the format 5d).
This program (for now) assumes that the file has been edited so that molcules of the same type are all listed together,
 and just need renumbering.  Also that the 1st molecule in the file has residue # 1.

Note: gromacs atom line format is: "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" (inc. velocities).

THIS SCRIPT WILL NOT WORK IF, IN RE-ORDERING, TWO DIFFERENT RESIDUES HAVING THE SAME RESIDUE NUMBER, ARE CONSECUTIVE. 

This is a danger in files where the same residue has been added repeatedly.
For files collated using genconf, and then reordered using reorder_genconf, this is a danger where there are more than 99999 residues, and is a full-blown bug when a particular residue type only occurs once in the 'building block' file (ie. input file to genconf).

To correct this bug, the script should be written quite differently, ie. it should know, for each resname, how many atoms are expected, and renumber residues accordingly.

Atom and residue names will be read only if they have characters in the set [A-Z] or [a-z] or [0-9] or '_', '+', '-'.

This script will now also deal with two (or more) gro files that have been concatenated simply (ie. using cat).
( as long as this still means that molecules of the same type are listed together...)
In this case the gro file will have more than one header/#atoms/boxdimensions line;
the script will only re-write the top header line, the top #atoms line (replacing with the new #atoms),
and the last box dimensions line. 
(added 9/4/2014)
'''

import os
import re
import sys
import shutil

def count_atoms(gro_file):
    
    f = open(gro_file,'r')

    atomno = 0
    for line in f:
        if re.match(r'[\s\d]{5}[\w\s\+-]{10}[\s\d]{5}\s*[\d\.-]+\s*[\d\.-]+\s*[\d\.-]+',line):
            atomno += 1
    f.close()
    return atomno

def renumber_residues(filename):
    # for when residues have been 'jumbled around'
    
    total_atomno = count_atoms(filename)
    
    f = open(filename, 'r')
    g = open('new_'+filename, 'w')
    
    current_resno = 0
    prev_lineres = -1
    atomno = 0
    atom100Kcount = 0
    res100Kcount = 0
    for line in f:
        if re.match(r'[\s\d]{5}[\w\s\+-]{10}[\s\d]{5}\s*[-\.\d]+\s*[-\.\d]+\s*[-\.\d]+', line):
            line_res = int(line[:5])
            if line_res != prev_lineres:
                prev_lineres = line_res
                if current_resno == 99999: # reset counter if 99999 has been reached (this is the max atom/ res number in a gro file)
                    current_resno = 0
                    res100Kcount += 1
                else: # otherwise increase res number by one
                    current_resno += 1
            if atomno == 99999:  #as previous two comments
            	atomno = 0
            	atom100Kcount += 1
            	if atom100Kcount%10 == 0:
            	    sys.__stdout__.write('Reached the {} million-th atom\r'.format(atom100Kcount/10))
            else:
	            atomno += 1    

            g.write('%5d' % current_resno + line[5:15] + '%5d' % atomno + line[20:])

        elif re.match(r'\s*\d+\n', line): # ie. if this is the totalAtoms line
            if atomno + atom100Kcount == 0:  # totalAtoms line is written only when found at start of file
                                        # ie. if two gro files have been concatenated, the second totalAtoms line will be ignored
                                        # added ALD 9/4/2014
                g.write(' %d\n' % total_atomno)
        else: #ie. the line is the title string or box dimensions
            # title/box dimensions line only written if found at start/end of file
            # this deals with the case that two gro files have been concatenated
            if (atomno + atom100Kcount == 0) or (atomno + atom100Kcount*100000 == total_atomno): 
                g.write(line)
    print 'Reordered {} atoms (which should be the same as total_atomno count: {})'.format(atomno + atom100Kcount*100000, total_atomno)
    print 'The file contains {} residues - check that this is what you expect'.format(current_resno + res100Kcount*100000)

    f.close()
    g.close()


def renumber_atomnumbers_only(filename):
    # renumber atoms but not residues
    
    total_atomno = count_atoms(filename)
    
    f = open(filename, 'r')
    g = open('new_'+filename, 'w')
    
    atomno = 0
    atom100Kcount = 0
    for i,line in enumerate(f):
        if re.match(r'[\s\d]{5}[\w\s\+-]{10}[\s\d]{5}\s*[-\.\d]+\s*[-\.\d]+\s*[-\.\d]+', line):
            if atomno == 99999:  #as previous two comments
            	atomno = 0
            	atom100Kcount += 1
            	if atom100Kcount%10 == 0:
            	    sys.__stdout__.write('Reached the {} million-th atom\r'.format(atom100Kcount/10))
            else:
	            atomno += 1    

            g.write(line[:15] + '%5d' % atomno + line[20:])

        elif re.match(r'\s*\d+\n', line): # ie. if this is the totalAtoms line
            if atomno + atom100Kcount == 0:  # totalAtoms line is written only when found at start of file
                                        # ie. if two gro files have been concatenated, the second totalAtoms line will be ignored
                                        # added ALD 9/4/2014
                g.write(' %d\n' % total_atomno)
        else: #ie. the line is the title string or box dimensions
            # title/box dimensions line only written if found at start/end of file
            # this deals with the case that two gro files have been concatenated
            if (atomno + atom100Kcount == 0) or (atomno + atom100Kcount*100000 == total_atomno): 
                g.write(line)
    print 'Reordered {} atoms (which should be the same as total_atomno count: {})'.format(atomno + atom100Kcount*100000, total_atomno)

    f.close()
    g.close()
    os.remove(filename)
    shutil.copyfile('new_'+filename, filename)

def renumber_residues_onemoleculetype(filename, n_atoms_per_mol):
    # for when residues have been 'jumbled around'
    
    total_atomno = count_atoms(filename)
    
    f = open(filename, 'r')
    g = open('new_'+filename, 'w')
    
    current_resno = 1
    atomno = 0
    atom100Kcount = 0
    res100Kcount = 0
    n_atoms_count = 0
    for line in f:
        if re.match(r'[\s\d]{5}[\w\s\+-]{10}[\s\d]{5}\s*[-\.\d]+\s*[-\.\d]+\s*[-\.\d]+', line):
            line_res = int(line[:5])
            if n_atoms_count >= n_atoms_per_mol:
                if current_resno == 99999: # reset counter if 99999 has been reached (this is the max atom/ res number in a gro file)
                    current_resno = 0
                    res100Kcount += 1
                else: # otherwise increase res number by one
                    current_resno += 1
                n_atoms_count = 1
            else:
                n_atoms_count += 1 
            if atomno == 99999:  #as previous two comments
            	atomno = 0
            	atom100Kcount += 1
            	if atom100Kcount%10 == 0:
            	    sys.__stdout__.write('Reached the {} million-th atom\r'.format(atom100Kcount/10))
            else:
	            atomno += 1   

            g.write('%5d' % current_resno + line[5:15] + '%5d' % atomno + line[20:])

        elif re.match(r'\s*\d+\n', line): # ie. if this is the totalAtoms line
            if atomno + atom100Kcount == 0:  # totalAtoms line is written only when found at start of file
                                        # ie. if two gro files have been concatenated, the second totalAtoms line will be ignored
                                        # added ALD 9/4/2014
                g.write(' %d\n' % total_atomno)
        else: #ie. the line is the title string or box dimensions
            # title/box dimensions line only written if found at start/end of file
            # this deals with the case that two gro files have been concatenated
            if (atomno + atom100Kcount == 0) or (atomno + atom100Kcount*100000 == total_atomno): 
                g.write(line)
    print 'Reordered {} atoms (which should be the same as total_atomno count: {})'.format(atomno + atom100Kcount*100000, total_atomno)
    print 'The file contains {} residues - check that this is what you expect'.format(current_resno + res100Kcount*100000)

    f.close()
    g.close()



if __name__ == "__main__":
    if len(sys.argv) == 2:
        renumber_atomnumbers_only(sys.argv[1])
    elif len(sys.argv) == 3:
        renumber_residues_onemoleculetype(sys.argv[1], int(sys.argv[2]))
    else:
        print ' Usage: python renumber_gro_file.py <jumbled_file_name>'
