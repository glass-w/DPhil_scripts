#!/bin/bash

wd=`pwd`
counter=0

TOPFILE=`ls *.top`

for i in 500 250 125 62.5 31.25;
do

echo $TOPFILE

	sed 's/posre_Protein_chain_A/posre_Protein_chain_A_'$i'/' $TOPFILE > tmp_$i  
	sed 's/posre_Protein_chain_B/posre_Protein_chain_B_'$i'/' tmp_$i > tmp2_$i  
	sed 's/posre_Protein_chain_C/posre_Protein_chain_C_'$i'/' tmp2_$i > tmp3_$i
	cp tmp3_$i topol_$i.top
		
	### Run Gromacs ####

	# Make .tpr file
	

	# If in first loop and no .cpt file yet:
	
	if [ $i = 500 ]
	then
		gmx_sse grompp -f reduce_restraints_1ns.mdp -p topol_$i.top -c long_5ns_npt_r0.gro -n index.ndx -o restraints_$i.tpr -maxwarn 2		

		gmx_sse mdrun -v -deffnm restraints_$i -pin on -nt 10 

	# once the first run has finished:
	
	else  
	
		# take the previous .tpr file and extend it by 1ns	
		
		gmx_sse convert-tpr -s $tpr_to_extend -extend 1000 -o restraints_$i.tpr

		gmx_sse mdrun -v -deffnm restraints_$i -cpi $prev_cpt_file -pin on -nt 10 

	fi

	#####

	# assign the current tpr and cpt to be the 'old' ones for the next run
	
	prev_cpt_file=`ls restraints_$i.cpt`
	tpr_to_extend=`ls restraints_$i.tpr`
	
done                              
