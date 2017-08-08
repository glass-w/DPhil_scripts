ORIG=`pwd`
i=0

for d in ./*/ ; do

        cd "$d"
	cd analysis

        BASHSCRIPTS=/biggin/b123/sedm5059/SCRIPTS/bash_scripts
	PYTHONSCRIPTS=/biggin/b123/sedm5059/SCRIPTS/python_scripts
	VMDSCRIPTS=/biggin/b123/sedm5059/SCRIPTS/vmd_scripts

	vmd -dispdev text -e $VMDSCRIPTS/measure_tilt.tcl -args confout_pbc_corrected.gro traj_pbc_corrected.xtc 2 50000 "$i"_ext8_beta3_tilt_new

	((i++))

 #       bash "$BASHSCRIPTS"/pbc_correct.sh
 #	 python $PYTHONSCRIPTS/calc_rmsd.py
	
	#mkdir input output         
         
        cd $ORIG

done

