ORIG=`pwd`

i=0
for d in ./*/ ; do 

	cd "$d"

	CURRDIR=`pwd`
	echo $CURRDIR
 
	cd "$CURRDIR/input" 
		
	mv confout.gro confout_from_short_equil.gro
	mv state.cpt state_from_short_equil.cpt

	gmx_sse grompp -f long_npt.mdp -c confout_from_short_equil.gro -t state_from_short_equil.cpt -p topol.top -o npt_R3r"$i"_long.tpr
	
	cp *.tpr $CURRDIR/output 

        let i++	
	
	cd $ORIG

done
