# A script to generate repeats, in the mdp file the gen-seed is random each time so this gives 4 distinct repeat simualtion setups

for i in {0..3} ; do

	gmx_sse grompp -f run.mdp -p topol.top -c confout_R2_system.gro -o npt_short_R2_r$i.tpr -maxwarn 1

done

