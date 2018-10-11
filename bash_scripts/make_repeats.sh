# A script to generate repeats, in the mdp file the gen-seed is random each time so this gives 4 distinct repeat simualtion setups

MDPFILE=`ls *.mdp`
TPRFILE_BN=`basename --suffix=.mdp -- *.mdp` # name the tpr file in line with the mdp file name 
TPRFILE_NAME="$TPRFILE_BN"

TOPFILE=`ls *.top`
GROFILE=`ls *.gro`
XTCFILE=`ls *.xtc`
NDXFILE=`ls *.ndx`
CPTFILE=`ls *.cpt`

count_0=`ls -1 *.ndx 2>/dev/null | wc -l`
if [ $count_0 = 0 ]
then
gmx_sse make_ndx -f $GROFILE -o index.ndx
fi

count=`ls -1 *.cpt 2>/dev/null | wc -l`
if [ $count != 0 ]
then 

for i in `seq 0 3` ; do

        gmx_sse grompp -f $MDPFILE -p $TOPFILE -c $GROFILE -n index.ndx -t $CPTFILE -o "$TPRFILE_NAME"_r"$i".tpr -maxwarn 2 

done

else

for i in `seq 0 3` ; do

	gmx_sse grompp -f $MDPFILE -p $TOPFILE -c $GROFILE -n index.ndx -o "$TPRFILE_NAME"_r"$i".tpr -maxwarn 2 
done

fi

#############################
#                           #
#   if using GROMACS 2018   #	
#                           # 
#############################


if [ $1 = "2018" ]; then

count_0=`ls -1 *.ndx 2>/dev/null | wc -l`
if [ $count_0 = 0 ]
then
gmx make_ndx -f $GROFILE -o index.ndx
fi

count=`ls -1 *.cpt 2>/dev/null | wc -l`
if [ $count != 0 ]
then

#for i in `seq 0 3` ; do
for i in `seq 1 3` ; do
        gmx grompp -f $MDPFILE -p $TOPFILE -c $GROFILE -n index.ndx -t $CPTFILE -o "$TPRFILE_NAME"_r"$i".tpr -maxwarn 2 -r $GROFILE
	rm mdout.mdp
done

else

#for i in `seq 0 3` ; do
for i in `seq 1 3` ; do
        gmx grompp -f $MDPFILE -p $TOPFILE -c $GROFILE -n index.ndx -o "$TPRFILE_NAME"_r"$i".tpr -maxwarn 2 -r $GROFILE
	rm mdout.mdp
done

fi

fi
