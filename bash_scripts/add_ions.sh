#!/usr/bin/env bash
MDPFILE=`ls *.mdp`
TOPFILE=`ls *.top`
GROFILE=`ls *.gro`
GROFILE_BN=`basename --suffix=.gro -- *.gro`
XTCFILE=`ls *.xtc`
NDXFILE=`ls *.ndx`
CPTFILE=`ls *.cpt`

# solvate system with water

gmx_sse solvate -cp $GROFILE -o "$GROFILE_BN"_solvated.gro # -p $TOPFILE

#mv vdwradii.dat vdwradii.dat.used

# check if we have a membrane in our system, if so then remove waters inside membrane.
in_mrbane_check=`grep "POPC" $GROFILE | wc -l`

if [ $in_mrbane_check != 0 ]
then

echo ""
echo "*** starting python script ***"
echo ""
echo "*** removing waters in bilayer... ***"
python ~/SCRIPTS/python_scripts/rem_water_at.py "$GROFILE_BN"_solvated.gro
echo ""
echo "*** waters in bilayer removed! ***"
echo ""
echo "*** python script complete ***"
echo ""

fi

Wnum=`grep "SOL     OW" "$GROFILE_BN"_solvated.gro | wc -l`

echo "SOL $Wnum" >> $TOPFILE

cp $TOPFILE topol_solvated.top
TOPFILE_SOLV="topol_solvated.top"

# add ions 

gmx_sse grompp -f $MDPFILE -p $TOPFILE_SOLV -c "$GROFILE_BN"_solvated.gro -o minimisation_solvated_system.tpr
echo 'SOL' | gmx_sse genion -s minimisation_solvated_system.tpr -p $TOPFILE_SOLV -conc 0.15 -nname CL -pname NA -neutral -o "$GROFILE_BN"_solvated.gro

# generate an energy minimisation .tpr file 

gmx_sse grompp -f $MDPFILE -c "$GROFILE_BN"_solvated.gro -p $TOPFILE_SOLV -o minim_"$GROFILE_BN"_solvated.tpr

mv $TOPFILE_SOLV "$GROFILE_BN"_solvated.top
rm $TOPFILE



# here you make a tpr file that is ready to be exectuted.
