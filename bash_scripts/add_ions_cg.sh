#!/usr/bin/env bash
MDPFILE=`ls *.mdp`
TOPFILE=`ls *.top`
GROFILE=`ls *.gro`
GROFILE_BN=`basename --suffix=.gro -- *.gro`

REQ_FILES_DIR='/biggin/b123/sedm5059/SCRIPTS/files_req'
cp "$REQ_FILES_DIR"/water-box-CG-303K-1bar.gro .

# solvate system with water

gmx_sse solvate -cp $GROFILE -cs water-box-CG-303K-1bar.gro -radius 0.21 -o "$GROFILE_BN"_solvated.gro

# check if we have a membrane in our system (assuming AT POPC membrane), if so then remove waters inside membrane.
in_mrbane_check=`grep "PO4" $GROFILE | wc -l`

if [ $in_mrbane_check != 0 ]
then

rm water-box-CG-303K-1bar.gro
echo ""
echo "*** starting python script ***"
python ~/SCRIPTS/python_scripts/rem_water_cg.py "$GROFILE_BN"_solvated.gro
echo ""
echo "*** python script complete ***"
echo ""

fi

Wnum=`grep "W        W" "$GROFILE_BN"_solvated.gro | wc -l`

echo "W $Wnum" >> $TOPFILE

cp $TOPFILE topol_solvated.top
TOPFILE_SOLV="topol_solvated.top"

# add ions 

gmx_sse grompp -f $MDPFILE -p $TOPFILE_SOLV -c "$GROFILE_BN"_solvated.gro -o minimisation_solvated_system.tpr
gmx_sse genion -s minimisation_solvated_system.tpr -p $TOPFILE_SOLV -conc 0.15 -nname CL -pname NA -neutral -o "$GROFILE_BN"_solvated.gro

rm minimisation_solvated_system.tpr

# generate an energy minimisation .tpr file 

gmx_sse grompp -f $MDPFILE -c "$GROFILE_BN"_solvated.gro -p $TOPFILE_SOLV -o minim_"$GROFILE_BN"_solvated.tpr

mv $TOPFILE_SOLV "$GROFILE_BN"_solvated.top
#rm $TOPFILE
