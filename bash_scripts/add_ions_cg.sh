MDPFILE=`ls *.mdp`
TOPFILE=`ls *.top`
GROFILE=`ls *.gro`
GROFILE_BN=`basename --suffix=.gro -- *.gro`
XTCFILE=`ls *.xtc`
NDXFILE=`ls *.ndx`
CPTFILE=`ls *.cpt`

REQ_FILES_DIR='/biggin/b123/sedm5059/SCRIPTS/files_req'
cp "$REQ_FILES_DIR"/water-box-CG-303K-1bar.gro .

# solvate system with water

gmx_sse solvate -cp $GROFILE -cs water-box-CG-303K-1bar.gro -radius 0.21 -o "$GROFILE_BN"_solvated.gro

rm water-box-CG-303K-1bar.gro
echo "starting python!!!!!"
python ~/SCRIPTS/python_scripts/rem_water_cg.py



Wnum=`grep "W " "$GROFILE_BN"_solvated.gro | wc -l`

echo "W $Wnum" >> $TOPFILE 

cp $TOPFILE topol_solvated.top
TOPFILE_SOLV="topol_solvated.top"

# add ions 

gmx_sse grompp -f $MDPFILE -p $TOPFILE_SOLV -c "$GROFILE_BN"_solvated.gro -o minimisation_solvated_system.tpr
gmx_sse genion -s minimisation_solvated_system.tpr -p $TOPFILE_SOLV -conc 0.15 -nname CL -pname NA -neutral -o "$GROFILE_BN"_solvated.gro

# generate an energy minimisation .tpr file 

gmx_sse grompp -f $MDPFILE -c "$GROFILE_BN"_solvated.gro -p $TOPFILE_SOLV -o minim_"$GROFILE_BN"_solvated.tpr

mv $TOPFILE_SOLV "$GROFILE_BN"_solvated.top
rm $TOPFILE



# here you make a tpr file that is ready to be exectuted.
