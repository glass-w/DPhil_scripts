MDPFILE=`ls *.mdp`
TOPFILE=`ls *.top`
GROFILE=`ls *.gro`
GROFILE_BN=`basename --suffix=.gro -- *.gro`
XTCFILE=`ls *.xtc`
NDXFILE=`ls *.ndx`
CPTFILE=`ls *.cpt`

# solvate system with water

gmx_sse solvate -cp $GROFILE -o "$GROFILE_BN"_solvated.gro -p $TOPFILE

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
