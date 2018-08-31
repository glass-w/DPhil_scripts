#!/usr/bin/env bash

################################################################################################################################
#															       # 
# Author: W. Glass 2018													       #
# Description: 														       #
#															       #
#	This script takes a .gro file without any water (this can be with or without a membrane) and generates                 # 
#	a solvated system with ions along with a .tpr file (assumed to be for energy minimisation).                            #
#															       #
#	1) gather info for input files.											       #
#	2) solvate the system and check if there is a membrane present, if not then skip the water removal python script.      #
#	3) if there is a membrane, start the python script (py2.7).							       #
#	4) add the new number of water molecules to the supplied .top file.						       #
#	5) add ions and create a .tpr file ready for minimisation, this is assuming the .mdp file in the current directory is  #
#	an energy minimisation script.											       #
#															       #
################################################################################################################################

MDPFILE=`ls *.mdp`
TOPFILE=`ls *.top`
GROFILE=`ls *.gro`
GROFILE_BN=`basename --suffix=.gro -- *.gro`
XTCFILE=`ls *.xtc`
NDXFILE=`ls *.ndx`
CPTFILE=`ls *.cpt`

# Path to python water removal script
PY_WATER_REM_PATH=~/SCRIPTS/python_scripts

# solvate system with water, don't read in topology file. Will add new waters after water removal script.

gmx_sse solvate -cp $GROFILE -o "$GROFILE_BN"_solvated.gro # -p $TOPFILE

# check if we have a membrane in our system (assuming AT POPC membrane), if so then remove waters inside membrane.
in_mrbane_check=`grep "POPC" $GROFILE | wc -l`

if [ $in_mrbane_check != 0 ]
then

echo ""
echo "*** starting python script ***"
echo ""
echo "*** removing waters in bilayer... ***"
python "$PY_WATER_REM_PATH"/rem_water_at_py36.py "$GROFILE_BN"_solvated.gro
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

gmx_sse grompp -f $MDPFILE -p $TOPFILE_SOLV -c "$GROFILE_BN"_solvated.gro -o temp.tpr
echo 'SOL' | gmx_sse genion -s temp.tpr -p $TOPFILE_SOLV -conc 0.15 -nname CL -pname NA -neutral -o "$GROFILE_BN"_solvated.gro

# generate an energy minimisation .tpr file 

gmx_sse grompp -f $MDPFILE -c "$GROFILE_BN"_solvated.gro -p $TOPFILE_SOLV -o minim_"$GROFILE_BN"_solvated.tpr

mv $TOPFILE_SOLV "$GROFILE_BN"_solvated.top

# Clean Up
rm $TOPFILE
rm temp.tpr # this is only needed to add the ions.

# You now have a .tpr file that is ready to be exectuted for energy minimisation.
