################################
#			       #
# A script to make a .tpr file #
#			       #
################################

# This script first checks if certain files are present, and if there are .ndx and .cpt files
# If there are no .ndx files, it creates one.
# If there is a .cpt file it uses this in grompp

MDPFILE=`ls *.mdp`
TPRFILE_BN=`basename --suffix=.mdp -- *.mdp` # name the tpr file in line with the mdp file name 
TPRFILE_NAME="$TPRFILE_BN"

TOPFILE=`ls *.top`
GROFILE=`ls *.gro`

xtc_count=`ls -1 *.xtc 2>/dev/null | wc -l`
if [ $xtc_count != 0 ]
then

XTCFILE=`ls *.xtc`

fi

ndx_count=`ls -1 *.ndx 2>/dev/null | wc -l`
if [ $ndx_count = 0 ]
then

gmx_sse make_ndx -f $GROFILE -o index.ndx

NDXFILE=`ls *.ndx`

fi

cpt_count=`ls -1 *.cpt 2>/dev/null | wc -l`
if [ $cpt_count != 0 ]
then

CPTFILE=`ls *.cpt`

gmx_sse grompp -f $MDPFILE -p $TOPFILE -c $GROFILE -n index.ndx -t $CPTFILE -o "$TPRFILE_NAME".tpr -maxwarn 2

elif [ $cpt_count = 0 ]
then

gmx_sse grompp -f $MDPFILE -p $TOPFILE -c $GROFILE -n index.ndx -o "$TPRFILE_NAME".tpr -maxwarn 2

fi
