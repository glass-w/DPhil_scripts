################################################################
#	A script to extract certian residues from Nav1.5       #
#	Used to extract the S1-4 D1 region of Nav1.5	       #
#                       William Glass Jul 17 		       #
################################################################

#### NOTE ####

# The indices used for make_ndx and trjconv will differ b/w systems
# remove echo in this case and use this script manually

#############

ORIG=`pwd`

GROFILE=`ls *.gro`
GROFILENAME=`ls *.gro | cut -d "." -f 1`
XTCFILE=`ls *.xtc`
TPRFILE=`ls *.tpr`

#REPEAT=${PWD##*/} # assumes the directory you're in is named as the repeat e.g r0, r1, r2 etc

#ROTAMER=`echo "$ORIG" | rev | cut -d/ -f2 | rev` # takes the parent directory name of the current directory

##### MAIN PROG #####

{ echo "ri 90-190"; echo "q"; } | gmx_sse make_ndx -f $GROFILE -o index.ndx

NDXFILE=`ls *.ndx`

{ echo "24"; echo "24"; } | gmx_sse trjconv -f $GROFILE -s $TPRFILE -n $NDXFILE -pbc mol -center -o "$GROFILENAME"_resid_90_190.gro
{ echo "24"; echo "24"; } | gmx_sse trjconv -f $XTCFILE -s $TPRFILE -n $NDXFILE -pbc mol -center -o "$GROFILENAME"_resid_90_190.xtc

