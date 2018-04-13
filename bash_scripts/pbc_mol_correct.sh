#!/bin/bash

# A script to remove artifacts from pbc in simulations

TPRFILE=$1
GROFILE=$2
XTCFILE=$3
GROUP=$4

ndx_count=`ls -1 *.ndx 2>/dev/null | wc -l`
if [ $ndx_count = 0 ]
then

gmx_sse make_ndx -f $GROFILE -o index.ndx

NDXFILE=`ls *.ndx`

else

NDXFILE=`ls *.ndx`

fi

#skip a number of frames
#trjconv_avx -f $XTCFILE -skip 10 -o $(basename $XTCFILE '.xtc')_skip10.xtc

# write each frame at certain ns time, e.g. 10ns = 10000ps
trjconv_avx -f $XTCFILE -dt 10000 -o $(basename $XTCFILE '.xtc')_dt10.xtc

SKIP_XTCFILE=$(basename $XTCFILE '.xtc')_skip10.xtc

echo $GROUP | gmx_sse trjconv -f $SKIP_XTCFILE -s $TPRFILE -n $NDXFILE -pbc mol -o "$(basename $SKIP_XTCFILE '.xtc')"_pbcmol_corrected_lipid_and_prot.xtc

echo $GROUP | gmx_sse trjconv -f $XTCFILE -s $TPRFILE -n $NDXFILE -pbc mol -o pbcmol_corrected_lipid_and_prot.xtc

echo $GROUP | gmx_sse trjconv -f $GROFILE -s $TPRFILE -n $NDXFILE -pbc mol -o "$(basename $GROFILE '.gro')"_pbcmol_corrected_lipid_and_prot.gro
