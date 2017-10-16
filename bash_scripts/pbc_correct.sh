#!/bin/bash

# A script to remove artifacts from pbc in simulations

TPRFILE=`ls *.tpr`
GROFILE=`ls *.gro`
XTCFILE=`ls *.xtc`
NDXFILE=`ls *.ndx`

gmx_sse make_ndx -f $GROFILE -o index.ndx

gmx_sse trjconv -f $XTCFILE -s $TPRFILE -n $NDXFILE -pbc mol -center -o traj_pbc_corrected.xtc

gmx_sse trjconv -f $GROFILE -s $TPRFILE -n $NDXFILE -pbc mol -center -o confout_pbc_corrected.gro


