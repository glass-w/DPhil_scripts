#!/bin/bash

# A script to remove artifacts from pbc in simulations

TPRFILE=`ls *.tpr`
GROFILE=`ls *.gro`
XTCFILE=`ls *.xtc`
NDXFILE=`ls *.ndx`
XTC_OUT=`basename --suffix=.xtc -- $XTCFILE`
GRO_OUT=`basename --suffix=.gro -- $GROFILE`

gmx_sse make_ndx -f $GROFILE -o index.ndx

gmx_sse trjconv -f $XTCFILE -s $TPRFILE -n $NDXFILE -pbc nojump -center -o "$XTC_OUT".xtc 

gmx_sse trjconv -f $GROFILE -s $TPRFILE -n $NDXFILE -pbc nojump -center -o "$GRO_OUT".gro


