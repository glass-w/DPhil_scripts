#!/bin/bash

# A script to remove jumping from pbc in simulations

TPRFILE=$1
GROFILE=$2
XTCFILE=$3
NDXFILE=$4
XTC_OUT=`basename --suffix=.xtc -- $XTCFILE`
GRO_OUT=`basename --suffix=.gro -- $GROFILE`

#gmx_sse make_ndx -f $GROFILE -o index.ndx

gmx_sse trjconv -f $XTCFILE -s $TPRFILE -n $NDXFILE -pbc nojump -center -o "nojump_$XTC_OUT".xtc 

gmx_sse trjconv -f $GROFILE -s $TPRFILE -n $NDXFILE -pbc nojump -center -o "nojump_$GRO_OUT".gro


