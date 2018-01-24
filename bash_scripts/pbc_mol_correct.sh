#!/bin/bash

# A script to remove artifacts from pbc in simulations

TPRFILE=`ls *.tpr`
GROFILE=`ls *.gro`
XTCFILE=`ls *.xtc`
NDXFILE=`ls *.ndx`

gmx_sse make_ndx -f $GROFILE -o index.ndx

trjconv_avx -f $XTCFILE -skip 10 -o $(basename $XTCFILE '.xtc')_skip10.xtc

SKIP_XTCFILE=$(basename $XTCFILE '.xtc')_skip10.xtc

gmx_sse trjconv -f $SKIP_XTCFILE -s $TPRFILE -n $NDXFILE -pbc mol -o "$(basename $SKIP_XTCFILE '.xtc')"_pbcmol_corrected_lipid_and_prot.xtc

gmx_sse trjconv -f $GROFILE -s $TPRFILE -n $NDXFILE -pbc mol -o "$(basename $GROFILE '.gro')"_pbcmol_corrected_lipid_and_prot.gro
