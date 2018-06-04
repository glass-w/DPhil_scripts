#!/bin/bash

#echo 0 | trjconv_sse -f $1 -s $2 -pbc mol -ur compact -o $(basename $1 '.xtc')_noPBC.xtc

#<<COMMENT1

TPR=$1
XTC=$2
SKIP=$3

if [ "$SKIP" -ge "1" ] 
then

trjconv_avx -f $XTC -skip $SKIP -o $(basename $XTC '.xtc')skip$SKIP.xtc
#cp $XTC $(basename $XTC '.xtc')skip$SKIP.xtc

echo 'system' | trjconv_avx -f $(basename $XTC '.xtc')skip$SKIP.xtc -pbc whole -s $TPR -o $(basename $XTC '.xtc')_noPBCWhole_skip$SKIP.xtc

echo 'system' | trjconv_avx -f $(basename $XTC '.xtc')_noPBCWhole_skip$SKIP.xtc -o $(basename $XTC '.xtc')_noPBCWhole_noJump_skip$SKIP.xtc -pbc nojump -s $TPR

echo 'protein' 'system' | trjconv_avx -f $(basename $XTC '.xtc')_noPBCWhole_noJump_skip$SKIP.xtc -o $(basename $XTC '.xtc')_noPBCWhole_noJump_Center_skip$SKIP.xtc -pbc atom -center -ur compact -s $TPR -n index.ndx

echo 'protein' 'system' | trjconv_avx -f $(basename $XTC '.xtc')_noPBCWhole_noJump_Center_skip$SKIP.xtc -o $(basename $XTC '.xtc')_noPBCComplete_skip$SKIP.xtc -fit progressive -s $TPR -n index.ndx

#rm Prod_output_noPBCWhole.xtc Prod_output_noPBCWhole_noJump.xtc Prod_output_noPBCWhole_noJump_Center.xtc

#COMMENT1

else

echo 'system' | trjconv_avx -f $XTC -pbc whole -s $TPR -o $(basename $XTC '.xtc')_noPBCWhole.xtc

echo 'system' | trjconv_avx -f $(basename $XTC '.xtc')_noPBCWhole.xtc -o $(basename $XTC '.xtc')_noPBCWhole_noJump.xtc -pbc nojump -s $TPR 

echo 'protein' 'system' | trjconv_avx -f $(basename $XTC '.xtc')_noPBCWhole_noJump.xtc -o $(basename $XTC '.xtc')_noPBCWhole_noJump_Center.xtc -pbc atom -center -ur compact -s $TPR 

echo 'protein' 'system' | trjconv_avx -f $(basename $XTC '.xtc')_noPBCWhole_noJump_Center.xtc -o $(basename $XTC '.xtc')_noPBCComplete.xtc -fit progressive -s $TPR

fi
