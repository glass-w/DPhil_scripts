#!/bin/bash

''' NOTE: this script does change the input .gro files in such a way that vmd will not recognise them,
 save a separate copy somewhere - this is to be fixed

 First input is the protein .gro file, second input is the membrane .gro file

 '''

box_increase=7

#count atoms

num1=$(sed -n '2p' $1)
num2=$(sed -n '2p' $2)
num=$((num1 + num2))

#get initial dimensions

x_y_z=$(tail -1 $2)

#split dimensions and add more to z

x=$(echo "$x_y_z" | awk '{print $1}')
y=$(echo "$x_y_z" | awk '{print $2}') 
old_z=$(echo "$x_y_z" | awk '{print $3}') 
z=$(echo $old_z + $box_increase | bc)

#apply new dimensions to protein and lipid files individually

#echo 'protein' | gmx_sse editconf -f $1 -o newbox.gro -c -rotate 0 -15 0 -box $x $y $z
echo 'protein' | gmx_sse editconf -f $1 -o newbox.gro -c -rotate 0 -15 -15 -box $x $y $z
gmx_sse editconf -f newbox.gro -translate 0 0 0 -o translated_newbox.gro

gmx_sse editconf -f $2 -o $2 -c -box $x $y $z

# Generate strong position restraints for inflateGRO later
echo 'protein' | gmx_sse genrestr -f newbox.gro -o strong_posre.itp -fc 100000 100000 100000

#delete box dimensions, title and number of atoms in protein file

sed -i '1,2d;$d' translated_newbox.gro #remove box dimensions from protein file
sed -i '1,2d' $2  		#remove title and num atoms from lipid file
sed -e "\$a$x $y $z" $2       #add new dimensions to end of lipid file

cat translated_newbox.gro $2 > combined.gro

#add title and number of atoms to top of .gro file

sed -i '1i '$num''  combined.gro
sed -i '1i Title' combined.gro

#renumber
gmx_sse genconf -f combined.gro -renumber -o combined.gro
