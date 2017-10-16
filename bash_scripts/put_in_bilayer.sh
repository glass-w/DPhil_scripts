#!/bin/bash

''' NOTE: this script does change the input .gro files in such a way that vmd will not recognise them, save a sepreate copy somewhere - this is to be fixed'''

box_increase=5

#count atoms

num1=$(sed -n '2p' beta3_tm_ext8.gro)
num2=$(sed -n '2p' popc_cutdown.gro)
num=$((num1 + num2))

#get initial dimensions

x_y_z=$(tail -1 popc_cutdown.gro)

#split dimensions and add more to z

x=$(echo "$x_y_z" | awk '{print $1}')
y=$(echo "$x_y_z" | awk '{print $2}') 
old_z=$(echo "$x_y_z" | awk '{print $3}') 
z=$(echo $old_z + $box_increase | bc)

#apply new dimensions to protein and lipid files individually

gmx_sse editconf -f beta3_tm_ext8.gro -o beta3_tm_ext8_newbox.gro -c -box $x $y $z -princ -rotate 0 90 0 
gmx_sse editconf -f popc_cutdown.gro -o popc_cutdown.gro -c -box $x $y $z

# Generate stronf position restraints for inflateGRO later
gmx_sse genrestr -f beta3_tm_ext8_newbox.gro -o strong_posre.itp -fc 100000 100000 100000

#delete box dimensions, title and number of atoms in protein file

sed -i '1,2d;$d' beta3_tm_ext8_newbox.gro #remove box dimensions from protein file
sed -i '1,2d' popc_cutdown.gro  		#remove title and num atoms from lipid file
sed -e "\$a$x $y $z" popc_cutdown.gro       #add new dimensions to end of lipid file

cat beta3_tm_ext8_newbox.gro popc_cutdown.gro > combined.gro

#add title and number of atoms to top of .gro file

sed -i '1i '$num''  combined.gro
sed -i '1i Title' combined.gro

#renumber 
gmx_sse genconf -f combined.gro -renumber -o combined.gro
