#!/bin/bash --login

#PBS -N Kir3PM
#PBS -l select=10
#PBS -l walltime=24:00:0
#PBS -A e460

cd $PBS_O_WORKDIR

module load gromacs/4.6.5

aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*
aprun -n 240 mdrun_mpi  -s md_prod_3x3_10to50us.tpr -deffnm md_prod3x3 -rdd 1.9 -cpi md_prod3x3.cpt -noappend
rm step*
rm dd*

