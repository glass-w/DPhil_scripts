#!/bin/bash --login

#PBS -N benchmarkNCORES
#PBS -l select=NNODES
#PBS -l walltime=00:15:0
#PBS -A e280-Biggin

cd $PBS_O_WORKDIR

module load gromacs/5.0.0

aprun -n NCORES mdrun_mpi -s pr_tictac_2msr0.tpr -stepout 25000 -deffnm NCORESn -rdd 1.9 -maxh 0.25 -noconfout -resethway -append -pin auto -rcon 0.1
 
