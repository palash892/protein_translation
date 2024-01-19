#!/bin/bash
#PBS -l nodes=jm4:ppn=32
#PBS -N at_0.2
#PBS -q batch1


cd $PBS_O_WORKDIR
module load libraries/fftw-3.3.8
module load compilers/nvcc/cuda-10

gmx='/home/palashb/programs/gromacs/gromacs_sph/bin/gmx'


for f in {1..16}
do
	mkdir seed_$f
    cd seed_$f
    mkdir log
    cp -rv ../../initial_configuration/pol$f".gro" .
    cp -rv ../{*.mdp,*.top,*.itp} .
    $gmx grompp -f prod_md.mdp -c pol$f".gro" -p topol.top -o prod.tpr -maxwarn 1 >& log/grompp_seed$f".log" &
	cd ../
done
wait

for f in {1..16}
do
    cd seed_$f
    $gmx mdrun -v -deffnm prod -nt 1 -nsteps 2000000 >& log/mdrun_seed$f".log" &
    rm *#
	cd ../
done
wait 
