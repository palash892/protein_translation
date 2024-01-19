#!/bin/bash
#PBS -l nodes=jm4:ppn=16
#PBS -N sw_80
#PBS -q batch1


cd $PBS_O_WORKDIR

module purge
module load libraries/fftw-3.3.8
module load compilers/nvcc/cuda-10
module load codes/gromacs/2019_jm9

export PYTHONPATH=$PYTHONPATH:/home/palashb/python_lib



for f in {17..32}
do
	mkdir seed_$f
	cd seed_$f
	cp -rv ../{*.top,*.itp,*.mdp} .
	cp -rv ../autorun.sh .
	cp -rv ../new_switch_no_overlap.py .
	cp -rv ../../../equilibration_dt_1e-4/seed_$f/{equi$f".gro",equi$f".cpt"} .
	bash autorun.sh $f >& job.log &
	cd ../
done
wait

