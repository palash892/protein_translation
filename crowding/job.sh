#!/bin/bash
#PBS -l nodes=jm5:ppn=16
#PBS -N pf_0.57
#PBS -q batch1


cd $PBS_O_WORKDIR
module load libraries/fftw-3.3.8
module load compilers/nvcc/cuda-10
module load codes/gromacs/2019_jm9


for f in {1..16}
do
	mkdir seed_$f
	cd seed_$f
	cp -rv ../{*.top,*.itp,*.mdp} .
	cp -rv ../../equilibration_dt_1e-4/seed_$f/{equi$f".gro",equi$f".cpt"} .
	mkdir log
	gmx grompp -f prod_sd_gpu_1.mdp -c equi$f".gro" -t equi$f".cpt" -p topol.top -o run1.tpr >& log/grompp_1st.log &
	cd ../
done
wait

for f in {1..16}
do
	cd seed_$f
	gmx mdrun -v -deffnm run1 -nsteps 50000 -nb gpu -nt 2 >& log/mdrun1.log &
	cd ../
done
wait

for f in {1..16}
do
	cd seed_$f
	gmx grompp -f prod_sd_gpu_2.mdp -c run1.gro -t run1.cpt -p topol.top -o run2.tpr >& log/grompp_2nd.log &
	cd ../
done
wait


for f in {1..16}
do
	cd seed_$f
	gmx mdrun -v -deffnm run2 -nsteps 10000000 -nb gpu -nt 2 >& log/mdrun2.log &
	cd ../
done


wait 
