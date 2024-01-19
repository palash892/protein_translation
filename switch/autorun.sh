gmx="/apps/codes/gromacs/2019/bin//gmx"
python='/home/palashb/miniconda3/bin/python'
nsims1=4
nsims2=125

nsteps=80000

arg1=$1


if [ -d log ]; then
	rm -rf log
fi

if [ -f cat.txt ]; then
    rm -rf cat.txt
fi


touch cat.txt
mkdir log



traj=""
$gmx grompp -f prod_sd_gpu.mdp -c equi$arg1".gro" -t equi$arg1".cpt" -o sd1.tpr -p topol.top >& log/log_grompp_md
$gmx mdrun -v -deffnm sd1 -nb gpu -nt 2 -nsteps $nsteps >& log/log_mdrun_md &
wait
echo 'c' >> cat.txt
traj="$traj sd1.xtc"

for((i=2;i<$nsims1;i++))
do
	$python switch.py $((i-1))
	echo python create the file
	

	$gmx grompp -f nrg_min_gpu_part.mdp -c next"$i".gro -o em"$i".tpr -p topol.top >& log/log_grompp_em"$i"
	$gmx mdrun -v -deffnm em"$i" -nb gpu -nt 2 >& log/log_mdrun_em"$i" &
	wait
	
	
	$gmx grompp -f prod_sd_gpu_1.mdp -c em"$i".gro -o sd"$i".tpr -p topol.top >& log/log_grompp_md"$i"
	$gmx mdrun -v -deffnm sd"$i" -nb gpu -nt 2 -nsteps $nsteps >& log/log_mdrun_md"$i" &
	wait

	echo 'c' >> cat.txt
	traj="$traj sd$((i)).xtc"

	rm *#

done

for((i=$nsims1;i<=$nsims2;i++))
do
	$python switch.py $((i-1))
	echo python create the file

	$gmx grompp -f nrg_min_gpu_part.mdp -c next"$i".gro -o em"$i".tpr -p topol.top  >& log/log_grompp_em"$i"
	$gmx mdrun -v -deffnm em"$i" -nb gpu -nt 2 >& log/log_mdrun_em"$i" &
	wait
	
	
	$gmx grompp -f prod_sd_gpu_2.mdp -c em"$i".gro -o sd"$i".tpr -p topol.top  >& log/log_grompp_md"$i"
	$gmx mdrun -v -deffnm sd"$i" -nb gpu -nt 2 -nsteps $nsteps >& log/log_mdrun_md"$i" &
	wait

	echo 'c' >> cat.txt
	traj="$traj sd$((i)).xtc"

	rm *#

done

$gmx trjcat -f $traj -o combined.xtc -settime yes < cat.txt

rm -rf *#
rm -rf em*

