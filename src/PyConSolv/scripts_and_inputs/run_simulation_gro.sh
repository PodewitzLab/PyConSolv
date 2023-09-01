gmx_mpi grompp -f simulation_gro.mdp -c npt.gro -t npt.cpt -p LIG_solv.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1 -nb gpu