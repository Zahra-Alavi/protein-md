#!/bin/bash

set -x

# Source the GROMACS environment
source /usr/local/gromacs/bin/GMXRC

# Run pdb2gmx
echo "1" | gmx pdb2gmx -f protein.pdb -o protein.gro -water spce

#Put the molecule in a box of wate
gmx editconf -f protein.gro -o protein_box.gro -c -d 2.0 -bt cubic
gmx solvate -cp protein.gro -cs spc216.gro -o protein_solv.gro -p topol.top

# Set salt concentration: 
gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o protein_ions.tpr
echo "13" | gmx genion -s protein_ions.tpr -o protein_ions.gro -p topol.top -pname NA -nname CL -conc 0.1

# Energy minimization and two-step equilibration: 
gmx grompp -f minim.mdp -c protein_ions.gro -p topol.top -o protein_em.tpr
gmx mdrun -v -deffnm protein_em
echo "10 0" | gmx energy -f protein_em.edr -o protein_potential.xvg

gmx grompp -f nvt.mdp -c protein_em.gro -r protein_em.gro -p topol.top -o protein_nvt.tpr
gmx mdrun -deffnm protein_nvt -nb gpu
echo "16 0" | gmx energy -f protein_nvt.edr -o protein_temperature.xvg

gmx grompp -f npt.mdp -c protein_nvt.gro -r protein_nvt.gro -t protein_nvt.cpt -p topol.top -o protein_npt.tpr
gmx mdrun -deffnm protein_npt -nb gpu
echo "18 0" | gmx energy -f protein_npt.edr -o protein_pressure.xvg
echo "24 0" | gmx energy -f protein_npt.edr -o protein_density.xvg

# MD run:
gmx grompp -f md.mdp -c protein_npt.gro -t protein_npt.cpt -p topol.top -o protein_md.tpr
gmx mdrun -deffnm protein_md -nb gpu


#Center the protein and account for periodicity: 
echo "1 0" | gmx trjconv -s protein_md.tpr -f protein_md.xtc -o protein_md_noPBC.xtc -pbc mol -center

# RMSD:
echo "4 4" | gmx rms -s protein_md.tpr -f protein_md_noPBC.xtc -o protein_rmsd.xvg -tu ns
echo "4 4" | gmx rms -s protein_em.tpr -f protein_md_noPBC.xtc -o protein_rmsd_xtal.xvg -tu ns

# Rg:
echo "1" | gmx gyrate -s protein_md.tpr -f protein_md_noPBC.xtc -o protein_gyrate.xvg

# RMSF per amino acid residue:
echo "1" | gmx rmsf -f protein_md_noPBC.xtc -s protein_md.tpr -o protein_rmsf.xvg -res

# Extract the protein trajectory:
echo "1" | gmx trjconv -s protein_npt.gro -f protein_md_noPBC.xtc -o protein_noPBC_nowater.xtc -e 500000

# Get the eigenvalues of the covariant matrix:
echo "1 1" | gmx covar -s protein_npt.gro -f protein_noPBC_nowater.xtc -o protein_eigenval.xvg -tu ns -ref -pbc

# Project the trajectory onto individual eigenvectors:
echo "1 1" | gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -proj protein_proj.xvg -first 1 -last 6 

# Plot eig1 vs eig2: 
echo "1 1" | gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -2d protein_2dproj.xvg -first 1 -last 2 

# Get the components of the first eigenvector: 
gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -comp protein_eigcomp.xvg -first 1 -last 6

# Get cosine content of PCs:
gmx analyze -f protein_proj.xvg -cc protein_cc.xvg -n 6

# Visualize the motion: 
echo "1 1" | gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -extr protein_extreme1.pdb -first 1 -last 1 -nframes 30

echo "1 1" | gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -filt protein_filter1.pdb -first 1 -last 1 -skip 100

echo "1 1" | gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -extr protein_extreme2.pdb -first 2 -last 2 -nframes 30

echo "1 1" | gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -filt protein_filter2.pdb -first 2 -last 2 -skip 100


# Repeat with C-alpha:

# Get the eigenvalues of the covariant matrix:
echo "3 3" | gmx covar -s protein_npt.gro -f protein_noPBC_nowater.xtc -o protein_eigenval_calpha.xvg -tu ns -ref -pbc

# Project the trajectory onto individual eigenvectors:
echo "3 3" | gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -proj protein_proj_calpha.xvg -first 1 -last 6 

# Plot eig1 vs eig2: 
echo "3 3" | gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -2d protein_2dproj_calpha.xvg -first 1 -last 2 

# Get the components of the first eigenvector: 
gmx anaeig -s protein_npt.gro -f protein_noPBC_nowater.xtc -comp protein_eigcomp_calpha.xvg -first 1 -last 6

# Get cosine content of PCs:
gmx analyze -f protein_proj_calpha.xvg -cc protein_cc_calpha.xvg -n 6

# Get the dihedrals
gmx chi -s protein_npt.gro -f protein_noPBC_nowater.xtc -corr protein_dihcorr.xvg

gmx chi -s protein_npt.gro -f protein_noPBC_nowater.xtc -rama





