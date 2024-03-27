# Protein Molecular Dynamics Simulations

This repository contains scripts and parameter files for performing molecular dynamics simulations of proteins using GROMACS. It uses amber03 force field (https://pubmed.ncbi.nlm.nih.gov/14531054/) to simulate the dynamics of a protein in a 80 Å x 80 Å x 80 Å cube, dissolved in 18,000 TIP3 waters with 0.1M NaCl ions. After performing energy minimization with steepest descent integrator, the system is gradually heated to 300 K. The leap frog integrator is then used to produce a 500 ns MD trajectory with a time-step of 2.5 fs under NPT condition.



## Prerequisites

- GROMACS 2020 or later
- A PDB file of the protein structure (protein.pdb)

  
## Files

- **protein_md.sh**: Bash script for running the entire molecular dynamics simulation workflow.
- **parameters/**: Folder containing the GROMACS parameter files
    - **ion.mdp**: Parameters for adding ions to the system.
    - **md.mdp**: Parameters for the molecular dynamics run.
    - **minim.mdp**: Parameters for energy minimization.
    - **npt.mdp**: Parameters for NPT equilibration.
    - **nvt.mdp**: Parameters for NVT equilibration.
 
## Usage

1. Place your **protein.pdb** file in the same directory as the **protein_md.sh** script.
2. Make the script executable:    `chmod +x protein_md.sh`
3. Run the script: `./protein_md.sh`

## Notes

- The script assumes that GROMACS is installed and the GMXRC script is located at `/usr/local/gromacs/bin/GMXRC`. Modify the source command in **protein_md.sh** if your GROMACS installation is located elsewhere.
- The salt concentration is set to $0.1 M$ in the script. Modify the -conc option in the genion command in protein_md.sh to change the concentration.

## Adapted from

- [MPIBPC CompBio1 Practical 4](http://www3.mpibpc.mpg.de/groups/de_groot/compbio1/p4/index.html#contents)
- [GROMACS MD Tutorial - Lysozyme in Water](http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html)


