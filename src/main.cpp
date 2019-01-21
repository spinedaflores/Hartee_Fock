#include<bits/stdc++.h> 
#include "molecule.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>  
#include <iomanip>
#include "blas_lapack_complex.h"
#include <stdio.h>      /* printf */
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <vector>
#include "basis_function.h"
#include "hartree_fock_engine.h"
#include "my_math_functions.h"
 
int main(int argc, char** argv)
{
	//input file read in and print general information
  Molecule mol(argv[1], 0);

  std::cout << "HARTREE FOCK PROGRAM" << std::endl;
  std::cout << "Author: Sergio D. Pineda Flores" << std::endl;
  std::cout << "Last updated: 1/03/2019" << std::endl;
  std::cout << std::endl;
  std::cout << "Number of atoms: " << mol.natom << std::endl;
  std::cout << "Input Cartesian coordinates (Bohr):" << std::endl;

  mol.print_geom_data();

  std::cout << "Beginning Hartree Fock" << std::endl;

	HartreeFockEngine hfeng(mol);

	clock_t t;
	t = clock();

	hfeng.create_matrices();

	t = clock() - t;
	printf ("evaluating 1 and 2 body integrals took %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

	std::cout << "Build the Orthogonalization Matrix S^(-1/2)\n";
	hfeng.build_ortho_mat();

	std::cout << "Build the Initial (Guess) Density\n";
	hfeng.compute_P_initial();

	std::cout << "Compute the Initial SCF Energy\n";
	double old_energy , new_energy;
		hfeng.compute_fock();
		hfeng.compute_P();
		old_energy = hfeng.compute_energy();
		printf("intial SCF energy: %f \n", old_energy);

  std::cout << "Beginning SCF procedure\n" << std::endl;
	t=clock();

	int max_iter = 200;
	int iter = 0;
	new_energy = 0.0;
	while ( iter < max_iter ){
		hfeng.compute_fock();
		hfeng.compute_P();
		new_energy = hfeng.compute_energy();
		printf("SCF energy: %f \n", new_energy);
		if( fabs(new_energy-old_energy) <10E-6 )
			break;
		old_energy = new_energy;
		iter ++;
	}
	std::cout << "~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout << "density matrix\n";
	hfeng.P_mat->print();
	std::cout << std::endl;
	std::cout << "C matrix\n";
	hfeng.C_mat->print();
	std::cout << std::endl;
	std::cout << "F matrix\n";
	hfeng.fock_mat->print();
	std::cout << std::endl;
	std::cout << "Hcore matrix\n";
	hfeng.H_core->print();
	std::cout << std::endl;
	//std::cout << "Compute the New Fock Matrix\n";
	//std::cout << "Build the New Density Matrix\n";
	//std::cout << "Compute the New SCF Energy\n";
	//std::cout << "Test for Convergence\n";
	
	t = clock() - t;
	printf ("SCF procedure took %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

	printf("Final Electronic Energy: %f Hartree\n", new_energy);
	
	std::cout << "Program terminated normally. \n";

  return 0;
}
