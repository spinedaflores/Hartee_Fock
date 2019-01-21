#ifndef HARTREE_FOCK_ENGINE_H 
#define HARTREE_FOCK_ENGINE_H
#include "molecule.h"
#include <vector>
#include <iostream>
#include "basis_function.h"
#include "my_math_functions.h"
#include "my_matrix.h"
#include "blas_lapack_complex.h"

class HartreeFockEngine 
{
	public:
	HartreeFockEngine(Molecule& mol);
	~HartreeFockEngine();

	//number of basis functions = number of H/He atoms
	int nao;
	//total number of molecular orbitals
	int nmo;
	//number of Occupied orbitals
	int o_nmo;
	//number of electrons
	int nel;
	//number of atoms
	int natom;
	//number of 2 electron integrals
	int num_2_el_int;
	//molecule being evaluated
	Molecule* Mol;

	//Matrices to perform hartree fock
	matrix<double>*  overlap_mat;
	matrix<double>*  kinetic_mat;
	matrix<double>*	 H_core;
	matrix<double>*	 fock_mat;
	//an array of potential matrices, one for each atom
	matrix<double>** pot_mat;
	//array of basis functions
	BasisFunction** 	basis_functions;	
	//array of unique 2-electron integrals
	double* two_el_integrals;
	//array of ordered indexes for 2-electron integrals
	int* ioff;
	//orthonalization matrix S^(-1/2)
	matrix<double>*  ortho_mat;
	//expansion coefficients for Molecular Orbitals (MO) matrix
	matrix<double>*	 C_mat;
	//Density matrix
	matrix<double>*	 P_mat;


	double E(int i, int j, int t, double Qx, double a, double b); 
	double overlap(double a, std::vector<int> lmn1, std::vector<double> A, double b, std::vector<int> lmn2, std::vector<double> B); 
	double S(BasisFunction& a, BasisFunction& b);

	double kinetic(double a, std::vector<int>lmn1, std::vector<double> A, double b, std::vector<int> lmn2, std::vector<double> B);
	double T(BasisFunction& a, BasisFunction& b);

	double R(int t, int u, int v, int n, double p, double PCx, double PCy, double PCz, double RPC);

	inline std::vector<double> gaussian_product_center(double a, std::vector<double> A, double b, std::vector<double> B)
	{
		return {(a*A[0]+b*B[0])/(a+b) , (a*A[1]+b*B[1])/(a+b) , (a*A[2]+b*B[2])/(a+b) };
	}

	inline double INDEX(int i,int j){
		 return (i>j) ? (ioff[i]+j) : (ioff[j]+i);
	}

	double nuclear_attraction(double a, std::vector<int> lmn1, std::vector<double> A, double b, std::vector<int> lmn2, std::vector<double> B, std::vector<double> C);
	double V(BasisFunction& a, BasisFunction& b, std::vector<double> C);

	double electron_repulsion(double a,std::vector<int> lmn1, std::vector<double> A,
														double b,std::vector<int> lmn2, std::vector<double> B,
														double c,std::vector<int> lmn3, std::vector<double> C,
														double d,std::vector<int> lmn4, std::vector<double> D);
 double ERI(BasisFunction& a, BasisFunction& b, BasisFunction& c, BasisFunction& d);

 void create_matrices(); 

 void build_ortho_mat();

 void compute_P_initial();
 void compute_P();

 void compute_fock();

 double compute_energy();
};



#endif
