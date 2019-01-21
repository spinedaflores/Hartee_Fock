#include "molecule.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include "blas_lapack_complex.h"
#include "my_math_functions.h"
#include "my_matrix.h"
#include <vector>

Molecule::Molecule(const char* const filename, int const q)
{
  charge = q;

  // open filename
  std::ifstream is(filename);
  assert(is.good());

  // read the number of atoms from filename
  is >> natom;

  // allocate space
  zvals = new int[natom];
  geom = new double* [natom];
  for(int i=0; i < natom; i++)
    geom[i] = new double[3];

  for(int i=0; i < natom; i++)
    is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

  is.close();
}

Molecule::~Molecule()
{
  delete[] zvals;
  for(int i=0; i < natom; i++)
    delete[] geom[i];
  delete[] geom;
}

//prints the geometry of the input molecule
void Molecule::print_geom_data()
{
  for(int i=0; i < natom; i++)
    printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);

	if (natom > 1){
		std::cout << "Interatomic distances (Bohr):\n";
		for(int i=0; i < natom; i++)
			for(int j=0; j < i; j++)
				printf("%d %d %8.5f\n", i, j, bond(i,j));
	}

	if (natom > 2){
		std::cout << "\nBond angles:\n";
		for(int i=0; i < natom; i++) {
			for(int j=0; j < i; j++) {
				for(int k=0; k < j; k++) {
					if(bond(i,j) < 4.0 && bond(j,k) < 4.0)
						printf("%2d-%2d-%2d %10.6f\n", i, j, k, angle(i,j,k)*(180.0/std::acos(-1.0)));
				}
			}
		}
	}

	if(natom > 3){
		std::cout << "\nOut-of-Plane angles:\n";
		for(int i=0; i < natom; i++) {
			for(int k=0; k < natom; k++) {
				for(int j=0; j < natom; j++) {
					for(int l=0; l < j; l++) {
						if(i!=j && i!=k && i!=l && j!=k && k!=l && bond(i,k) < 4.0 && bond(k,j) < 4.0 && bond(k,l) < 4.0)
								printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, oop(i,j,k,l)*(180.0/acos(-1.0)));
					}
				}
			}
		}
	
		std::cout << "\nTorsional angles:\n\n";
		for(int i=0; i < natom; i++) {
			for(int j=0; j < i; j++) {
				for(int k=0; k < j; k++) {
					for(int l=0; l < k; l++) {
						if(bond(i,j) < 4.0 && bond(j,k) < 4.0 && bond(k,l) < 4.0)
							printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, torsion(i,j,k,l)*(180.0/acos(-1.0)));
					}
				}
			}
		}

	}

	//atomic number to masses  
  double an2masses []= {0.0,1.008,4.003,6.941,9.012,10.811,12.011,14.007,15.999,18.998,20.18,22.99,24.305,26.982,28.086,30.974,32.065,35.453,39.948,39.098,40.078,44.956,47.867,50.942,51.996,54.938,55.845,58.933,58.693,63.546,65.39,69.723,72.64,74.922,78.96,79.904,83.8,85.468,87.62,88.906,91.224,92.906,95.94,98,101.07,102.906,106.42,107.868,112.411,114.818,118.71,121.76,127.6,126.905,131.293,132.906,137.327,138.906,140.116,140.908,144.24,145,150.36,151.964,157.25,158.925,162.5,164.93,167.259,168.934,173.04,174.967,178.49,180.948,183.84,186.207,190.23,192.217,195.078,196.967,200.59,204.383,207.2,208.98,209,210,222,223,226,227,232.038,231.036,238.029,237,244,243,247,247,251,252,257,258,259,262,261,262,266,264,277,268};


  /* find the center of mass (COM) */
  double M = 0.0;
  for(int i=0; i < natom; i++)
    M += an2masses[(int) zvals[i]];
 
  double xcm=0.0;
  double ycm=0.0;
  double zcm=0.0;
  double mi;
  for(int i=0; i < natom; i++) {
    mi = an2masses[(int) zvals[i]];
    xcm += mi * geom[i][0];
    ycm += mi * geom[i][1];
    zcm += mi * geom[i][2];
  }
  xcm /= M;
  ycm /= M;
  zcm /= M;
  printf("\nMolecular center of mass: %12.8f %12.8f %12.8f\n", xcm, ycm, zcm);
 
//  translate(-xcm, -ycm, -zcm);

	//creating a column majored matrix whose elements can be accessed via I[row_index + (column_index * rows_per_column)]
	matrix <double> I(3,3);
	
  for(int i=0; i < natom; i++) 
	{
    mi = an2masses[(int) zvals[i]];
    I(0,0) += mi * (geom[i][1]*geom[i][1] + geom[i][2]*geom[i][2]);
    I(1,1) += mi * (geom[i][0]*geom[i][0] + geom[i][2]*geom[i][2]);
    I(2,2) += mi * (geom[i][0]*geom[i][0] + geom[i][1]*geom[i][1]);
    I(0,1) += mi * geom[i][0]*geom[i][1];
    I(0,2) += mi * geom[i][0]*geom[i][2];
    I(1,2) += mi * geom[i][1]*geom[i][2];
  }
 
  I(1,0) = I(0,1);
  I(2,0) = I(0,2);
  I(2,1) = I(1,2);


//	I.print();
	//check that I is a square matrix
	assert(I.row==I.col);

  // allocate data for eigenvalue solver
	double* data = I.begin();
	int n = I.row;
  char Nchar='N';
  double *eigReal=new double[n];
  double *eigImag=new double[n];
  double *vl,*vr;
  int one=1;
  int lwork=6*n;
  double *work=new double[lwork];
  int info;

  // calculate eigenvalues using the DGEEV subroutine
  dgeev_(&Nchar,&Nchar,&n,data,&n,eigReal,eigImag,vl,&one,vr,&one,work,&lwork,&info);
	
  // check for errors
  if (info!=0){
    std::cout << "Error: dgeev returned error code " << info << std::endl;
  }

  // find the principal moments
  std::cout << "\nPrincipal moments of inertia (amu *Bohr^2):\n";
  for (int i=0;i<n;i++){
    std::cout << eigReal[i] << " ";
  }
  std::cout << std::endl;
 
  // classify the rotor 
  if(natom == 2)
		std::cout << "\nMolecule is diatomic.\n";
  else if(eigReal[0] < 1e-4)
		std::cout << "\nMolecule is linear.\n";
  else if((fabs(eigReal[0] - eigReal[1]) < 1e-4) && (fabs(eigReal[1] - eigReal[2]) < 1e-4))
    std::cout << "\nMolecule is a spherical top.\n";
  else if((fabs(eigReal[0] - eigReal[1]) < 1e-4) && (fabs(eigReal[1] - eigReal[2]) > 1e-4))
    std::cout << "\nMolecule is an oblate symmetric top.\n";
  else if((fabs(eigReal[0] - eigReal[1]) > 1e-4) && (fabs(eigReal[1] - eigReal[2]) < 1e-4))
    std::cout << "\nMolecule is a prolate symmetric top.\n";
  else 
		std::cout << "\nMolecule is an asymmetric top.\n"; 

  // deallocate
  delete [] eigReal;
  delete [] eigImag;
  delete [] work;

}

//translates the entire molecule 
void Molecule::translate(double x, double y, double z)
{
  for(int i=0; i < natom; i++) {
     geom[i][0] += x;
     geom[i][1] += y;
     geom[i][2] += z;
  }
}

//returns the bond length between atom1 and atom2
double Molecule::bond(int atom1, int atom2)
{
  double sum_squares(0);
  
  for (int i=0; i<3; i++)
    sum_squares += pow(geom[atom2][i]-geom[atom1][i],2); 

  return std::sqrt(sum_squares);
}
 
//returns a component (cart=x,y,z) of the unit vector pointing from atom_i to atom_j
double Molecule::unit(int cart, int i, int j)
{
  double R_ji = bond(j, i);
  double e_ji = (geom[j][cart] - geom[i][cart]) / R_ji;

  return e_ji;
}

//returns the angle formed by atoms i,j,k
double Molecule::angle(int atom_i, int atom_j, int atom_k)
{
  int dim(3);
  int incr(1);

  double e_ji[dim];
  double e_jk[dim]; 

  for(int i=0; i<dim; i++)
    e_ji[i] = unit(i, atom_j,atom_i);

  for(int i=0; i<dim; i++)
    e_jk[i] = unit(i, atom_j,atom_k);

  double cos_phi_ijk;
  
  cos_phi_ijk = ddot_(&dim, &(e_ji[0]), &incr, &(e_jk[0]), &incr); 

  return acos(cos_phi_ijk);  
}

//out of plane angle
double Molecule::oop(int atom_i,int atom_j,int atom_k,int atom_l)
{
  int dim(3);
  int incr(1);

  double e_kj[dim];
  double e_kl[dim]; 
  double e_ki[dim];  
  // cross product of e_kj and e_kl
  double e_crss_prdct[dim];

  for(int i=0; i<dim; i++)
  {
    e_kj[i] = unit(i, atom_k,atom_j);
    e_kl[i] = unit(i, atom_k,atom_l);
    e_ki[i] = unit(i, atom_k,atom_i);
  }

  cross_product( &(e_kj[0]), &(e_kl[0]), &(e_crss_prdct[0]) );

  double theta = ddot_(&dim, &(e_crss_prdct[0]), &incr, &(e_ki[0]), &incr)/sin(angle(atom_j, atom_k, atom_l));

  if(theta < -1.0) theta = asin(-1.0);
  else if(theta > 1.0) theta = asin(1.0);
  else theta = asin(theta);

  return theta;
}

double Molecule::torsion(int atom_i, int atom_j, int atom_k, int atom_l)
{
  int dim(3);
  int incr(1);
  
  double e_ij[dim];
  double e_jk[dim];
  double e_kl[dim];
  //cross products of the above vectors
  double e_ij_jk[dim];
  double e_jk_kl[dim];

  for(int i=0; i<dim; i++)
  {
    e_ij[i] = unit(i, atom_i,atom_j); 
    e_jk[i] = unit(i, atom_j,atom_k);
    e_kl[i] = unit(i, atom_k,atom_l);
  }

  cross_product( &(e_ij[0]), &(e_jk[0]), &(e_ij_jk[0]) );
  cross_product( &(e_jk[0]), &(e_kl[0]), &(e_jk_kl[0]) );

  double denominator = (sin(angle(atom_i,atom_j,atom_k)) * sin(angle(atom_j,atom_k,atom_l)));

  double numerator = ddot_(&dim, &(e_ij_jk[0]), &incr, &(e_jk_kl[0]), &incr);
  double tau = numerator/denominator;

  if(tau < -1.0) tau = acos(-1.0);
  else if(tau > 1.0) tau = acos(1.0);
  else tau = acos(tau);

  if(numerator < 0.0)
    return -tau;
  else
    return tau;    
}



