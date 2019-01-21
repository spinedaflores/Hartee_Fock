#include "basis_function.h"
#include <stdio.h>      /* printf */
#include <math.h>       /* pow, erf */
#include "my_math_functions.h"
#include <iostream>
#include <stdlib.h>     /* exit, EXIT_FAILURE */


BasisFunction::BasisFunction(std::vector<double> origin_info,
														 std::vector<int> shell_info,
														 std::vector<double> exps_info,
														 std::vector<double> coefs_info)
{


//  printf ("stopping here to debug\n");
//  exit (EXIT_FAILURE);
	origin 	= origin_info;
	shell 	= shell_info;
	exps 		= exps_info;
	coefs		= coefs_info;

	normalize();
}

BasisFunction::~BasisFunction()
{
}

double BasisFunction::norm_prim(double exponent, int l, int m, int n)
{
	return sqrt(pow(2, 2*(l+m+n)+1.5)*
									pow(exponent,l+m+n+1.5)/
									fact2(2*l-1)/fact2(2*m-1)/
									fact2(2*n-1)/pow(M_PI,1.5));
}

void BasisFunction::normalize()
{
	// Routine to normalize the basis functions, in case they
	// do not integrate to unity.
	int l,m,n,L;
	l = shell[0]; 
	m = shell[1];
	n = shell[2];

	L = l+m+n;

	int num_exps = exps.size();
	norm.resize(num_exps);
	// norm is a list of length equal to number primitives
	// normalize primitives first (PGBFs)
	for (int i=0; i < num_exps; i++)
		norm[i] = norm_prim(exps[i], l,m,n);

	// now normalize the contracted basis functions (CGBFs)
	// Eq. 1.44 of Valeev integral whitepaper
	double prefactor = pow(M_PI,1.5)*\
				 fact2(2*l - 1)*fact2(2*m - 1)*fact2(2*n - 1)/pow(2.0,L);

	double N = 0.0;
	for (int ia = 0; ia < num_exps; ia++){
		for (int ib = 0; ib < num_exps; ib++){
					N += norm[ia]*norm[ib]*coefs[ia]*coefs[ib]/\
									 pow(exps[ia] + exps[ib],L+1.5);
		}
	}

	N *= prefactor;
	N = pow(N,-0.5);
	for (int ia = 0; ia < num_exps; ia++)
		coefs[ia] *= N;


}
