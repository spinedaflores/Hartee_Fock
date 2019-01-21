#ifndef BLAS_LAPACK_COMPLEX_H 
#define BLAS_LAPACK_COMPLEX_H 
#include <random>
#include <complex>

extern "C" 
{
	extern int 												dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
	extern double 										ddot_(int*,double*,int*,double*,int*);
	extern std::complex<double> 			zdotc_(int*,std::complex<double>*,int*,std::complex<double>*,int*);
	extern void												daxpy_(int*,double*,double*,int*,double*,int*);
	extern void												dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
	extern void												dscal_(int*,double*,double*,int*);
	extern void 											dcopy_(int*,double*,int*,double*,int*);

}



#endif

