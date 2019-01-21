#include <math.h>       /* tgamma */
#include "my_math_functions.h"
#include <complex>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>

// produces cross product of vector1 and vector2 and stores in vector3
void cross_product( const double* const vec1, const double* const vec2, double* const vec3)
{
	double a1,a2,a3,b1,b2,b3;

	a1 = *(vec1);
	a2 = *(vec1+1);
	a3 = *(vec1+2);
	
	b1 = *(vec2);
	b2 = *(vec2+1);
	b3 = *(vec2+2);

	vec3[0] = 	a2*b3 - a3*b2; 
	vec3[1] = -(a1*b3 - a3*b1);
	vec3[2] = 	a1*b2 - a2*b1;

}

long long int fact( long long int n)
{
	return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}


long long int fact2( long long int n)
{
	if(n==-1)
		return 1;
	else{
		long long int step1 = fact(n);
		return fact(step1);
	}
}


double boys(double n, double T)
{
	if (T==0){
		return 1.0/(2*n + 1.0);
	}
	else{
		double g1 = boost::math::tgamma(0.5 + n,0);
		double g2 = boost::math::gamma_q(0.5 + n, T);

		double incomplete_gamma = g1*g2;
		double gamma = g1;

		return 0.5*pow(T, -0.5-n)*( gamma - incomplete_gamma); 
	}

}

