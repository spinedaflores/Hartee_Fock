#ifndef BASIS_FUNCTION_H
#define BASIS_FUNCTION_H

#include <vector>

/*
		A class that contains all our basis function data
		Attributes:
		origin: array/list containing the coordinates of the Gaussian origin
		shell:  tuple of angular momentum
		exps:   list of primitive Gaussian exponents
		coefs:  list of primitive Gaussian coefficients
		norm:   list of normalization factors for Gaussian primitives
*/

class BasisFunction
{
	public:
	std::vector<int> shell;
	std::vector<double> origin;
	std::vector<double> exps;
	std::vector<double> coefs;
	std::vector<double> norm;

	BasisFunction(std::vector<double> origin_info,
								std::vector<int> shell_info,
								std::vector<double> exps_info,
								std::vector<double> coefs_info);

	~BasisFunction();

	double norm_prim(double exponent, int l, int m, int n);
	void normalize();


};

#endif
