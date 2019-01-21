# Hartee_Fock
Hartree Fock code with parallel features

Notes about current implementation of code:
	*for now only systems of 1s orbitals should be run. Generalization will be done later
	*be careful of numerical precision issues that can arise...
		-it can occur due to the number of 2 electron integrals becoming very large
		-the orthogonalization procedure can become unstable when there is near linear dependence in the basis set


current characteristic to be kept in mind 
	*eignvectors lose normalilzation in conversion from C' to C
	*ordering of eigenvalues from LAPACK isn't always from most negative to most positive
