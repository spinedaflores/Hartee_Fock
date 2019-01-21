#include "hartree_fock_engine.h"
#include <math.h>       /* pow, erf */
#include <stdio.h>
#include <iterator>

HartreeFockEngine::HartreeFockEngine(Molecule& mol)
{
	Mol = &mol;
	natom = Mol->natom;
	nao = Mol->natom;
	nmo = nao;
	nel = 6;
	o_nmo = nel/2;
	
	int pair = nao*(nao+1)/2;
	num_2_el_int = pair*(pair+1)/2;

	overlap_mat	 = new matrix<double>(nao,nao);
	kinetic_mat	 = new matrix<double>(nao,nao);
	H_core 		 	 = new matrix<double>(nao,nao);
	fock_mat 		 = new matrix<double>(nao,nao);

	pot_mat			 = new matrix<double>* [natom];
	for (int i=0; i<natom; i++)
		pot_mat[i] = new matrix<double>(nao,nao);	

	two_el_integrals = new double [num_2_el_int];
	ortho_mat	 = new matrix<double>(nao,nao);
	C_mat	 = new matrix<double>(nao,nmo);
	P_mat	 = new matrix<double>(nao,nao);

//building list of basis_functions;
// for now I can only create 1s basis functions (sto-3g)
	double x,y,z;
	std::vector<int> 		myShell  = {0,0,0}; // p-orbitals would be (1,0,0) or (0,1,0) or (0,0,1), etc.
	std::vector<double> myExps   = {3.42525091, 0.62391373, 0.16885540};
	std::vector<double> myCoefs  = {0.15432897, 0.53532814, 0.44463454};

	basis_functions = new BasisFunction* [nao]; 	
	for (int i=0; i<nao; i++){
		x = mol.geom[i][0];	
		y = mol.geom[i][1];	
		z = mol.geom[i][2];
		std::vector<double> myOrigin = {x,y,z};

		basis_functions[i] = new BasisFunction(myOrigin, myShell, myExps, myCoefs);
	}

//	std::vector<double> myOrigin1 = {0.0, 0.0, 0.0};
//	std::vector<double> myOrigin2 = {0.0, 0.0, 1.4};
//  std::cout << "creating basis function" << std::endl;
//	BasisFunction a = BasisFunction(myOrigin1, myShell, myExps, myCoefs);
//	BasisFunction b = BasisFunction(myOrigin2, myShell, myExps, myCoefs);
//	std::cout << "testing overlap S(a,a) " << S(a,a) << std::endl;
//	std::cout << "testing overlap S(a,b) " << S(a,b) << std::endl;
//	std::cout << "testing overlap S(b,b) " << S(b,b) << std::endl;
//	std::cout << "testing kinetic T(a,a) " << T(a,a) << std::endl;
//	std::cout << "testing kinetic T(a,b) " << T(a,b) << std::endl;
//	std::cout << "testing kinetic T(b,b) " << T(b,b) << std::endl;
//	std::cout << "testing potential V(a,a,C) on site 1 " << V(a,a,{0,0,0}) << std::endl;
//	std::cout << "testing potential V(a,b,C) on site 1 " << V(a,b,{0,0,0}) << std::endl;
//	std::cout << "testing potential V(b,b,C) on site 1 " << V(b,b,{0,0,0}) << std::endl;
//	std::cout << "testing potential V(a,a,C) on site 2 " << V(a,a,{0,0,1.4}) << std::endl;
//	std::cout << "testing potential V(a,b,C) on site 2 " << V(a,b,{0,0,1.4}) << std::endl;
//	std::cout << "testing potential V(b,b,C) on site 2 " << V(b,b,{0,0,1.4}) << std::endl;
//	std::cout << "testing two_electron integral (aa|aa) " << ERI(a,a,a,a) << std::endl;
//	std::cout << "testing two_electron integral (aa|bb) " << ERI(a,a,b,b) << std::endl;
//	std::cout << "testing two_electron integral (ba|aa) " << ERI(b,a,a,a) << std::endl;
//	std::cout << "testing two_electron integral (ba|ba) " << ERI(b,a,b,a) << std::endl;


}

HartreeFockEngine::~HartreeFockEngine()
{
	delete  overlap_mat;
	delete  kinetic_mat;	 
	delete  H_core; 		 
	delete  fock_mat; 		 

	for (int i=0; i<natom; i++)
		delete  pot_mat[i];	
	delete  pot_mat;

	for (int i=0; i<nao; i++)
		delete basis_functions[i];
	delete [] basis_functions;

	delete [] two_el_integrals;
	delete ortho_mat;
	delete C_mat;
	delete P_mat;

	delete [] ioff;
}
/*
Recursive definition of Hermite Gaussian coefficients.
Returns a float.
a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
i,j: orbital angular momentum number on Gaussian 'a' and 'b'
t: number nodes in Hermite (depends on type of integral, 
	 e.g. always zero for overlap integrals)
Qx: distance between origins of Gaussian 'a' and 'b'
*/
double HartreeFockEngine::E(int i, int j, int t, double Qx, double a, double b) 
{
    double p = a + b;
    double q = (a*b)/p;
    if ( (t < 0) || (t > (i + j)) ){
        //out of bounds for t  
        return 0.0;
		}
    else if (i == j && i== t && i == 0){
        //base case
        return exp(-q*Qx*Qx); // K_AB
		}
    else if (j == 0){
        //decrement index i
        return (1/(2*p))*E(i-1,j,t-1,Qx,a,b) - \
               (q*Qx/a)*E(i-1,j,t,Qx,a,b)    + \
               (t+1)*E(i-1,j,t+1,Qx,a,b);
		}
    else{
        //decrement index j
        return (1/(2*p))*E(i,j-1,t-1,Qx,a,b) + \
               (q*Qx/b)*E(i,j-1,t,Qx,a,b)    + \
               (t+1)*E(i,j-1,t+1,Qx,a,b);
		}
}
/*
Evaluates overlap integral between two Gaussians
Returns a float.
a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
for Gaussian 'a'
lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
B:    list containing origin of Gaussian 'b'
*/
double HartreeFockEngine::overlap(double a, std::vector<int> lmn1, std::vector<double> A, double b, std::vector<int> lmn2, std::vector<double> B)
{ 
    double S1 = E(lmn1[0], lmn2[0], 0, A[0]-B[0], a, b); // X
    double S2 = E(lmn1[1], lmn2[1], 0, A[1]-B[1], a, b); // Y
    double S3 = E(lmn1[2], lmn2[2], 0, A[2]-B[2], a, b); // Z
    return S1*S2*S3*pow(M_PI/(a+b),1.5);
}

/*
Evaluates overlap between two contracted Gaussians
Returns float.
Arguments:
a: contracted Gaussian 'a', BasisFunction object
b: contracted Gaussian 'b', BasisFunction object
*/    
double HartreeFockEngine::S(BasisFunction& a, BasisFunction& b)
{
	double s = 0.0;
	for (int ia=0; ia<a.coefs.size(); ia++){
		for (int ib=0; ib<b.coefs.size(); ib++){
			s += a.norm[ia]*b.norm[ib]*a.coefs[ia]*b.coefs[ib]*\
							 overlap(a.exps[ia],a.shell,a.origin,
							 b.exps[ib],b.shell,b.origin);
			
		}
	}

	return s;
}
/*
Evaluates kinetic energy integral between two Gaussians
Returns a float.
a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0)) for Gaussian 'a'
lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
B:    list containing origin of Gaussian 'b'
*/

double HartreeFockEngine::kinetic(double a, std::vector<int> lmn1, std::vector<double> A, double b, std::vector<int> lmn2, std::vector<double> B)
{
    int l1 = lmn1[0];
		int m1 = lmn1[1];
		int n1 = lmn1[2];

    int l2 = lmn2[0];
		int m2 = lmn2[1];
		int n2 = lmn2[2];

    double term0 = b*(2.0*(l2+m2+n2)+3.0)*\
                            overlap(a, {l1,m1,n1}, A,b, {l2,m2,n2} ,B);

    double term1 = -2.0*pow(b,2.0)*\
                           (overlap(a, {l1,m1,n1}, A,b, {l2+2,m2,n2}, B) +
                            overlap(a, {l1,m1,n1}, A,b, {l2,m2+2,n2}, B) +
                            overlap(a, {l1,m1,n1}, A,b, {l2,m2,n2+2}, B));

    double term2 = -0.5*\
												(l2*(l2-1)*overlap(a, {l1,m1,n1}, A,b, {l2-2,m2,n2}, B) +
												 m2*(m2-1)*overlap(a, {l1,m1,n1}, A,b, {l2,m2-2,n2}, B) +
												 n2*(n2-1)*overlap(a, {l1,m1,n1}, A,b, {l2,m2,n2-2}, B));

//		printf("term0 %f term1 %f term2 %f \n",term0,term1,term2);	
    return term0+term1+term2;
}

/*
Evaluates kinetic energy between two contracted Gaussians
Returns float.
Arguments:
a: contracted Gaussian 'a', BasisFunction object
b: contracted Gaussian 'b', BasisFunction object
*/
double HartreeFockEngine::T(BasisFunction& a, BasisFunction& b)
{
  double t = 0.0;
	for (int ia=0; ia<a.coefs.size(); ia++){
		for (int ib=0; ib<b.coefs.size(); ib++){
			t += a.norm[ia]*b.norm[ib]*a.coefs[ia]*b.coefs[ib]*\
							 kinetic(a.exps[ia],a.shell,a.origin,
							 b.exps[ib],b.shell,b.origin);
			
		}
	}

  return t;
}

/*
Returns the Coulomb auxiliary Hermite integrals 
Returns a float.
Arguments:
t,u,v:   order of Coulomb Hermite derivative in x,y,z
(see defs in Helgaker and Taylor)
n:       order of Boys function 
PCx,y,z: Cartesian vector distance between Gaussian 
composite center P and nuclear center C
RPC:     Distance between P and C
*/
double HartreeFockEngine::R(int t, int u, int v, int n, double p, double PCx, double PCy, double PCz, double RPC)
{
    double T = p*RPC*RPC;
    double val = 0.0;
    if( t==u && t==v && t==0){
        val += pow(-2*p,n)*boys(n,T);
		}
    else if (t == u && t == 0){
        if (v > 1)
            val += (v-1)*R(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC);
        val += PCz*R(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC);
		}
    else if (t == 0){
        if (u > 1)
            val += (u-1)*R(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC);
        val += PCy*R(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC);
		}
    else{
        if (t > 1)
            val += (t-1)*R(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC);
        val += PCx*R(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC);
		}
    return val;
}

/*
Evaluates kinetic energy integral between two Gaussians
Returns a float.
a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
for Gaussian 'a'
lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
B:    list containing origin of Gaussian 'b'
C:    list containing origin of nuclear center 'C'
*/    
double HartreeFockEngine::nuclear_attraction(double a, std::vector<int> lmn1, std::vector<double> A, double b, std::vector<int> lmn2, std::vector<double> B, std::vector<double> C)
{
    int l1 = lmn1[0];
		int m1 = lmn1[1];
		int n1 = lmn1[2];

    int l2 = lmn2[0];
		int m2 = lmn2[1];
		int n2 = lmn2[2];

    double p = a + b;
    std::vector<double> P = gaussian_product_center(a,A,b,B); // Gaussian composite center
    //double RPC = np.linalg.norm(P-C)

	  double accum = 0.0;
    for (int i = 0; i < P.size(); i++) {
				double diff = P[i]-C[i];
				if (abs(diff) < 10E-15)
					diff = 0.0;
        accum += pow( (diff),2);
    }
    double RPC = sqrt(accum);

    double val = 0.0;
    for (int t=0; t <(l1+l2+1); t++){
        for (int u=0; u <(m1+m2+1); u++){
            for (int v=0; v <(n1+n2+1); v++){
                val += E(l1,l2,t,A[0]-B[0],a,b) * \
                       E(m1,m2,u,A[1]-B[1],a,b) * \
                       E(n1,n2,v,A[2]-B[2],a,b) * \
                       R(t,u,v,0,p,P[0]-C[0],P[1]-C[1],P[2]-C[2],RPC);
						}
				}
		}
    val *= 2*M_PI/p; 
    return val;
}

/*
Evaluates overlap between two contracted Gaussians
Returns float.
Arguments:
a: contracted Gaussian 'a', BasisFunction object
b: contracted Gaussian 'b', BasisFunction object
C: center of nucleus
*/
double HartreeFockEngine::V(BasisFunction& a, BasisFunction& b, std::vector<double> C)
{
	double v = 0.0;
	for (int ia=0; ia<a.coefs.size(); ia++){
			for (int ib=0; ib<b.coefs.size(); ib++){
					v += a.norm[ia]*b.norm[ib]*a.coefs[ia]*b.coefs[ib]*\
									 nuclear_attraction(a.exps[ia],a.shell,a.origin,
									 b.exps[ib],b.shell,b.origin,C);
			}
	}

	return -v;
}

/*
Evaluates kinetic energy integral between two Gaussians
Returns a float.
a,b,c,d:   orbital exponent on Gaussian 'a','b','c','d'
lmn1,lmn2
lmn3,lmn4: int tuple containing orbital angular momentum
for Gaussian 'a','b','c','d', respectively
A,B,C,D:   list containing origin of Gaussian 'a','b','c','d'
*/
double HartreeFockEngine::electron_repulsion(double a,std::vector<int> lmn1, std::vector<double> A,
																						 double b,std::vector<int> lmn2, std::vector<double> B,
																						 double c,std::vector<int> lmn3, std::vector<double> C,
																						 double d,std::vector<int> lmn4, std::vector<double> D)
{
	int l1 = lmn1[0]; int l2 = lmn2[0]; int l3 = lmn3[0]; int l4 = lmn4[0];
	int m1 = lmn1[1]; int m2 = lmn2[1]; int m3 = lmn3[1]; int m4 = lmn4[1];
	int n1 = lmn1[2]; int n2 = lmn2[2]; int n3 = lmn3[2]; int n4 = lmn4[2];

	double p = a+b; // composite exponent for P (from Gaussians 'a' and 'b')
	double q = c+d; // composite exponent for Q (from Gaussians 'c' and 'd')
	double alpha = p*q/(p+q);
	std::vector<double> P = gaussian_product_center(a,A,b,B); // A and B composite center
	std::vector<double> Q = gaussian_product_center(c,C,d,D); // C and D composite center
//    RPQ = np.linalg.norm(P-Q)
	double accum = 0.0;
	for (int i = 0; i < P.size(); i++) {
			double diff = P[i]-Q[i];
			if (abs(diff) < 10E-15)
				diff = 0.0;
			accum += pow( (diff),2);
	}
	double RPQ = sqrt(accum);

	double val = 0.0;

		for (int t=0; t<(l1+l2+1); t++){
			for (int u=0; u<(m1+m2+1); u++){
				for (int v=0; v<(n1+n2+1); v++){
					for (int tau=0; tau<(l3+l4+1); tau++){
						for (int nu=0; nu<(m3+m4+1); nu++){
							for (int phi=0; phi<(n3+n4+1); phi++){
													val += E(l1,l2,t	,A[0]-B[0],a,b) * \
																 E(m1,m2,u	,A[1]-B[1],a,b) * \
																 E(n1,n2,v	,A[2]-B[2],a,b) * \
																 E(l3,l4,tau,C[0]-D[0],c,d) * \
																 E(m3,m4,nu ,C[1]-D[1],c,d) * \
																 E(n3,n4,phi,C[2]-D[2],c,d) * \
																 pow(-1,tau+nu+phi) * \
																 R(t+tau,u+nu,v+phi,0,\
																		 alpha,P[0]-Q[0],P[1]-Q[1],P[2]-Q[2],RPQ);
							}
						}
					}	
				}
			}
		}

	val *= 2*pow(M_PI,2.5)/(p*q*sqrt(p+q));
	return val;
}
/*
Evaluates overlap between two contracted Gaussians
Returns float.
Arguments:
a: contracted Gaussian 'a', BasisFunction object
b: contracted Gaussian 'b', BasisFunction object
c: contracted Gaussian 'b', BasisFunction object
d: contracted Gaussian 'b', BasisFunction object
*/
double HartreeFockEngine::ERI(BasisFunction& a, BasisFunction& b, BasisFunction& c, BasisFunction& d)
{
	double eri = 0.0;
	for (int ja=0; ja<a.coefs.size(); ja++){
		for (int jb=0; jb<b.coefs.size(); jb++){
			for (int jc=0; jc<c.coefs.size(); jc++){
				for (int jd=0; jd<d.coefs.size(); jd++){
					eri += a.norm[ja]*b.norm[jb]*c.norm[jc]*d.norm[jd]*\
									 a.coefs[ja]*b.coefs[jb]*c.coefs[jc]*d.coefs[jd]*\
									 electron_repulsion(a.exps[ja],a.shell,a.origin,\
																			b.exps[jb],b.shell,b.origin,\
																			c.exps[jc],c.shell,c.origin,\
																			d.exps[jd],d.shell,d.origin);
				}
			}	
		}
	}

	return eri;
}


void HartreeFockEngine::create_matrices()
{

	std::cout << "creating overlap, kinetic, potential matrices \n";

  for (int i=0; i<nao; i++){
    for(int j=0; j<=i; j++){
			if(i!=j){
				(*overlap_mat)(i,j) = (*overlap_mat)(j,i) = S(*basis_functions[i],*basis_functions[j]);
				(*kinetic_mat)(i,j) = (*kinetic_mat)(j,i) = T(*basis_functions[i],*basis_functions[j]);
				for (int k=0; k<natom; k++){
					(*(pot_mat[k]))(i,j) = (*(pot_mat[k]))(j,i) = V(*basis_functions[i],*basis_functions[j],
																												 {Mol->geom[k][0],Mol->geom[k][1],Mol->geom[k][2]});
				}
			}
			else{
				(*overlap_mat)(i,j) = S(*basis_functions[i],*basis_functions[j]);
				(*kinetic_mat)(i,j) = T(*basis_functions[i],*basis_functions[j]);
				for (int k=0; k<natom; k++){
					(*(pot_mat[k]))(i,j) = V(*basis_functions[i],*basis_functions[j],
																  {Mol->geom[k][0],Mol->geom[k][1],Mol->geom[k][2]});
				}
			}
    }   
  }

	std::cout << "creating H_core matrix \n";
	int			N(H_core->size());
	double	alpha(1.0);
	int			INCX(1);
	int			INCY(1);
	daxpy_(&N,&alpha,kinetic_mat->begin(),&INCX,H_core->begin(),&INCY);
	for (int i=0; i<natom; i++)
		daxpy_(&N,&alpha,(pot_mat[i])->begin(),&INCX,H_core->begin(),&INCY);

//	std::cout << "printing overlap matrix\n";
//	overlap_mat->print();
//	std::cout << "printing kinetic matrix\n";
//	kinetic_mat->print();
//	std::cout << "printing potential matrices\n";
//	for (int i=0; i<natom; i++){
//		std::cout << "	matrix number "<< i << std::endl;
//		(pot_mat[i])->print();
//	}
//	std::cout << "printing H_core matrix\n";
//	H_core->print();

	std::cout << "performing 2 electron integrals \n";

	int index(0);
  for (int i=0; i<nao; i++){
    for(int j=0; j<=i; j++){
      for (int k=0; k<=i; k++){
        if(i==k){
          for(int l=0; l<=j; l++){
//          	printf("	(%d%d|%d%d) \n",i,j,k,l);
						two_el_integrals[index] = ERI(*basis_functions[i],
																					*basis_functions[j],
																					*basis_functions[k],
																					*basis_functions[l]);
						index ++;
          }   
        }   
        else{
          for(int l=0; l<=k; l++){
//          	printf("	(%d%d|%d%d) \n",i,j,k,l);
						two_el_integrals[index] = ERI(*basis_functions[i],
																					*basis_functions[j],
																					*basis_functions[k],
																					*basis_functions[l]);
						index ++;
          }   
        }   
      }   
    }   
  }

	ioff = new int[num_2_el_int];
	ioff[0] = 0;
	for(int i=1; i < num_2_el_int; i++)
		ioff[i] = ioff[i-1] + i;

}

void HartreeFockEngine::build_ortho_mat()
{
	// for dcopy
	int 		SIZE = overlap_mat->size();
	int 		INCX = 1;
	// for dgeev
	char 		JOBVL = 'N';
	char 		JOBVR = 'V';
	int 		N = nao;
	matrix<double> LAMBDA(nao,nao);
	double* A =  LAMBDA.begin();
	int 		LDA = nao;
	double 	WR [nao];
	double 	WI [nao];
	int 		LDVL = 1;
	double* VL; // = new double[N*LDVL];
	int			LDVR = nao;
	//Matrix containing right eigenvectors
	matrix<double> ReigV(LDVR,N);
	double* VR = ReigV.begin();
	int 		LWORK = 6*N;
	double  WORK [LWORK];
	int 		INFO; 
	
	//deep copy the overlap matrix into A, so that the original matrix is untouched
	dcopy_(&SIZE, overlap_mat->begin(),&INCX,A,&INCX);
	// calcuatie the right eignvectors associated with the the overlap matrix
	dgeev_(&JOBVL,&JOBVR,&N,A,&LDA,WR,WI,VL,&LDVL,VR,&LDVR,WORK,&LWORK,&INFO);

  // check for errors
  if (INFO!=0){
    std::cout << "Error: dgeev returned error code " << INFO << std::endl;
  }

//	std::cout << "printing right eigenvectors\n";
//	ReigV.print();
//
//	std::cout << "printing LAMBDA matrix\n";
//	LAMBDA.print();
//
  for (int i=0;i<nao;i++){
    LAMBDA(i,i) = pow( LAMBDA(i,i),-0.5);
  }
//	std::cout << "priting LAMDA^-1/2\n";
//	LAMBDA.print();

	char 		TRANSA = 'N';
	char 		TRANSB = 'N';
	int 		M = nao;
	int 		K = nao;
	double 	ALPHA = 1.0;
					A = ReigV.begin();
					LDA = nao;
	double*	B = LAMBDA.begin();
	int 		LDB = nao;
	double 	BETA = 0.0;
	matrix<double> intermediate_matrix(N,N);
	double*	C = intermediate_matrix.begin();
	int 		LDC = nao;

	// (Right Eigenvectors).(LAMBDA^-0.5)
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

					TRANSB	= 'T';
					A = intermediate_matrix.begin();
					B = ReigV.begin();
					C = ortho_mat->begin();
	// (intermediate_matrix).~(Right Eigenvectors)
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

//	std::cout << "printing ortho_mat\n";
//	ortho_mat->print();

}	

void HartreeFockEngine::compute_P_initial()
{
	char 		TRANSA = 'T';
	char 		TRANSB = 'N';
	int 		M = nao;
	int	 		N = nao;
	int 		K = nao;
	double 	ALPHA = 1.0;
	double*	A = ortho_mat->begin();
	int			LDA = nao;
	double*	B = H_core->begin();
	int 		LDB = nao;
	double 	BETA = 0.0;
	matrix<double> intermediate_matrix(nao,nao);
	double*	C = intermediate_matrix.begin();
	int 		LDC = nao;

	// ~S^(-1/2).H_core
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

					TRANSA = 'N';
					A = intermediate_matrix.begin();
					B = ortho_mat->begin();
					C = fock_mat->begin();		
	// Intermediate_matrix.S^(-1/2)
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

//	std::cout << "printing F' initial \n";
//	fock_mat-> print();

	char 		JOBVL = 'N';
	char 		JOBVR = 'V';
	double 	WR [nao];
	double 	WI [nao];
	int 		LDVL = 1;
	double* VL; // = new double[N*LDVL];
	int			LDVR = nao;
	//Matrix containing right eigenvectors
	matrix<double> ReigV(LDVR,N);
	double* VR = ReigV.begin();
	int 		LWORK = 6*N;
	double 	WORK [LWORK];
	int 		INFO; 

	int SIZE = fock_mat->size();
	int INCX = 1;
	double F [SIZE];
	// deep copy of fock matrix
	dcopy_(&SIZE, fock_mat->begin(),&INCX,F,&INCX);
	// Eigensystem of initial Fock Matrix
	dgeev_(&JOBVL,&JOBVR,&N,F,&LDA,WR,WI,VL,&LDVL,VR,&LDVR,WORK,&LWORK,&INFO);

					TRANSA = TRANSB = 'N';
					M = N = K = LDA = LDB = LDC = nao;
					ALPHA = 1.0; BETA = 0.0;
					A = ortho_mat->begin();
					B = ReigV.begin();
					C = C_mat->begin();		
	// S^(-1/2).C'
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

////	for (int m=0; m<o_nmo; m++){
////		double accum = 0.0;
////		for (int nu = 0; nu < nao; nu++) {
////				accum += pow( (*C_mat)(nu,nmo-1-m),2);
////		}
////		double scale = 1.0/sqrt(accum);
////		int incr=1;	
////
////		dscal_(&nao,&scale,&(*C_mat)(0,nmo-1-m),&incr);
////	}
///	(*C_mat)(0,1)= 0.57735;
///	(*C_mat)(1,1)= 0.816497;
///	std::cout << "printing C0 \n";
///	C_mat->print();

	std::fill(P_mat->begin(), P_mat->end(), 0.0);
	//make more efficient later with BLAS calls and taking advantage of sym
	for(int mu=0; mu<nao; mu++){
		for(int nu=0; nu<nao; nu++){
			for(int m=0; m < o_nmo; m++){
				(*P_mat)(mu,nu) += 2*(*C_mat)(mu,m)*(*C_mat)(nu,m);  
//				(*P_mat)(mu,nu) += 2*(*C_mat)(mu,nmo-1-m)*(*C_mat)(nu,nmo-1-m);  
			}
		}
	}

//	std::cout << "printing initial density matrix\n";
//	P_mat->print();
}

void HartreeFockEngine::compute_P()
{
	char 		TRANSA = 'T';
	char 		TRANSB = 'N';
	int 		M = nao;
	int	 		N = nao;
	int 		K = nao;
	double 	ALPHA = 1.0;
	double*	A = ortho_mat->begin();
	int			LDA = nao;
	double*	B = fock_mat->begin();
	int 		LDB = nao;
	double 	BETA = 0.0;
	matrix<double> intermediate_matrix(nao,nao);
	double*	C = intermediate_matrix.begin();
	int 		LDC = nao;

	// ~S^(-1/2).F
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

					TRANSA = 'N';
					A = intermediate_matrix.begin();
					B = ortho_mat->begin();
	//second intermediate matrix
	int SIZE = fock_mat->size();
	double F [SIZE];
	// Intermediate_matrix.S^(-1/2)
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,F,&LDC);

//	std::cout << "printing F compute_P' \n";
//	fock_mat-> print();

	char 		JOBVL = 'N';
	char 		JOBVR = 'V';
	double  WR [N];
	double  WI [N];
	int 		LDVL = 1;
	double* VL; // = new double[N*LDVL];
	int			LDVR = nao;
	//Matrix containing right eigenvectors
	matrix<double> ReigV(LDVR,N);
	double* VR = ReigV.begin();
	int 		LWORK = 6*N;
	double WORK [LWORK];
	int 		INFO; 
	
	//computing C' from F
	dgeev_(&JOBVL,&JOBVR,&N,F,&LDA,WR,WI,VL,&LDVL,VR,&LDVR,WORK,&LWORK,&INFO);

					TRANSA = TRANSB = 'N';
					M = N = K = LDA = LDB = LDC = nao;
					ALPHA = 1.0; BETA = 0.0;
					A = ortho_mat->begin();
					B = ReigV.begin();
					C = C_mat->begin();		
	//comput C = S^(-1/2)C'
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

////	for (int m=0; m<o_nmo; m++){
////		double accum = 0.0;
////		for (int nu = 0; nu < nao; nu++) {
////				accum += pow( (*C_mat)(nu,nmo-1-m),2);
////		}
////		double scale = 1.0/sqrt(accum);
////		int incr=1;	
////
////		dscal_(&nao,&scale,&(*C_mat)(0,nmo-1-m),&incr);
////	}

	std::fill(P_mat->begin(), P_mat->end(), 0.0);
	//make more efficient later with BLAS calls and taking advantage of sym
	for(int mu=0; mu<nao; mu++){
		for(int nu=0; nu<nao; nu++){
			for(int m=0; m < o_nmo; m++){
				(*P_mat)(mu,nu) += 2*(*C_mat)(mu,m)*(*C_mat)(nu,m);  
//				(*P_mat)(mu,nu) += 2*(*C_mat)(mu,nmo-1-m)*(*C_mat)(nu,nmo-1-m);  
			}
		}
	}

//	std::cout << "printing F @ end of compute_P end \n";
//	fock_mat-> print();
}

void HartreeFockEngine::compute_fock()
{
//reset fock matrix to zeros
std::fill(fock_mat->begin(), fock_mat->end(), 0.0);

int ij, kl, ijkl, ik, jl, ikjl;
	for(int i=0; i < nao; i++){
		for(int j=0; j < nao; j++) {
			(*fock_mat)(i,j) += (*H_core)(i,j);
			for(int k=0; k < nao; k++){
				for(int l=0; l < nao; l++) {
					ij = INDEX(i,j);
					kl = INDEX(k,l);
					ijkl = INDEX(ij,kl);
					ik = INDEX(i,k);
					jl = INDEX(j,l);
					ikjl = INDEX(ik,jl);
 
					(*fock_mat)(i,j) += (*P_mat)(k,l) *\
															(two_el_integrals[ijkl] - 0.5*two_el_integrals[ikjl]);
				}
			}
		}
 	}

//std::cout << "F matrix @ compute_fock\n";
//fock_mat->print();
//std::cout << std::endl;

}

double HartreeFockEngine::compute_energy()
{
	double energy = 0.0;
	double nn_repulsion = 0.0;
	
//std::cout << "F matrix for energy\n";
//fock_mat->print();
//std::cout << std::endl;

	for(int mu=0; mu<nao; mu++){
		for(int nu=0; nu<nao; nu++){
			energy += (*P_mat)(mu,nu) * ( (*H_core)(nu,mu) + (*fock_mat)(nu,mu) );
		}
	}

	return 0.5*energy + nn_repulsion; 

}





