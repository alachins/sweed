//these three parameters affect precision of optimization.
//m is not recommended to be higher than 20
//decreasing factr, pgtol will increase precision
//factr is the multiple of machine precision that result will be
//pgtol is size of gradient on exit
#define MVAL 12
#define FACTR 1.0e6
#define PGTOL 1.0e-3


//nbd is a vector of integers of dimension numpars.
//nbd[i]=0 if there are no bounds for parameter i, 
//      =1 if there are only lower bounds
//      =2 if there are both lower/upper bounds
//      =3 if there is only upper bound                
//or send nbd=NULL is equivalent to nbd=2 for all parameters
//
//noisy=0 => no output, noisy=1 => one line of output, noisy < 99 some output,
//noisy>=100 probably too much
//
//dfun is derivative function or send NULL to use numerical derivative
//(getgradient function)
double findmax_bfgs(int numpars, double *invec, double (*fun)(const double x[]),
		    void (*dfun)(const double x[], double y[]),
		    double *lowbound, double *upbound,
		    int *nbd, int noisy);

void getgradient(int npar, const double invec[], const int need_gradient[], 
		 double outvec[], double(*func)(const double []), 
		 const double* lowbound, const double* upbound);
