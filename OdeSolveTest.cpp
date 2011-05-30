
#include "GridsCommon.hpp"

#include "intel_ode.h"

void rhs(int *n,double *t,double *y,double *f){
	f[0] = 77.27*(y[1] + y[0]*(1.0 - 8.375*10e-6*y[0] - y[1]));
	f[1] = 1.0/77.27 * (y[2] - (1.0+y[0])*y[1]);
	f[2] = 0.161*(y[0]-y[2]);
}

int main(){

	const int IPAR_SIZE=128;
	int ipar[IPAR_SIZE];
	for(int i=0;i<IPAR_SIZE;i++)ipar[i]=0;
	ipar[0] = 0;
	ipar[1] = 0;//scheme auto-choose
	ipar[2] = 1;//0 - exit at the end ; 1 - exit after every step
	ipar[3] = 0;//0 - autocalc Jacobi matrix; 1 - user-defined Jacobi matrix
	ipar[4] = 0;//0 - don't freeze Jacobi matrix; 1 - freeze Jacobi matrix
	ipar[5] = 0;
	ipar[6] = 0;
	ipar[7] = 0;

	int N=3;
	double t=0;
	double t_end=800;

	const int DPAR_SIZE = max(13*N,(7+2*N)*N);	//As ODE system has size n=2, than the size of dpar array is equal to
	                                            //max{13*n,(7+2*n)*n}=max{26,22}=26. More details on dpar array can be
	                                            //found in the Manual
	double y[N];
	double dpar[DPAR_SIZE];
	int kd[N];

	//double hm=1.e-12; /* minimal step size for the methods */
	//double ep=1.e-6;  /* relative tolerance. The code cannot guarantee the requested accuracy for ep<1.d-9 */
	//double tr=1.e-3;  /* absolute tolerance */
	//double h=1.e-7;

	double hm=1e-3; /* minimal step size for the methods */
	double ep=1e-6;  /* relative tolerance. The code cannot guarantee the requested accuracy for ep<1.d-9 */
	double tr=1e-3;  /* absolute tolerance */
	double h=1.0;

	int ierr=0;


	y[0]=0.01;
	y[1]=0.01;
	y[2]=0.01;
	do{
		cout << t << "    " << y[0] << "    " << y[1]  << "    " << y[2]  << endl;
		dodesol(ipar,&N,&t,&t_end,y,rhs,NULL,&h,&hm,&ep,&tr,dpar,kd,&ierr);
	}while(t<t_end);


	return 0;
}
