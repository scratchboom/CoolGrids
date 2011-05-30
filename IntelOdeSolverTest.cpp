/*******************************************************************************
!                              INTEL CONFIDENTIAL
!   Copyright(C) 2007-2008 Intel Corporation. All Rights Reserved.
!   The source code contained  or  described herein and all documents related to
!   the source code ("Material") are owned by Intel Corporation or its suppliers
!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!   suppliers and licensors. The Material contains trade secrets and proprietary
!   and  confidential  information of  Intel or its suppliers and licensors. The
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way without Intel's prior express written permission.
!   No license  under any  patent, copyright, trade secret or other intellectual
!   property right is granted to or conferred upon you by disclosure or delivery
!   of the Materials,  either expressly, by implication, inducement, estoppel or
!   otherwise.  Any  license  under  such  intellectual property  rights must be
!   express and approved by Intel in writing.
!
!******************************************************************************
!
!  This example gives the solution of initial value problem for the Van der
!  Pol equation:
!
!        y'-1.d6*[(1-y*y)*y'+1.d6*y=0,  0<t<160,    y(0)=2,  y'(0)=0.
!
!*******************************************************************************/

#include <stdio.h>
#include <time.h>
#include "math.h"
#include "intel_ode.h"

extern void example_v_d_p(int*,double*,double*);
extern void rhs_v_d_p(int*,double*,double*,double*);
extern void jacmat_v_d_p(int*,double*,double*,double*);

int main(void)
{
	int n, ierr, i;


	int kd[2], ipar[128];// It is higly recommended to declare ipar array of size 128 for compatibility with future versions of ODE solvers */
	double t, t_end, h, hm, ep, tr;// As ODE system has size n=2, than the size of dpar array is equal to max{13*n,(7+2*n)*n}=max{26,22}=26. More details on dpar array can be found in the Manual

	double y[2], dpar[26];
	clock_t time_begin,time_end;

    //global parameter settings suitable for all 6 dodesol routines
	hm=1.e-12; // minimal step size for the methods
	ep=1.e-6;  // relative tolerance. The code cannot guarantee the requested accuracy for ep<1.d-9
	tr=1.e-3;  // absolute tolerance

    /****************************** dodesol ********************************/
	for (i=0;i<128;i++) ipar[i]=0;/* Please don't forget to initialize ipar array with zeros before the first call to dodesol routines */

	t=0.e0;
	h=1.e-7;

	example_v_d_p(&n,&t_end,y);//setting size of the system n, end of integration interval t_end, and initial value y at t=0 */

	time_begin=clock();


	dodesol(ipar,&n,&t,&t_end,y,rhs_v_d_p,jacmat_v_d_p,&h,&hm,&ep,&tr,dpar,kd,&ierr);// universal solver
	time_end=clock();

	if(ierr!=0)
	{
		printf("\n========================\n");
		printf("DODESOL C example FAILED\n");
		printf("dodesol routine exited with error code %4d\n",ierr);
		return -1;
	}

	printf("\ndodesol results\n\n");
	printf("ipar[1]=%4d, ipar[3]=%4d\n",ipar[1],ipar[3]);
	printf("t=%5.1f\n",t);
	printf("Solution      y1=%17.14f,   y2=%17.14f\n",y[0],y[1]);
	printf("--------------------------------------------------------\n");
	printf("CPU time=%f seconds\n", ((double)(time_end-time_begin))/CLOCKS_PER_SEC);
	printf("========================================================\n\n");
	if(fabs(y[0]-1.878e0)+fabs(y[1]+0.7436e0)>1.e-2)
		printf("Solution seems to be inaccurate. Probably example FAILED\n");






/*************************** dodesol_rkm9st *****************************/

/* Please don't forget to initialize ipar array with zeros before the first
call to dodesol routines */
	for (i=0;i<128;i++) ipar[i]=0;

	t=0.e0;
	h=1.e-7;

/* setting size of the system n, end of integration interval t_end, and
initial value y at t=0 */
	example_v_d_p(&n,&t_end,y);

	time_begin=clock();
	/* explicit solver */
	dodesol_rkm9st(ipar,&n,&t,&t_end,y,rhs_v_d_p,&h,&hm,&ep,&tr,dpar,&ierr);
	time_end=clock();

	if(ierr!=0)
	{
		printf("\n========================\n");
		printf("DODESOL C example FAILED\n");
		printf("dodesol_rkm9st routine exited with error code %4d\n",ierr);
		return -1;
	}

	printf("\ndodesol_rkm9st results\n\n");
	printf("t=%5.1f\n",t);
	printf("Solution      y1=%17.14f,   y2=%17.14f\n",y[0],y[1]);
	printf("--------------------------------------------------------\n");
	printf("CPU time=%f seconds\n", ((double)(time_end-time_begin))/CLOCKS_PER_SEC);
	printf("========================================================\n\n");
	if(fabs(y[0]-1.878e0)+fabs(y[1]+0.7436e0)>1.e-2)
	{
		printf("Solution seems to be inaccurate. Probably, example FAILED...\n");
		return -1;
	}
/*************************** dodesol_mk52lfn *****************************/

/* Please don't forget to initialize ipar array with zeros before the first
call to dodesol routines */
	for (i=0;i<128;i++) ipar[i]=0;

	t=0.e0;
	h=1.e-7;

/* setting size of the system n, end of integration interval t_end, and
initial value y at t=0 */
	example_v_d_p(&n,&t_end,y);

	time_begin=clock();
/* implicit solver with automatic numerical Jacobi matrix computations */
	dodesol_mk52lfn(ipar,&n,&t,&t_end,y,rhs_v_d_p,&h,&hm,&ep,&tr,dpar,kd,&ierr);
	time_end=clock();

	if(ierr!=0)
	{
		printf("\n========================\n");
		printf("DODESOL C example FAILED\n");
		printf("dodesol_mk52lfn routine exited with error code %4d\n",ierr);
		return -1;
	}

	printf("\ndodesol_mk52lfn results\n\n");
	printf("t=%5.1f\n",t);
	printf("Solution      y1=%17.14f,   y2=%17.14f\n",y[0],y[1]);
	printf("--------------------------------------------------------\n");
	printf("CPU time=%f seconds\n", ((double)(time_end-time_begin))/CLOCKS_PER_SEC);
	printf("========================================================\n\n");
	if(fabs(y[0]-1.878e0)+fabs(y[1]+0.7436e0)>1.e-2)
	{
		printf("Solution seems to be inaccurate. Probably, example FAILED...\n");
		return -1;
	}
/*************************** dodesol_mk52lfa *****************************/

/* Please don't forget to initialize ipar array with zeros before the first
call to dodesol routines */
	for (i=0;i<128;i++) ipar[i]=0;

	t=0.e0;
	h=1.e-7;

/* setting size of the system n, end of integration interval t_end, and
initial value y at t=0 */
	example_v_d_p(&n,&t_end,y);

	time_begin=clock();
/* implicit solver with user-defined Jacobi matrix computations */
	dodesol_mk52lfa(ipar,&n,&t,&t_end,y,rhs_v_d_p,jacmat_v_d_p,&h,&hm,&ep,&tr,dpar,kd,&ierr);
	time_end=clock();

	if(ierr!=0)
	{
		printf("\n========================\n");
		printf("DODESOL C example FAILED\n");
		printf("dodesol_mk52lfa routine exited with error code %4d\n",ierr);
		return -1;
	}

	printf("\ndodesol_mk52lfa results\n\n");
	printf("t=%5.1f\n",t);
	printf("Solution      y1=%17.14f,   y2=%17.14f\n",y[0],y[1]);
	printf("--------------------------------------------------------\n");
	printf("CPU time=%f seconds\n", ((double)(time_end-time_begin))/CLOCKS_PER_SEC);
	printf("========================================================\n\n");
	if(fabs(y[0]-1.878e0)+fabs(y[1]+0.7436e0)>1.e-2)
	{
		printf("Solution seems to be inaccurate. Probably, example FAILED...\n");
		return -1;
	}




/*************************** dodesol_rkm9mkn *****************************/

/* Please don't forget to initialize ipar array with zeros before the first
call to dodesol routines */
	for (i=0;i<128;i++) ipar[i]=0;

	t=0.e0;
	h=1.e-7;

/* setting size of the system n, end of integration interval t_end, and
initial value y at t=0 */
	example_v_d_p(&n,&t_end,y);

	time_begin=clock();
/* hybrid solver with automatic numerical Jacobi matrix computations */
	dodesol_rkm9mkn(ipar,&n,&t,&t_end,y,rhs_v_d_p,&h,&hm,&ep,&tr,dpar,kd,&ierr);
	time_end=clock();

	if(ierr!=0)
	{
		printf("\n========================\n");
		printf("DODESOL C example FAILED\n");
		printf("dodesol_rkm9mkn routine exited with error code %4d\n",ierr);
		return -1;
	}

	printf("\ndodesol_rkm9mkn results\n\n");
	printf("t=%5.1f\n",t);
	printf("Solution      y1=%17.14f,   y2=%17.14f\n",y[0],y[1]);
	printf("--------------------------------------------------------\n");
	printf("CPU time=%f seconds\n", ((double)(time_end-time_begin))/CLOCKS_PER_SEC);
	printf("========================================================\n\n");
	if(fabs(y[0]-1.878e0)+fabs(y[1]+0.7436e0)>1.e-2)
	{
		printf("Solution seems to be inaccurate. Probably, example FAILED...\n");
		return -1;
	}
/*************************** dodesol_rkm9mka *****************************/






/* Please don't forget to initialize ipar array with zeros before the first
call to dodesol routines */
	for (i=0;i<128;i++) ipar[i]=0;

	t=0.e0;
	h=1.e-7;

/* setting size of the system n, end of integration interval t_end, and initial
value y at t=0 */
	example_v_d_p(&n,&t_end,y);

	time_begin=clock();
/* hybrid solver with user-defined Jacobi matrix computations */
	dodesol_rkm9mka(ipar,&n,&t,&t_end,y,rhs_v_d_p,jacmat_v_d_p,&h,&hm,&ep,&tr,dpar,kd,&ierr);
	time_end=clock();

	if(ierr!=0)
	{
		printf("\n========================\n");
		printf("DODESOL C example FAILED\n");
		printf("dodesol_rkm9mka routine exited with error code %4d\n",ierr);
		return -1;
	}

	printf("\ndodesol_rkm9mka results\n\n");
	printf("t=%5.1f\n",t);
	printf("Solution      y1=%17.14f,   y2=%17.14f\n",y[0],y[1]);
	printf("--------------------------------------------------------\n");
	printf("CPU time=%f seconds\n", ((double)(time_end-time_begin))/CLOCKS_PER_SEC);
	printf("========================================================\n\n");
	if(fabs(y[0]-1.878e0)+fabs(y[1]+0.7436e0)>1.e-2)
	{
		printf("Solution seems to be inaccurate. Probably, example FAILED...\n");
		return -1;
	}
	printf("\n========================\n");
	printf("DODESOL C example successfully PASSED through all steps of computations\n");
	return 0;
}

/*************** Data for Van der Pol equations ****************/
void example_v_d_p(int*n,double*t_end,double*y)
/* The routine initializes the size of the system n, the end of
integration interval t_end, and inital data y at t=0.0 */
{
	*n=2;
	*t_end=160.e0;

	y[0]=2.e0;
	y[1]=0.e0;
}

/************* Right hand side of Van der Pol equations ***********/
void rhs_v_d_p(int*n,double*t,double*y,double*f)
{
	double c;

	c=1.e0-y[0]*y[0];

	f[0]=y[1];
	f[1]=(c*y[1]-y[0])*1.0e6;
}

/******* analytical Jacoby matrix for Van der Pol equations *******/
void jacmat_v_d_p(int*n,double*t,double*y,double*a)
{
/* Please make sure that Jacobi matrix is stored in column-wise order:
a[j*n+i]=df(i)/dx(j) */
	a[0]=0.e0;
	a[1]=-1.e6*(1.e0+2.e0*y[0]*y[1]);
	a[2]=1.e0;
	a[3]=1.e6*(1.e0-y[0]*y[0]);
}
/********************* End of C code example **********************/
