#pragma once

void solveODE(){


	int ipar[128];
	int n=5;
	double t;
	double t_end;

	double* y;



	dodesol(ipar,&n,double*,double*,double*,void*,void*,\
			    			 double*,double*,double*,double*,double*,int*,int*);



}
