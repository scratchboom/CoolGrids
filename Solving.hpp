#pragma once

//solves f(x)=0
template<typename F>
double solve(F f,double x0,double dx){
	double x=x0;
	double xNext=x0;
	double f_;
	int i=0;
    do{
		x=xNext;
		f_= (f(x+dx)-f(x-dx))/2.0/dx;
		xNext=x - f(x)/f_;
		i++;
		if(i>100)std::cout << "iteration " << i <<std::endl;
		if(i>100)return xNext;
	}while(abs(x-xNext)>0.01*dx);

    return xNext;
};
