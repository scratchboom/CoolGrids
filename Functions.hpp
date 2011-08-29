#pragma once


double sqr(double x){
	return x*x;
}

double sqr(double x0,double x1){
	return sqr(x0)+sqr(x1);
}

double sqr(double x0,double x1,double x2){
	return sqr(x0)+sqr(x1)+sqr(x2);
}

double length(double x,double y){
	return sqrt(sqr(x,y));
}

double length(double x,double y,double z){
	return sqrt(sqr(x,y,z));
}

double cube(double x){
	return x*x*x;
}

double quad(double x){
	return sqr(sqr(x));
}


const double SQRT_LN_2=0.8325546111576977563531646;
//http://www4d.wolframalpha.com/Calculate/MSP/MSP475119ge38fab96e21ac0000571gibbhh4e7aa5b?MSPStoreType=image/gif&s=49&w=401&h=200&cdf=Coordinates&cdf=Tooltips
double gauss(double x,double x0,double w){
    return exp(-sqr((x-x0)*2.0*SQRT_LN_2/w));
}

double gaussStep(double x,double x0,double w){
    if(x<x0)return gauss(x,x0,w);
    else return 1.0;
}

double gaussWindow(double x,double x0,double x1,double w){
    if(x<x0)return gauss(x,x0,w);
    else if(x>x1)return gauss(x,x1,w);
    else return 1.0;
}

double max(double d1,double d2){
	return d1>d2 ? d1 : d2;
}


double max(double d1,double d2,double d3){
    return max(d1,max(d2,d3));
}

double max(double d1,double d2,double d3,double d4,double d5,double d6){
    return max(max(d1,d2,d3),max(d4,d5,d6));
}

double maxAbs(double d1,double d2){
    return max(abs(d1),abs(d2));
}

double maxAbs(double d1,double d2,double d3){
    return max(abs(d1),abs(d2),abs(d3));
}

double maxAbs(double d1,double d2,double d3,double d4,double d5,double d6){
    return max(abs(d1),abs(d2),abs(d3),abs(d4),abs(d5),abs(d6));
}


double min(double d1,double d2){
	return d1<d2 ? d1 : d2;
}


double min(double d1,double d2,double d3){
    return min(d1,min(d2,d3));
}

double min(double d1,double d2,double d3,double d4,double d5,double d6){
    return min(min(d1,d2,d3),min(d4,d5,d6));
}

double minAbs(double d1,double d2){
    return min(abs(d1),abs(d2));
}

double minAbs(double d1,double d2,double d3){
    return min(abs(d1),abs(d2),abs(d3));
}

double minAbs(double d1,double d2,double d3,double d4,double d5,double d6){
    return min(abs(d1),abs(d2),abs(d3),abs(d4),abs(d5),abs(d6));
}


double avg(double x1,double x2){
	return (x1+x2)/2.0;
}

double avg(double x1,double x2,double x3){
	return (x1+x2+x3)/3.0;
}

double avg(double x1,double x2,double x3,double x4){
	return (x1+x2+x3+x4)/4.0;
}

bool isEveryNth(double n,int N){
	return ((int)n)%N==0;
}

bool between(double x,double a,double b){
    return (x>=a) && (x<=b);
}

double sign(double a){
	if(a>0)return 1;
	else if(a<0)return -1;
	else return 0;
}

double sign(double a,double b){
	return abs(a)*sign(b);
}

double minMod(double y,double z){
	return 0.5 * (sign(y)+sign(z)) * min(abs(y),abs(z));
}

double cosPulse(double x){
	if(abs(x) < 0.5*M_PI)return cos(x);
	else return 0.0;
}


double сosPulse(double x){
	if(abs(x) < 0.5*M_PI)return cos(x);
	else return 0.0;
}

double onePlusCosPulse(double x){
	if(abs(x) < M_PI)return 0.5*(1.0+cos(x));
	else return 0.0;
}

