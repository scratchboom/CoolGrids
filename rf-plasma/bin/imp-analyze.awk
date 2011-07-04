#!/usr/bin/awk -f
BEGIN {
	printf "enter <freq> <n>: ";
}
{
	h=6./40000;
	tau=0.500500525e-12;
	c=3e8;
	
	f=$1;
	n=$2;
	
	T=1/f;
	l=c/f;
	print "------------------";
	print "FREQ: " f "\tn: " n;
	print "";
	print "points: " l/h;
	print "lambda: " l;
	print "time steps: " T/tau;
	print "period: " T;
	print "";
	print "MK: " l/h*n;
	print "LEN: " l*n;
	print "N: " T/tau*n;
	print "T: " T*n;
	printf "\n\nenter <freq> <n>: ";
}
