#!/bin/sh -e

LEN=`head -3 parskhgs.cfg | tail -1 | perl -ne 's/d/e/i; /\s*([0-9.\-+ed]*)/i; print $1;'`
MK=`head -1 parskhgs.cfg | perl -ne 's/d/e/i; /\s*([0-9.\-+ed]*)/i; print $1;'`

parse-results.sh '
# Foreach data point
if (m==1 && N>1) { intS += Ex*Hy*dT*dx2; }
if (m==points-1 && N>1) { intS -= Ex*Hy*dT*dx2; }
 
Ee += ne*dx3*Et*evlt + ne*dx3*cme*(vx*vx+vz*vz)/2;
Ef += (e0*(Ex*Ex+Ez*Ez) + mu0*(Hy*Hy))*dx3/2;
Ei += ne*dx3*10*evlt;

Px += dx3*ne*cme*vx;
Pz += dx3*ne*cme*vz;

' '
# Foreach data line
if (N>1) {
  print (N-1) "\t" pT "\t" Ef "\t" Ef-pEf "\t" Ee "\t" intS "\t" Ei;
}
Ef=Ee=Px=Pz=0;
dT=T-pT;
pT=T;
pEf
' '
# BEGIN
len='"$LEN"'
points='"$MK"'
cme=9.10938215e-31;
cze=evlt=1.602176487e-19;
e0=8.854187817e-12;
PI=3.14159265;
mu0=4*PI*1e-7;
dx=len/points;
dx2=dx*dx;
dx3=dx*dx*dx;
' '
# END
'
