#!/usr/bin/awk -f
BEGIN{
    line=-2;
}
/^ *$/{line=-1;}
{
    if (line>0) {
	if ((line-1)%3==0) { u[1] = $1; u[2] = $2; u[3] = $3; }
	if ((line-1)%3==1) { u[4] = $1; u[5] = $2; u[6] = $3; }
	if ((line-1)%3==2) { u[7] = $1; }
	if ((line-1)%3==2) {
	    print u[1] " " u[2] " " u[3] " " u[4] " " u[5] " " u[6] " " u[7];
	}
	line++;
    }
    if (line==0) {
	print "TIME: " $1;
	print "Z: " $2;
	line++;
    }
    if (line==-1) line++;
}
