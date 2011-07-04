#!/bin/bash -e

if [ -z "$1" ]; then
    echo "usage: `basename $0` PLOTDIR"
fi

mkdir $1
mkdir $1/1
mkdir $1/2
mkdir $1/3
mkdir $1/4
mkdir $1/5
mkdir $1/6
mkdir $1/7

((n=0))
sed 's/^ *//g' | normalize.awk | while read l; do
    if echo $l | grep 'TIME:' > /dev/null; then
	time=`echo $l | sed 's/TIME://g;s/ //g'`
	((n=0))
    else
	if echo $l | grep 'Z:' > /dev/null; then
	    z=`echo $l | sed 's/Z://g;s/ //g'`
	else
	    ((cnt=1))
	    for v in $l; do
		printf "$n\t$v\n" >> $1/$cnt/$time.txt
		((cnt++))
	    done
	    ((n++))
	fi
    fi
    
done 