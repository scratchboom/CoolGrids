#!/bin/sh -e

if [ -z "$1" ]; then
	echo "usage: `basename $0` fort.501" >&2
	exit 1
fi

for p in `perl -e 'print join " ", map {"A$_"} map {$a=$_;map {"$a$_"} 1..6} 1..6;'` B4; do
	echo $p
	parse-results.sh -l N -A '{print '$p';}' < $1 \
		| awk 'BEGIN{getline;mn=mx=$1}{if(mn>$1){mn=$1}if(mx<$1){mx=$1}}END{printf("%+.3e\n%+.3e\n",mn,mx);}'
done | grep -v 0.000e+00
