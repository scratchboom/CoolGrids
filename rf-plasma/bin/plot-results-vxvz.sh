#!/bin/sh -e

usage () {
    echo "plot vz(vx) for each result in RESULTS and for each m in POINTS"
    echo "usage: `basename $0` RESULTS [POINTS] [XRANGE] [YRANGE]"
    echo "       RANGE is 'from:to' expression"
}

if [ -z "$1" ]; then
    usage 1>&2
    exit 1
fi

RESULTS=$1
if [ -z "$2" ]; then ARGS="160 165 170 175 180 185 190 195 200 205 210 215 220 225 230 235 240 245"; else ARGS="$2"; fi
if [ -z "$3" ]; then XRANGE="-4e6:8e6"; else XRANGE="$3"; fi
if [ -z "$4" ]; then YRANGE="-8e8:7e8"; else YRANGE="$4"; fi

for r in $RESULTS; do (
	for m in $ARGS; do
	    parse-results.sh 'if(m=='$m') {print vx "\t" vz}' < $r/results.txt \
		| plt -f -t 'vz(vx)(m'$m')' \
		$r/'vz(vx)(m'$m').png' \
		'set xrange ['"$XRANGE"'];
                 set yrange['"$YRANGE"'];
                 set arrow 1 from -3e8,3e8 to 3e8,3e8 nohead lt 0;
                 set arrow 2 from -3e8,-3e8 to 3e8,-3e8 nohead lt 0;
                 set arrow 3 from 3e8,-3e8 to 3e8,3e8 nohead lt 0;
                 set arrow 4 from -3e8,-3e8 to -3e8,3e8 nohead lt 0' \
		     '' \
		     'linespoints' &
	done
	wait
) & done
wait
