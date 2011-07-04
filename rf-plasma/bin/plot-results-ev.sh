#!/bin/sh -e

usage () {
    echo "plot p(N) (m=M) for each p in PARAMS and foreach result in RESULTS"
    echo "usage: `basename $0` M RESULTS [PARAMS]"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    usage 1>&2
    exit 1
fi

M="$1"
RESULTS="$2"
PARAMS="$3"

if [ -z "$PARAMS" ]; then
    PARAMS="ne vx vz Et Ex Hy Ez"
fi

for p in $PARAMS; do (
    for r in $RESULTS; do (
	parse-results.sh 'if(m=='"$M"') {print N "\t" '"$p"';}' < $r/results.txt \
	    | plt -t "$p(N)" "$r/$p(m$M).png"
    ) & done
    wait
) & done
wait
