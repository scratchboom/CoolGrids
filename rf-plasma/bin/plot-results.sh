#!/bin/sh -e

usage () {
    echo "plot p(m) (N=N) for each p in PARAMS and foreach result in RESULTS"
    echo "usage: `basename $0` N RESULTS [PARAMS] [gnuplot-cols.tcl optional arguments]"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    usage 1>&2
    exit 1
fi

N="$1"
RESULTS="$2"
PARAMS="$3"
shift || true
shift || true
shift || true

if [ -z "$PARAMS" ]; then
    PARAMS="ne vx vz Et Ex Hy Ez"
fi

for p in $PARAMS; do (
    for r in $RESULTS; do (
	parse-results.sh 'if(N=='"$N"') {print m "\t" '"$p"';} else if(N>'"$N"') {exit}' < $r/results.txt \
	    | tee "$r/$p(N$N).txt" | plt -t "$p(m)" "$r/$p(N$N).png" "$@" 
    ) & done
    wait
) & done
wait
