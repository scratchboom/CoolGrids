#!/bin/sh -e

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "usage: `basename $0` OUTPUT_FILE EXTRA [[-t PLOT_TITLE] PLOT_FORMAT...]"
    exit 1
fi

OUTFILE=$1
EXTRA=$2
shift
shift

INPUT=
for f in *.out; do
    if echo $f | egrep 'N[0-9]+' 2>&1 >/dev/null; then
		TIME=$(echo $f | perl -e '<> =~ /N0*([0-9]+)/; $N=$1; %h = map {split /\t/, $_} split("\n", `cat times.out`); print $h{$N};')
		N=$(echo $f | perl -e '<> =~ /N0*([0-9]+)/; print $1')
		while [ -n "$1" ]; do
			TITLE=`echo $f | egrep -o 'N[0-9]+'`
			while [ -n "$1" ] && echo "$1" | grep '^-' >/dev/null 2>&1; do
				case $1 in
					-t)
						TITLE="$2"
						shift
						;;
				esac
				shift
			done
			if [ -z "$1" ]; then echo "bad format" >&2; exit 1; fi
			FMT=`echo $1 | sed "s/@TIME\b/$TIME/g" | sed "s/@N\b/$N/g"`
			INPUT="$INPUT $f $FMT $TITLE"
			shift
		done
    fi
done

if [ -z "$INPUT" ]; then
	echo "No input files *.out" >&2
	exit 1
fi

gnuplot-cols.tcl "$INPUT" "$OUTFILE" "$EXTRA"
