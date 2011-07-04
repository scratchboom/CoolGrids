#!/bin/sh -e

if [ -z "$1" ]; then
    echo "usage: `basename $0` FORTRAN_OUTPUT_FILE" >&2
fi

CODE=""
while read n; do
    CODE="$CODE"' if (N == '"$n"') { print N "\t" T } '
done

parse-results.sh ' ' "$CODE" < $1
