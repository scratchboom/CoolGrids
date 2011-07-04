#!/bin/sh -e

usage () {
    echo "save results to RESULT folder"
    echo "usage: `basename $0` RESULT"
}

if [ -z "$1" ]; then
    usage 1>&2
    exit 1
fi

RESULT="$1"

mkdir $RESULT
cp bconst bconist parskhgs $RESULT/
mv results.txt $RESULT/
