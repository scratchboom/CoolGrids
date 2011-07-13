#!/bin/sh -e
LEN="$1"
MK="$2"
EXT=png
EXTRA="$3 set terminal $EXT size 500,350; set key off; set logscale y; set xrange [0:$LEN];"
CS="$4"
ZMOVE='(($1/'"$MK"'*'"$LEN"')-@TIME*'"$CS"')'
Z='($1/'"$MK"'*'"$LEN"')'
mplot.sh ne.$EXT "$EXTRA"' set ylabel "Концентрация Ne, м^-3";' "$Z"':($2)'
mplot.sh Et.$EXT "$EXTRA"' set ylabel "Температура Te, эВ";' "$Z"':($5)'
