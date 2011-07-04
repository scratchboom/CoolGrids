#!/bin/sh -e
MP=`dirname $0`/mplot.sh
$MP field-form.png 'set xrange [-0.3:0.1]; set yrange[-1e7:1e7];' '(($1/2000*6)-@TIME*3e8):($6)'
$MP ne.png ' ' 1:2
