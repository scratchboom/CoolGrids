#/bin/sh -e
log="$1"
png="$2"
extra="$3"
logcomplex.sh EVAL 4 < $log | sort -u | cut -f 1,2 | plt -t "Eigen Values" "$png" "$extra" "" points
