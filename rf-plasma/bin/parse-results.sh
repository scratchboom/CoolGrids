#!/bin/sh -e

usage() {
    echo "parse fortran output file (results.txt)"
    echo "usage: `basename $0` AWK_CODE_FOREACH_DATAPOINT [AWK_CODE_FOREACH_DATALINE] [BEGIN_CODE] [END_CODE]"
}

getparams='
  ne = $1;
  getline; vx = $1;
  getline; vz = $1;
  getline; Et = $1;
  getline; Ex = $1;
  getline; Hy = $1;
  getline; Ez = $1;
'

getlineparams='
  getline; N = $1;
  getline; T = $1;
  getline; Z = $1;
'

while [ ! -z "$1" ] && (echo "$1" | egrep '^-' >/dev/null 2>&1); do 
    case $1 in
		-p)
			getparams=`echo "$2" | perl -e '$_ = <>; chomp; s/,/ = \\$1; getline; /g; print; print " = \\$1;";'`
			shift
			;;
		-l)
			getlineparams=`echo "$2" | perl -e '$_ = <>; chomp; foreach (split /,/) { print "getline; $_ = \\$1; " }'`
			shift
			;;
		-A)
			tmp=m,`perl -e 'print join ",", map {"A$_"} map {$a=$_;map {"$a$_"} 1..6} 1..6;'`,B4
			getparams=`echo "$tmp" | perl -e '$_ = <>; chomp; s/,/ = \\$1; getline; /g; print; print " = \\$1;";'`
			;;
    esac
    shift
done

if [ -z "$1" ]; then
    usage 1>&2
    exit 1
fi

code="$1"
Ncode="$2"
Bcode="$3"
Ecode="$4"

plain.pl | awk '
BEGIN {
  '"$Bcode"'
}
/^$/{
  m=0;
  '"$getlineparams"'
  getline;
  '"$Ncode"'
}
{
  m++;
  '"$getparams"'
  '"$code"'
}
END {
  '"$Ecode"'
}
'
