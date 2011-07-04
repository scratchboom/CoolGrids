#!/bin/sh -e

usage() {
    echo "parse fortran log file, filter and format output for given variables"
    echo "usage: `basename $0` MODE [VARIABLE]..."
}

all() {
    #tr [:space:] ' ' | \
    grep -o "$ptrn" | sed 's/ //g' | sed 's/=/=	/g'
}

ptrn=""

add_pattern() {
    if [ -z "$ptrn" ]; then
	ptrn="\\(\\b$1 *= *[0-9.edED+-]*\\b\\)"
    else
	ptrn="$ptrn\\|\\(\\b$1 *= *[0-9.edED+-]*\\b\\)"
    fi
}

if [ -z "$1" ]; then
    usage 1>&2
    exit 1
fi

mode=$1
shift

for p in "$@"; do
    add_pattern $1
    shift
done

case $mode in

    all)
	all
	;;

    di*)
	di=`echo $mode | sed s/^di//g`
	for p in N m; do add_pattern $p; done
	all | perl -e '
use strict;
my $n;
my $m;
my $N;
my $v; 
my $flag = 0;
my $head = "m";
my $headprinted = 0;
my $str = "";

sub doprint {
  if ($str ne "") {
    if (!$headprinted) {
      print $head."\n";
      $headprinted = 1;
    }
    print "$m\t$str\n";
    $str="";
  }
} 

while (<>) {
  chomp;
  ($n,$v) = split /=	/;
  if ($n eq "N") {
    $N = $v;
    if ($flag) {
      doprint;
      exit 0;
    }
    if ($N == '$di') { $flag = 1; }
  }
  elsif ($flag) {
    if ($n eq "m") { doprint; $m = $v; } 
    else {
      if (!$headprinted) {
        $head = $head . "\t$n";
      }
      $str = $str . "\t$v";
    }
  }
}
doprint;
'
	;;

    ev*)
	ev=`echo $mode | sed s/^ev//g`
	for p in N m; do add_pattern $p; done
	all | perl -e '
use strict;
my $head = "N";
my $headprinted = 0;
my $str = "";
my $n;
my $m;
my $N;
my $v; 
sub doprint {
  if ($str ne "") {
    if (!$headprinted) {
      print $head."\n";
      $headprinted = 1;
    }
    print $N."\t".$str."\n";
    $str="";
  }
} 

while (<>) {
  chomp;
  ($n,$v) = split /=	/;
  if ($n eq "N") { doprint; $N = $v; }
  elsif ($n eq "m") { doprint; $m = $v;}
  elsif ($m == '$ev') {
    if (!$headprinted) {
      $head = $head . "\t$n";
    }
    $str = $str . "\t$v";
  }
}
doprint;
'
	;;

    *)
	echo "unknown mode: $mode" 1>&2
	;;

esac
