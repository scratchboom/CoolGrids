#!/bin/sh -e

usage() {
    echo "find PARAM in STDIN and then find '(x, y)' pair VALUES_COUNT times and repeat all this procedure"
    echo "usage: `basename $0` PARAM VALUES_COUNT"
}

if [ -z "$1" ]; then
    usage 1>&2
    exit
fi

param="$1"
count="$2"

perl -e '
my $cnt = 0;
while (<>) {
  chomp;
  
  if ($cnt == 0) {
    if ($_ =~ /\b\Q'$param'\E\b/) {
      $cnt = '$count';
    }
  }
  else {
    while (m/\(([0-9.DE+-]*), ?([0-9.DE+-]*)\)/g) {
      print $1 . "\t" . $2 . "\t" . ('$count' - $cnt + 1) . "\n";
      $cnt--;
      if ($cnt == 0) {
        #print "\n";
        last;
      }
    }
  }
}
' \
| tr 'D' 'E'
