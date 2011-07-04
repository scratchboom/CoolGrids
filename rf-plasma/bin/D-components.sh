#!/bin/sh -e

usage() {
    echo "create "
    echo "usage: `basename $0 OUTPUTDIR`"
}

if [ -z "$1" ]; then
    usage 1>&2
    exit
fi

DIR="$1"

mkdir $DIR
logreal.sh D 44 | awk '
BEGIN{
  /*8-13,15-20,22-27,29-34,36-41,43-48*/
  name[8] = "1.1";
  name[16] = "2.2";
  name[17] = "2.3";
  name[23] = "3.2";
  name[24] = "3.3";
  name[29] = "4.1";
  name[30] = "4.2";
  name[31] = "4.3";
  name[32] = "4.4";
  name[20] = "2.6";
  name[26] = "3.5";
  name[38] = "5.3";
  name[44] = "6.2";
}
name[$2]!~/$^/{
  print name[$2] "\t" $1
}
' | perl -e '
use strict;
my $N = 1;
my $m = 0;
my %files=();
while (<>) {
  chomp;
  my ($name, $value) = split /\t/;
  if ($name eq "1.1") { $m++; }
  if ($m > 1000) { $N++; $m = 1; }
  if ($m < 150 or $m > 160 or $m % 2 == 1) { next; }
  if (not defined $files{"$m.$name"}) {
    open ($files{"$m.$name"}, ">", "'"$DIR"'/$m.$name");
  }
  my $fh = $files{"$m.$name"};
  print $fh "$N\t$value\n";
}
foreach my $k (keys %files) {
  close $files{$k};
}
'
