#!/bin/sh -e

for n in $*; do
	ac_slice N=$n < /dev/null
	last=$n
done > out.$$.awk < /dev/null

ac_slice e=`expr $last + 1` >> out.$$.awk
parse-results.sh "`cat out.$$.awk`" | perl -e 'use strict;
my %fh = ();
while (<STDIN>) {
  chomp;
  my @fld = split /\t/;
  my $key = shift @fld;
  if (!defined($fh{$key})) {
    my $fn = $ARGV[0] . sprintf(".N%06d.out", $key);
    open my $tmp, ">$fn";
    $fh{$key} = $tmp;
  }
  my $tmp = $fh{$key};
  print $tmp join("\t", @fld)."\n";
}
foreach (values %fh) { close $_; }
' `date +%Y%m%d-%H%M%S`

rm -f out.$$.awk
