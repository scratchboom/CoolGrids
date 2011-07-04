#!perl

use strict;

my $dir = $ARGV[0];

# exec "mkdir $dir";

my $l = -2;
my $t;
my $z;
my @u;

while (my $s = <>) {
    chomp $s;
    if ($s =~ /^ *$/) {
	$l = 0;
    }
    elsif ($l>0) {
	if (($l-1)%3==0) { ($u[1], $u[2], $u[3]) = split(/ +/, $s); }
	if (($l-1)%3==1) { ($u[4], $u[5], $u[6]) = split(/ +/, $s); }
	if (($l-1)%3==2) { ($u[7]) = split(/ +/, $s); }
	if (($l-1)%3==2) {
	    
	}
	$l++;
    }
    elsif ($l==0) {
	($t, $z) = split(/ +/, $s);
	$l++;
	print "$t $z";
    }
}
