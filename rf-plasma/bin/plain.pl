#!/usr/bin/perl

while (<>) {
    chomp;
    if ($_ =~ /^\s*$/) { print "\n"; next; }
    foreach $v (split /\s+/) {
	if ($v eq "") { next; }
	print "$v\n";
    }
}
