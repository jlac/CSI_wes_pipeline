#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my @line = ();
my $inbed = $ARGV[0];
my $outbed = $ARGV[1];

open G, ">$outbed";

open (H, $inbed);
while (<H>) {
	chomp;
  	last if m/^$/;
	@line = split;
	if ($line[0] =~ m'#') {
		print G "$line[0]\t$line[1]\t$line[2]\t$line[3]\tCNVtype\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\n";
	}
	else {
		print G "$line[0]\t$line[1]\t$line[2]\t$line[3]\t";
		if ($line[3]<2) {
			print G "DEL\t";
		}
		else {
			print G "DUP\t";
		}
		print G "$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\n";
	}
}