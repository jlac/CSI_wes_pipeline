#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my @line=();
my $keyfile = $ARGV[0];
my $samplefile = $ARGV[1];
my $headerfile = $ARGV[2];
my $number='';

open C, ">$samplefile";
open D, ">$headerfile";

open H, "<$keyfile";
while (<H>){
	chomp;
  	last if m/^$/;
  	@line = split;
	if ($line[0] ne 'Exome_ID') {
		$number = ((split 'ATCH', $line[4])[1]);
		if ($number < 29) {
			print C "$line[0]\n";
			print D "$line[0]\t$line[1]\n";
		}
		else {
			print C "$line[1]\n";
			print D "$line[1]\t$line[1]\n";
		}
	}
}