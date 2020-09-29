#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT
$dir = $ARGV[0];
my @AllFiles = ();
my @files = ();

#open directory (SWITCH CHROMOSOME ARM HERE)
opendir $dir, "." or die "couldn't open directory\n";
@AllFiles = readdir(DIR);
closedir DIR;

#find correct input files
for ($i = 0; $i < @AllFiles; $i++){
  if ($AllFiles[$i] =~ m/divergence/){
    push @DistFiles, $AllFiles[$i];
  }
}

my $infile = $ARGV[0];
my $fixed = $ARGV[0] . '.tmp';
my $filtered = $ARGV[1];

open C, ">$filtered";

my @line = ();

my $cmd = '';
$cmd = 'awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' ' . $infile . ' > ' . $fixed;
system($cmd);

open G, "<$fixed";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split "\t", $_;
	if ($line[0] =~ m'AnnotSV') {
		print C "$_\n";
	}
	else {
		next if ($line[41] eq '-');
		next if ($line[7] eq 'full');
		next if (($line[1] ne 'X') && ($line[42] eq $line[5]));
#		if ($line[28] < 0.01) {
#			if ($line[36] < 0.01) {
				if (($line[43] eq '-') || ($line[43] > 0.75)) {
					if ((($line[5] == 0) && (($line[50] < 1) && ($line[51] < 1))) || (($line[5] == 1) && ($line[51] < 1)) || ($line[5] > 1)) {
#						if ((($line[6] eq 'DUP') && ($line[32] < 0.01)) || (($line[6] eq 'DEL') && ($line[34] < 0.01))) {
#							if ($line[39] < 0.01) {
								if ($line[78] > 3) {
#								if (($line[73] eq 'yes') || ($line[74] eq 'yes') || (($line[6] eq 'DEL') && (($line[69] ne '-')||($line[69] > 0.8)))) {#|| (($line[6] eq 'DEL') && ($line[70] > 1)) || (($line[6] eq 'DUP') && ($line[71] > 1)) || ($line[72] > 1) || ($line[66] > 90)  || ($line[61] != '-')  || (($line[6] eq 'DEL') && ($line[59] eq '3'))  || (($line[6] eq 'DUP') && ($line[60] eq '3')) || ($line[78] eq '4') || ($line[78] eq '5')) {
									print C "$_\n";
								}
#							}
#						}
#					}
#				}
			}
		}
	}
}

$cmd = 'rm ' . $fixed;
system($cmd);
