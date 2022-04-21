#!/usr/bin/perl -w
use strict;

my $ini=shift or die "$0 [lachesis.ini]\n";
my $lachesis="/public/agis/ruanjue_group/wuzhichao/software/LACHESIS/src/bin/Lachesis";
#"/public/agis/ruanjue_group/lialun/sw/LACHESIS/bin/Lachesis";

my %param;
open INI,"$ini" or die $!;
while(<INI>){
	my @e=split /\s+/,$_;
	next unless @e>2;
	$param{$e[0]}=$e[2];
}
close INI;
exit if( -f "$param{'OUTPUT_DIR'}/REPORT.txt" );
my $cmd="$lachesis $ini";
print "$cmd\n";
system("$cmd");

$cmd="rm -r $param{'OUTPUT_DIR'}/cached_data";
system("$cmd");


