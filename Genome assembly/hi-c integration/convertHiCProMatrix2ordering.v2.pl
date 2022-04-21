#!/usr/bin/perl -w
use strict;
die "perl $0 [LACHESIS.main.dir] [groupN] [binsize] [out] [hicPro.bed] [hicPro.matrix] [hicPro.matrix2]\n" unless @ARGV>=6;

my $la=shift @ARGV;
my $Ng=shift @ARGV;
$Ng--;
my @chr=(0..$Ng);
my $binsize=shift @ARGV;
my $out=shift @ARGV;
my $bed=shift @ARGV;
my @matrix=@ARGV;
my %order;
my %pass;
my @genomeCTGorder;
for my $chr (@chr){
	my $n=0;
	my $file="$la/group$chr.ordering"; 
	#contig_ID(local)       contig_name     contig_rc       orientation_Q_score     gap_size_after_contig
	#33      ctg155_pilon_pilon      1       299.42  .
	print "$file\n";
	open GG,"$file" or die $!;
	while(<GG>){
		next if /^#/;
		my @e=split;
		$pass{$e[1]}=1;
		$order{$chr}{$n}{'ctg'}=$e[1];
		$order{$chr}{$n}{'r'}=$e[2];
		$n++;
	}
	close GG;
}

print "Reading $bed, bin with size smaller than 0.6x$binsize will be removed\n";
#bed file, HiCPro  defined bin
#contigName start end bin_id
#ctg1000_pilon   0       19340   1
my %bin_coord;
open BED,"$bed" or die $!;
while(<BED>){
	my @e=split;
	next if abs ($e[2]-$e[1])<$binsize*0.6;
	#$bin_coord{$e[0]}{$e[1]}{'end'}=$e[2];
	$bin_coord{$e[0]}{$e[1]}{'bin'}=$e[3];
}
close BED;

#make chrom breaks(chrom length
#convert bin ID
#defined bin coordinates
#bin length of chromosomes composed by each contigs
my %new_bin_coord;
#my %raw_bin;
open CHRB,">$out.chrom.breaks"or die $!;
my $x=0;
foreach my $chr (@chr){
	print CHRB "$x\n";
	foreach my $n (sort {$a<=>$b} keys %{$order{$chr}}){
		my $ctg=$order{$chr}{$n}{'ctg'};
		my $reverse=$order{$chr}{$n}{'r'};
		my @bin_start=sort {$a<=>$b} keys %{$bin_coord{$ctg}};
		@bin_start=reverse @bin_start if $reverse;
		foreach my $bin_start (@bin_start){
			my $raw_bin=$bin_coord{$ctg}{$bin_start}{'bin'};
			#$raw_bin{$x}=$raw_bin; #keys in %raw_bin are new bins
			$new_bin_coord{$raw_bin}=$x; #keys in %new_bin_coord are raw bins
			$x+=0.5;
		}
	}
}
print CHRB  "$x\n";
close CHRB;
%bin_coord=();
exit unless keys %new_bin_coord;

#matrix contact=========
#1       1       827
my %contact;
foreach my $m (@matrix){
	open MM,"$m" or die $!;
	while(<MM>){
		my @e=split;
		$contact{$e[0]}{$e[1]}+=$e[2];
	}
	close MM;
	print "Finished reading $m\n";
}

open HM,">$out.heatmap.txt" or die $!;
print HM "X\tY\tZ\n";
my @bins=sort keys %new_bin_coord;
foreach my $x_raw (@bins){
	foreach my $y_raw (@bins){
		next unless (defined $new_bin_coord{$y_raw} && defined $new_bin_coord{$x_raw} );
		next if (! defined $contact{$x_raw}{$y_raw} && ! defined $contact{$y_raw}{$x_raw} );
		my $contact= defined $contact{$x_raw}{$y_raw} ? $contact{$x_raw}{$y_raw} : $contact{$y_raw}{$x_raw};
		my $y=$new_bin_coord{$y_raw};
		my $x=$new_bin_coord{$x_raw};
		print HM "$x\t$y\t$contact\n";
		#print HM "$y\t$x\t$contact\n" if $x ne $y;;
	}
}
close HM;

