#!/usr/bin/perl -w
use strict;
die "perl $0 [chrNum] [path/to/REPORT.txt] [path/to/REPORT.txt] ... 
\n" if (@ARGV<2 or $ARGV[0]=~/\D/);
my $chr_num=shift;

my @header=("CLUSTER_MIN_RE_SITES","CLUSTER_MAX_LINK_DENSITY","CLUSTER_NONINFORMATIVE_RATIO","ORDER_MIN_N_RES_IN_TRUNK","ORDER_MIN_N_RES_IN_SHREDS","Nchr","Nclust","Lclust","Norder","Lorder");
print "Run\t",join("\t",@header),"\tScore";
my @scaf=reverse (1..$chr_num);
#my $lastScaf=(sort {$a<=>$b} @scaf)[-1];
foreach (@scaf){
	print "\tScaffold_$_";
}
print "\n";
for(@ARGV){
	my @x=(split /\//,$_);
	pop @x;
	my $run=pop @x;
	print "$run";
	readReport($_);
}

sub meanQscore{
	my $reportfile=shift;
	my $dir=$reportfile;
	$dir=~s/REPORT.txt//;
#contig_ID(local)       contig_name     contig_rc       orientation_Q_score     gap_size_after_contig
#27      ctg163_pilon    0       917.039 .
#50      ctg276_pilon    1       1111.42 .
	my ($sum,$n);
	for(00..$chr_num-1){
		my $odfile="$dir/main_results/group$_.ordering";
		open OD,"$odfile" or die $!;
		while(<OD>){
			next if /^#/;
			my @e=split;
			next unless @e>=4;
			$sum+=$e[3];
			$n++;
		}
		close OD;
	}
	my $mean=sprintf "%.1f",$sum/$n;
	print "\t$mean";
	return 1;
}
sub readReport{
	my $file=shift;
	my %param;
	open RP,"$file" or die $!;
	while(<RP>){
		next if /^#/;
		next unless /\S+/;
		my @e=split /[\s=]+/,$_;
		if($e[0] eq "CLUSTER_MIN_RE_SITES" or $e[0] eq "CLUSTER_MAX_LINK_DENSITY" or $e[0] eq "CLUSTER_NONINFORMATIVE_RATIO" or $e[0] eq "ORDER_MIN_N_RES_IN_TRUNK" or $e[0] eq "ORDER_MIN_N_RES_IN_SHREDS"){
			$param{$e[0]}=$e[1];
		}elsif(/N\sclusters.*\s(\d+)\n/){
			$param{'Nchr'}=$1;
		}elsif(/Number of contigs in clusters:\s+(\d+)\s+/){
			$param{'Nclust'}=$1;
		}elsif(/Length of contigs in clusters:\s+(\d+)\s+/){
			$param{'Lclust'}=$1;
		}elsif(/Number of contigs in orderings:\s+(\d+)\s/){
			$param{'Norder'}=$1;
		}elsif(/Length of contigs in orderings:\s+(\d+)\s/){
			$param{'Lorder'}=$1;
		}elsif(/\|\s+([0-9]+)\s+\|\s+(\d+)\s+\|\s+(\d+)\s+\|/){
			#|     58   |       3   |     6705543 |
			my $n=$1;$n++;
			$param{$n}=$3;
		}
	}
	close RP;
	for my $i (@header){
		print "\t$param{$i}";
	}
	meanQscore($file);
	print ";";
	foreach my $n (@scaf){
		if(defined $param{$n}){
			print "\t$param{$n}";
		}else{
			print "\tNA";
		}
	}
	print "\n";
	return 1;
}
#N clusters (derived):   26
#N non-singleton clusters:       26
#N orderings found:      26
