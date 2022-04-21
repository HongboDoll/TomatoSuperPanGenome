#!/usr/bin/perl -w
use strict;

die "perl [nchr] [restricion.seq] [ref.fa] [bam.dir] [bam] [bam] ... 

nchr:\tNumber of chromosomes
restricion.seq:\tRestriction site, mobI:GATC, HindIII:I do not remmber
ref.fa:\tgenom.fa to be scaffolding 
bam.dir:\tdirector of the bam files
bam:\tbam files, multi bam files supportive
" unless @ARGV>=5;

my $nchr=shift;
my $res=shift;
my $ref=shift;
my $bam_dir=shift;
my @bam=@ARGV;

sub print_baseline{
	my $out=shift;
	print qq{
SPECIES = B45
OUTPUT_DIR = Lachesis/$out

DRAFT_ASSEMBLY_FASTA = $ref 
SAM_DIR = $bam_dir
SAM_FILES = };
print join(" ",@bam),"\n";
print qq{RE_SITE_SEQ = $res
USE_REFERENCE = 0
SIM_BIN_SIZE = 0
REF_ASSEMBLY_FASTA = test_case/hg19/Homo_sapiens_assembly19.fasta

BLAST_FILE_HEAD = test_case/draft_assembly/assembly

DO_CLUSTERING = 1
DO_ORDERING   = 1
DO_REPORTING  = 1

OVERWRITE_GLM = 0
OVERWRITE_CLMS = 0

CLUSTER_N = $nchr
CLUSTER_CONTIGS_WITH_CENS = -1
};
}

#random_value(min,max,number_of_value);
my $num_test=60000;
my @CLUSTER_MIN_RE_SITES=22; #3307 is over #1-4988
my @CLUSTER_MAX_LINK_DENSITY=2; #0.0008-28.9
my @CLUSTER_NONINFORMATIVE_RATIO=2; #1.0004-24.3
my @ORDER_MIN_N_RES_IN_TRUNK=10; #1-5131
my @ORDER_MIN_N_RES_IN_SHREDS=10; #1-6489

my $out=0;
my $dir=1;
my $while1=0;
my %pass_param;
pass_param("01.pass.tested.param.txt") if -f "01.pass.tested.param.txt";
while($out<$num_test){
	my $p1=int random_prarm(@CLUSTER_MIN_RE_SITES);
	my $p2=random_prarm(@CLUSTER_MAX_LINK_DENSITY);
	my $p3=random_prarm(@CLUSTER_NONINFORMATIVE_RATIO);
	my $p4=int random_prarm(@ORDER_MIN_N_RES_IN_TRUNK);
	my $p5=int random_prarm(@ORDER_MIN_N_RES_IN_SHREDS);
	my $p=join("-",($p1,$p2,$p3,$p4,$p5));
	$while1++;
	if ($while1==100){
		print STDERR "Random selected 100 times same parameter, Exit\n";last;
	}
	next if defined $pass_param{$p};
	$pass_param{$p}=1;
	$out++;
	print STDERR "Param$out\n";
	$dir++ if ($out/1000 > $dir);
	system("mkdir run$dir run$dir/Lachesis") if ! -d "run$dir";
	open OUT,">run$dir/Lachesis.greedy.$out.ini" or die $!;
	select OUT;
	print_baseline("run$out");
	print "CLUSTER_MIN_RE_SITES = $p1
CLUSTER_MAX_LINK_DENSITY = $p2
CLUSTER_NONINFORMATIVE_RATIO = $p3
CLUSTER_DRAW_HEATMAP = 0
CLUSTER_DRAW_DOTPLOT = 0
ORDER_MIN_N_RES_IN_TRUNK = $p4
ORDER_MIN_N_RES_IN_SHREDS = $p5\n";
	print "ORDER_DRAW_DOTPLOTS = 0\nREPORT_EXCLUDED_GROUPS = -1\nREPORT_QUALITY_FILTER = 1\nREPORT_DRAW_HEATMAP = 0";
	close OUT;
	$while1=0;
	
}

sub pass_param{
	my $par_file=shift;
	open PM,"$par_file" or die $!;
	while(<PM>){
		next if /^\D/;
		my @e=split;
		my $p=join("-",@e[1..5]);
		$pass_param{$p}=1;
	}
	close PM;
}

sub random_prarm{
	my ($min,$max)=@_;
	my $r=rand($max);
	$r=sprintf "%.4f", $r if $r=~/\./;
	$r+=$min if $r<$min;
	return $r;
}
