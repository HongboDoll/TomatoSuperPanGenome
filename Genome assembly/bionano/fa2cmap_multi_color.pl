# $Id: fa2cmap_multi_color.pl 8271 2018-12-11 21:57:23Z twang $
#!/usr/bin/perl -w
#################################################################################
# File: fa2cmap_multi_color.pl                                                  #
# Date: 10/06/2016                                                              #
# Purpose: Transform fasta file to BioNano cmap file format (Color-aware)       #
#                                                                               #
# Author: Xiang Zhou, Computational Biologist                                   #
# Email : xzhou@bionanogenomics.com                                             #
# Affiliation: Research Department, BioNano Genomics Inc.                       #
#                                                                               #
# Usage:                                                                        #
#   fa2cmap_multi_color.pl [options] <Args>                                     #
# Options:                                                                      #
#   -h : This help message                                                      #
#   -v : Print program version information                                      #
#   -i : Input fasta file (Required)                                            #
#   -k : Input key file (If provided, use as contigIDs rather than sequentially)#
#   -o : Output folder (Default: the same as the input file)                    #
#   -e : Name or sequence of the enzyme (or a combination of name and sequence  #
#        in the format of name:sequence) followed by channel # (Can be multiple)#
#   -m : Filter: Minimum sites (Integer, default: 0)                            #
#   -M : Filter: Minimum size (kb) (Integer, default: 0)                        #
#   -B : Output BED file of Nbase gaps (Default: OFF)                           #
#   -g : Minimum N gap size (bp) when generating Nbase gap file (Default: 1000) #
#   -W : For web use only, and must be the first option (Default: OFF)          #
#   -S : Add an additional column, Strand, to the cmap file (Default: OFF)      #
#                                                                               #
# NOTE: CMAP index is 1-based, and is color-aware.                              #
#################################################################################

use strict;
use warnings;
use POSIX;
use File::Spec;
use File::Basename;
use File::Spec::Functions;
#use List::MoreUtils qw(uniq);
#use Getopt::Std;
use File::Path qw(make_path);
use Getopt::Long qw(:config no_ignore_case);
#use File::Slurp qw(edit_file_lines);
#$Getopt::Std::STANDARD_HELP_VERSION = 1;

sub Init;
sub Usage;
sub Version;
sub Find_enzymes;
sub Print_cmap_header;
sub Generate_cmap;
sub Get_and_set_enzyme_sequence;
sub is_enzyme_name;
sub is_nt_sequence;
sub is_enzyme_name_and_sequence;
sub Uniq;
sub Get_distance;
sub Get_N50;
sub Get_N90;
sub StdDev;
sub getNGaps;
sub printNGaps;

my %enzyme = (
	"BspQI" => "GCTCTTC",
	"BbvCI" => "CCTCAGC",
	"BsmI"  => "GAATGC",
	"BsrDI" => "GCAATG",
	"bseCI" => "ATCGAT",
	"DLE1" => "CTTAAG",
	"BssSI" => "CACGAG"
	
	# You can add more enzymes here ...
	
);
my (@enzymes, @enzyme_sequences, @color, @color_uniq);
my ($min_labels, $min_length, $minNGapSize) = (0, 0, 1000);
my ($FASTA, $CMAP, $KEY, $BED, $SUMMARY, $STATUS, $ERROR, $KEY_IN);

# Autoflush for files
my $old_fh = select($STATUS);
local $| = 1;
select($old_fh);

my $num_args = $#ARGV;
my $argvs = "";
for(my $i = 0; $i <= $num_args; $i++){
	$argvs .= " $ARGV[$i]";
}
my ($help, $version, $web, $strand, $keyfile, $input, $output, $nbase_bed, $bedfile);
my ($strand_header, $strand_header_type) = ("", "");
my (@loci, @loci_tmp);
my (%color, %color_to_enzyme_sequence, %strand);
my ($count, $num_cmaps, $total_nicks) = (0, 0, 0);
my @num_nicks;  # For each channel
my @num_nicks_global;  # For each channel
my @density;  # For each channel
my @total_label_distance;  # For each channel
my @label_distances_by_channel;  # Two-dimensional array containing label distances of each channel
my $density_global;
my ($total_length_original, $total_length_filtered, $total_label_distance_filtered) = (0, 0, 0);
my (@length_filtered, @label_distance_filtered);
my $nGaps = {};  # hashref variable for NGap begin/end with minNGapSize - Zhanyang

my ($command, $seq, $seq_length, $tmp);
my ($fastaHeader, $fastaHeader_previous);
my ($A, $C, $G, $T, $N, $global_N, $global_GC, $global_ACGT, $global_ACGTN);
my ($N_percentage, $GC_percentage, $global_N_percentage, $global_GC_percentage);

Init();
my ($filename, $filename_key, $filename_bed, $filename_summary, $filename_status, $filename_error) = ($input) x 6;

for(my $i = 0; $i < @color_uniq; $i++){
	@num_nicks_global[$color_uniq[$i]] = 0;
}
for(my $i = 0; $i < @enzymes; $i++){
	$tmp .= "_$enzymes[$i]";
}

# If output folder is not defined, use the input folder
if($output){
	#if(substr($output, -1) eq "/"){
	#	$filename = $output . (split("/", $input))[-1];
	#	$filename_key = $output . (split("/", $input))[-1];
	#}
	#else{
	#	$filename = $output . "/" . (split("/", $input))[-1];
	#	$filename_key = $output . "/" . (split("/", $input))[-1];
	#}
	
	# TODO: Convert ".." to absolute path!
	my ($nothing, $out_dir) = fileparse($output."/");
	#mkdir($out_dir);
	make_path($out_dir) unless (-d $out_dir);
	#print "\nOUTPUT_DIR\n", $out_dir, "\n\n";
	
	$filename = basename($filename);
	$filename = catfile($out_dir, $filename);
	$filename_key = basename($filename_key);
	$filename_key = catfile($out_dir, $filename_key);
	$filename_bed = basename($filename_bed);
	$filename_bed = catfile($out_dir, $filename_bed);
	$filename_summary = basename($filename_summary);
	$filename_summary = catfile($out_dir, $filename_summary);
	$filename_status = catfile($out_dir, "status.txt");
	$filename_error = catfile($out_dir, "error.txt");
}
else{
	my ($nothing, $in_dir) = fileparse($input);
	$filename_status = catfile($in_dir, "status.txt");
	$filename_error = catfile($in_dir, "error.txt");
}
$filename =~ s/(\S+)\.\w+$/$1${tmp}_${min_length}kb_${min_labels}labels.cmap/;
$filename_key =~ s/(\S+)\.\w+$/$1${tmp}_${min_length}kb_${min_labels}labels_key.txt/;
$filename_summary =~ s/(\S+)\.\w+$/$1${tmp}_${min_length}kb_${min_labels}labels_summary.txt/;
$filename_bed =~ s/(\S+)\.\w+$/$1_min${minNGapSize}_nbase.bed/;

open($STATUS, ">$filename_status") || die ("ERROR: Can't open $filename_status: $!\n");
printf $STATUS("0%% done\n");
my $num_contigs = 0;
my $ct = 0;
open($FASTA, $input) || die ("ERROR: Can't open $input: $!\n");
while(my $line = <$FASTA>){
	if($ct == 0){
		$ct++;
		#chomp $line;
		#$line =~ s/\r//g;
		#if( $line ne "" && !($line =~ /^[>ACGTNacgtn]/) ){
		if(! ($line =~ /^>/) ){
			open($ERROR, ">$filename_error") || die ("ERROR: Can't open $filename_error: $!\n");
			print ("ERROR: Invalid input file!\n");
			print $ERROR("ERROR: Invalid input file!\n");
			close($ERROR);
			Usage();
		}
	}
	if($line =~ /^>/){
		$num_contigs++;
	}
}
close($FASTA);

my %table;
if( defined($keyfile) ){
	open($KEY_IN, $keyfile) || die ("ERROR: Can't open $keyfile: $!\n");
	while(my $line = <$KEY_IN>){
		if($line =~ /^#/){
			next;
		}
		else{
			chomp $line;
			$line =~ s/\r//g;
			
			my @x = split("\t", $line);
			$table{$x[1]} = $x[0];
		}
	}
	close($KEY_IN);
}

open($KEY, ">$filename_key") || die ("ERROR: Can't open $filename_key: $!\n");
print $KEY("# CMAP = ", File::Spec -> rel2abs($filename), "\n");
print $KEY("# filter: Minimum Sites = $min_labels\n");
print $KEY("# filter: Minimum Size (kb) = $min_length\n");
print $KEY("CompntId\tCompntName\tCompntLength\n");

open($SUMMARY, ">$filename_summary") || die ("ERROR: Can't open $filename_summary: $!\n");
print $SUMMARY("#PARAMS\n");
print $SUMMARY("Version:\t", '$Id: fa2cmap_multi_color.pl 8271 2018-12-11 21:57:23Z twang $', "\n");
print $SUMMARY("Command:\t", "$0", $argvs, "\n");
print $SUMMARY("Input file:\t", File::Spec -> rel2abs($input), "\n");
print $SUMMARY("Output file:\t", File::Spec -> rel2abs($filename), "\n");
print $SUMMARY("ID keyfile:\t", File::Spec -> rel2abs($filename_key), "\n");
if($nbase_bed){
	# Do not re-generate Nbase gap file!
	if(-e $filename_bed){
		$nbase_bed = 0;
	}
	print $SUMMARY("BED file:\t", File::Spec -> rel2abs($filename_bed), "\n");
}
print $SUMMARY("Minimum sites filter:\t$min_labels\n");
print $SUMMARY("Minimum size filter (kbp):\t$min_length\n");
for(my $i = 0; $i < @enzymes; $i++){
	print $SUMMARY("Enzyme ", $i+1, " (Channel $color[$i]):\t$enzymes[$i]($enzyme_sequences[$i])\n");
}
print $SUMMARY("\n");

print $SUMMARY("#RESULT\n");
print $SUMMARY("CMapId\tChannel\tLength(Mb)\tDensity(sites/100kbp)\tGC%\tN%\n");
#print $SUMMARY("------------------------------------------------------------------\n");

Print_cmap_header($filename);

open($FASTA, $input) || die ("ERROR: Can't open $input: $!\n");
while(my $line = <$FASTA>){
	chomp $line;
	$line =~ s/\r//g;
	
	if($line =~ /^>/){
		$fastaHeader_previous = $fastaHeader;
		$fastaHeader = substr($line, 1);
		
		if($count != 0){
			$seq_length = length($seq);
			
			# IF the sequence length is 0, ignore this sequence!!!
			if(!$seq_length){
				$count++;
				next;
			}
			
			$A = ($seq =~ tr/A/A/);
			$C = ($seq =~ tr/C/C/);
			$G = ($seq =~ tr/G/G/);
			$T = ($seq =~ tr/T/T/);
			$N = ($seq =~ tr/N/N/);
			$global_N += $N;
			$global_GC += ($C+$G);
			$global_ACGT += ($A+$C+$G+$T);
			$global_ACGTN += ($A+$C+$G+$T+$N);
			
			$total_length_original += $seq_length;
			
			@loci_tmp = ();
			%color = ();
			for(my $i = 0; $i < @color_uniq; $i++){
				@num_nicks[$color_uniq[$i]] = 0;
			}
			for(my $i = 0; $i < @enzyme_sequences; $i++){
				my @nicks_tmp = Find_enzymes($seq, $enzyme_sequences[$i], $color[$i]);
				push(@loci_tmp, @nicks_tmp);
			}
			
			# Remove duplicated values!
			@loci = Uniq(@loci_tmp);
			
			# Filter by $min_labels and $min_length
			if(scalar(@loci) >= $min_labels && $seq_length >= $min_length * 1000){
				if(%table){
					Generate_cmap($filename, $table{$fastaHeader_previous}, \@loci, $seq_length);
					print $KEY("$table{$fastaHeader_previous}\t$fastaHeader_previous\t", $seq_length, "\n");
				}
				else{
					Generate_cmap($filename, $count, \@loci, $seq_length);
					print $KEY("$count\t$fastaHeader_previous\t", $seq_length, "\n");
				}
				
				my $length_tmp = sprintf("%.3f", $seq_length/1000000);
				
				if($A+$C+$G+$T == 0){
					$GC_percentage = sprintf("%.1f", 0);
				}
				else{
					$GC_percentage = sprintf("%.1f", ($C+$G)/($A+$C+$G+$T)*100);
				}
				if($N == 0){
					$N_percentage = sprintf("%.1f", 0);
				}
				else{
					$N_percentage = sprintf("%.1f", $N/$seq_length*100);
				}
				
				for(my $i = 0; $i < @color_uniq; $i++){
					$density[$count] = sprintf("%.1f", $num_nicks[$color_uniq[$i]]/$seq_length * 100000);
					print $SUMMARY("$count\t$color_uniq[$i]\t$length_tmp\t$density[$count]\t$GC_percentage\t$N_percentage\n");
				}
				
				# NOTE: No duplicate sites!!!
				$total_nicks += @loci;
				$total_length_filtered += $seq_length;
				push(@length_filtered, $seq_length);
			}
			# Count Ngaps for this chromosome (Zhanyang Zhu)
			if($nbase_bed){
				$nGaps->{$count} = getNGaps($seq, $minNGapSize);
			}
		}
		
		# Generate the status bar
		if($count % ceil($num_contigs/100) == 0){
			printf $STATUS("%d%% done\n", $count/$num_contigs*100);
		}
		
		# Reset
		$seq = "";
		$count++;
	}
	else{
		$seq .= uc($line);
	}
}

$seq_length = length($seq);

# IF the sequence length is not 0!!!
if($seq_length){
	$A = ($seq =~ tr/A/A/);
	$C = ($seq =~ tr/C/C/);
	$G = ($seq =~ tr/G/G/);
	$T = ($seq =~ tr/T/T/);
	$N = ($seq =~ tr/N/N/);
	$global_N += $N;
	$global_GC += ($C+$G);
	$global_ACGT += ($A+$C+$G+$T);
	$global_ACGTN += ($A+$C+$G+$T+$N);

	$total_length_original += $seq_length;
	
	@loci_tmp = ();
	%color = ();
	for(my $i = 0; $i < @color_uniq; $i++){
		@num_nicks[$color_uniq[$i]] = 0;
	}
	for(my $i = 0; $i < @enzyme_sequences; $i++){
		my @nicks_tmp = Find_enzymes($seq, $enzyme_sequences[$i], $color[$i]);
		push(@loci_tmp, @nicks_tmp);
	}
	
	# Remove duplicated values!
	@loci = Uniq(@loci_tmp);
	
	if(scalar(@loci) >= $min_labels && $seq_length >= $min_length * 1000){
		if(%table){
			Generate_cmap($filename, $table{$fastaHeader}, \@loci, $seq_length);
			print $KEY("$table{$fastaHeader}\t$fastaHeader\t", $seq_length, "\n");
		}
		else{
			Generate_cmap($filename, $count, \@loci, $seq_length);
			print $KEY("$count\t$fastaHeader\t", $seq_length, "\n");
		}
		
		my $length_tmp = sprintf("%.3f", $seq_length/1000000);
		
		if($A+$C+$G+$T == 0){
			$GC_percentage = sprintf("%.1f", 0);
		}
		else{
			$GC_percentage = sprintf("%.1f", ($C+$G)/($A+$C+$G+$T)*100);
		}
		if($N == 0){
			$N_percentage = sprintf("%.1f", 0);
		}
		else{
			$N_percentage = sprintf("%.1f", $N/$seq_length*100);
		}
		
		for(my $i = 0; $i < @color_uniq; $i++){
			$density[$count] = sprintf("%.1f", $num_nicks[$color_uniq[$i]]/$seq_length * 100000);
			print $SUMMARY("$count\t$color_uniq[$i]\t$length_tmp\t$density[$count]\t$GC_percentage\t$N_percentage\n");
		}
		
		# NOTE: No duplicate sites!!!
		$total_nicks += @loci;
		$total_length_filtered += $seq_length;
		push(@length_filtered, $seq_length);
	}
	
	if($nbase_bed){
		# Count Ngaps for the last chromosome (Zhanyang Zhu)
		$nGaps->{$count} = getNGaps($seq, $minNGapSize);
	}
}
close($FASTA);
close($KEY);

if($nbase_bed){
	printNGaps($nGaps, $filename_bed, $minNGapSize);
}

#$command = "$^X -pi -e 's|N/A|$num_cmaps|' $filename";
#print "Running command:\n$command\n\n";
#system("$command");
#edit_file_lines {s |N/A|$num_cmaps| } $filename;

$command = "$^X -i.bak -p -e \"s/N\\/A/$num_cmaps/\" $filename";
#print "Running command:\n$command\n\n";
system("$command");
unlink("$filename.bak");

#print $SUMMARY("\n================================Summary================================\n");
print $SUMMARY("\n#SUMMARY\n");
print $SUMMARY("Minimum sites filter:\t$min_labels\n");
print $SUMMARY("Minimum size filter (kbp):\t$min_length\n");
print $SUMMARY("Total sequence entries processed:\t$count\n");
print $SUMMARY("Total maps generated:\t$num_cmaps\n");

print $SUMMARY("Total length of the original sequence entries (bp):\t$total_length_original\n");
print $SUMMARY("Total length of the filtered maps (bp):\t$total_length_filtered\n");
print $SUMMARY("Map N50 (bp):\t", Get_N50(\@length_filtered, $total_length_filtered), "\n");
print $SUMMARY("Map N90 (bp):\t", Get_N90(\@length_filtered, $total_length_filtered), "\n");
$global_GC_percentage = sprintf("%.1f", $global_GC/$global_ACGT*100);
$global_N_percentage = sprintf("%.1f", $global_N/$global_ACGTN*100);
print $SUMMARY("Global GC (%):\t$global_GC_percentage\n");
print $SUMMARY("Global N (%):\t$global_N_percentage\n\n");

foreach my $channel ( sort {$a <=> $b} @color_uniq ){
	my $distance_tmp = sprintf("%.1f", Get_N50(\@{$label_distances_by_channel[$channel]}, $total_label_distance[$channel]));
#	print $SUMMARY("Site to site distance N50 for channel $channel\t$distance_tmp\n");
	print $SUMMARY("Channel $channel site to site distance N50 (bp):\t$distance_tmp\n");
}
my $label_distance_global = sprintf("%.1f", Get_N50(\@label_distance_filtered, $total_label_distance_filtered));
#print $SUMMARY("Overall site to site distance N50\t$label_distance_global\n");
print $SUMMARY("\n");

foreach my $channel ( sort {$a <=> $b} @color_uniq ){
	my $distance_stddev_tmp = sprintf("%.1f", StdDev(\@{$label_distances_by_channel[$channel]}));
#	print $SUMMARY("Standard deviation of site to site distance for channel $channel\t$distance_stddev_tmp\n");
	print $SUMMARY("Channel $channel standard deviation of site to site distance (bp):\t$distance_stddev_tmp\n");
}
my $distance_stddev_global = sprintf("%.1f", StdDev(\@label_distance_filtered));
#print $SUMMARY("Standard deviation of overall site to site distance\t$distance_stddev_global\n");
print $SUMMARY("\n");

foreach my $channel ( sort {$a <=> $b} @color_uniq ){
	my $density_tmp;
	my @numLabel_filtered = GetNumLabelByDist(\@{$label_distances_by_channel[$channel]});
	my $density_saphyr;
	my $density_irys;

	if($total_length_filtered){
		$density_tmp = sprintf("%.1f", $num_nicks_global[$channel]/$total_length_filtered * 100000);
		$density_saphyr = sprintf("%.1f", $numLabel_filtered[1]/$total_length_filtered * 100000);
		$density_irys = sprintf("%.1f", $numLabel_filtered[2]/$total_length_filtered * 100000);
	}
	else{
		$density_tmp = sprintf("%.1f", 0);
		$density_saphyr = sprintf("%.1f", 0);
		$density_irys = sprintf("%.1f", 0);
	}
	print $SUMMARY("Channel $channel site density (sites/100kbp):\t$density_tmp\n");
	print $SUMMARY("Channel $channel estimated label density (labels/100kbp) for Saphyr:\t$density_saphyr\n");
	print $SUMMARY("Channel $channel estimated label density (labels/100kbp) for Irys:\t$density_irys\n\n");
}
if($total_length_filtered){
	$density_global = sprintf("%.1f", $total_nicks/$total_length_filtered * 100000);
}
else{
	$density_global = sprintf("%.1f", 0);
}
#print $SUMMARY("Global site frequency (sites/100kb)\t$density_global\n");

#print $SUMMARY("=======================================================================\n");
close($SUMMARY);

sleep(1);
printf $STATUS("%d%% done\n", $count/$num_contigs*100);
close($STATUS);
exit 0;
######################################################################
#                           Subroutines                              #
######################################################################
sub Init{
=COMMENT
	my $opt_string = 'hvBi:o:k:s:e:m:M:g:WS';
	if(!getopts("$opt_string", \%opts)){
		print("ERROR: Invalid parameter(s)!\n");
		Usage();
	}
=cut
	my @tmp;
	my $ret = GetOptions(
		'help|h|?'        => \$help,
		'version|v'       => \$version,
		'input|i=s'       => \$input,
		'keyfile|k=s'     => \$keyfile,
		'output|o=s'      => \$output,
		'enzyme|e=s{1,18}' => \@tmp,
		'min_labels|m:i'  => \$min_labels,
		'min_length|M:i'  => \$min_length,
		'nbaseBED|B'      => \$nbase_bed,
		'minSizeNGap|g:i' => \$minNGapSize,
		'Web|W'           => \$web,
		'Strand|S'        => \$strand,
	);
	
	if(!$ret){
		print("ERROR: Missing or invalid parameter(s)!\n");
		Usage();
	}
	
	Usage() if $help;
	Version() if $version;
	
	if(!$input){
		print("ERROR: Missing parameter(s)!\n");
		Usage();
	}
	
	if(!@tmp){
		print("ERROR: Missing parameter(s)!\n");
		Usage();
	}
	
	# Total number of values must be even
	if(@tmp%2 != 0){
		print("ERROR: There must be even numbers of -e parameters!\n");
		Usage();
	}
	
	for(my $i = 1; $i < @tmp; $i += 2){
		# The values at odd positions must be positive integers
		if($tmp[$i] =~ /^\d+$/){
			$color[($i-1)/2] = $tmp[$i];
		}
		else{
			print("ERROR: Invalid parameter(s)!\n");
			Usage();
		}
		
		# The values at even positions must be either name or sequence
		if(is_enzyme_name_and_sequence($tmp[$i-1])){
			#ZZhu: when name:sequence is provided
			my @name_seq = split(":", $tmp[$i-1]);
			$enzyme_sequences[($i-1)/2] = uc($name_seq[1]);
			$enzymes[($i-1)/2] = uc($name_seq[0]);
		}
		elsif(is_enzyme_name(\%enzyme, $tmp[$i-1])){
			$enzyme_sequences[($i-1)/2] = Get_and_set_enzyme_sequence(\%enzyme, \$tmp[$i-1]);
			$enzymes[($i-1)/2] = uc($tmp[$i-1]);
		}
		elsif(is_nt_sequence($tmp[$i-1])){
			$enzyme_sequences[($i-1)/2] = uc($tmp[$i-1]);
			$enzymes[($i-1)/2] = uc($tmp[$i-1]);
		}
		else{
			print("ERROR: Invalid parameter(s)!\n");
			Usage();
		}
		
		push( @{$color_to_enzyme_sequence{$tmp[$i]}}, $enzyme_sequences[($i-1)/2] );
	}
	
	@color_uniq = sort {$a <=> $b} (Uniq(@color));
	
	if($minNGapSize < 1){
		$minNGapSize = 1;
	}
	
	if($strand){
		($strand_header, $strand_header_type) = ("\tStrand", "\tchar");
	}
	
	# @enzymes:          enzyme names in the same order as the input
	# @enzyme_sequences: enzyme sequences in the same order as the input
	# @color:      colors of each enzyme in the same order as the input
	# @color_uniq: a unique set of colors of all enzymes in ascending order
	
	#print "[@enzymes]", "\n";
	#print "[@enzyme_sequences]", "\n";
	#print "[@color]", "\n";
	#print "[@color_uniq]", "\n";
}

sub Usage{
	if($web){
		exit 1;
	}
	
	print << "EOF";

Usage: $^X $0 [options] <Args>
Options:
  -h : This help message
  -v : Print program version information
  -i : Input fasta file (Required)
  -k : Input key file (If provided, use as contigIDs rather than sequentially)
  -o : Output folder (Default: the same as the input file)
  -e : Name or sequence of the enzyme (or a combination of name and sequence
       in the format of name:sequence) followed by channel # (Can be multiple)
  -m : Filter: Minimum sites (Integer, default: 0)
  -M : Filter: Minimum size (kb) (Integer, default: 0)
  -B : Output BED file of Nbase gaps (Default: OFF)
  -g : Minimum N gap size (bp) when generating Nbase gap file (Default: 1000)
  -W : For web use only, and must be the first option (Default: OFF)
  -S : Add an additional column, Strand, to the cmap file (Default: OFF)

NOTE: CMAP index is 1-based, and is color-aware.
EOF
	exit 1;
}

sub Version{
	# $Id: fa2cmap_multi_color.pl 8271 2018-12-11 21:57:23Z twang $
	my $REV = '$Id: fa2cmap_multi_color.pl 8271 2018-12-11 21:57:23Z twang $';
	($REV) = $REV =~ /\$Id: (.*)\$/;
	if($REV){
		print $REV, "\n";
	}
	else{
		print "No version # was found!", "\n";
	}
	exit 0;
}

sub Find_enzymes{
	my ($seq, $enzyme, $channel) = @_;
	my @result;
	
	# HashTable: $color{ nicking_site }[0] = total_num_of_colors
	# HashTable: $color{ nicking_site }[1 .. total_num_of_colors_at_this_location] = $channel
	# HashTable: $strand{$nicking_site} = "+" or "-"
	
	# Find the enzymes in the forward strand, staring from the first nucleotide!!!
	my $current_loc = index($seq, $enzyme, 0);
	while ($current_loc != -1){
		# Add the current location into @result IFF
		# 1) it is not in @result, or
		# 2) it is in @result, and the current color is not in the color array.
		if( !defined($color{$current_loc+1}[0]) || !($channel ~~ @{$color{$current_loc+1}}[1 .. $color{$current_loc+1}[0]]) ){
			# This is color-aware!
			push(@result, $current_loc+1);
			
			# Increase the color count
			if( !defined($color{$current_loc+1}) ){
				$color{$current_loc+1}[0] = 1;
			}
			else{
				$color{$current_loc+1}[0]++;
			}
			
			# Store the current color on the current nicking site in the color array (unique, no order!):
			# $color{ nicking_site }[1 .. total_num_of_colors_at_this_location]
			$color{$current_loc+1}[ $color{$current_loc+1}[0] ] = $channel;
		}
		
		$strand{$current_loc+1} = "+";
		
		$current_loc = index($seq, $enzyme, $current_loc + 1);
	}

	my $enzyme_rc = reverse($enzyme);
	$enzyme_rc =~ tr/ACGTUN/TGCAAN/;
	
	# Find the rc(enzymes) in the forward strand, staring from the first nucleotide!!!
	$current_loc = index($seq, $enzyme_rc, 0);
	while ($current_loc != -1){
		# Add the current location into @result IFF
		# 1) it is not in @result, or
		# 2) it is in @result, and the current color is not in the color array.
		if( !defined($color{$current_loc+1}[0]) || !($channel ~~ @{$color{$current_loc+1}}[1 .. $color{$current_loc+1}[0]]) ){
			# This is color-aware!
			push(@result, $current_loc+1);
			
			if( !defined($color{$current_loc+1}) ){
				$color{$current_loc+1}[0] = 1;
			}
			else{
				$color{$current_loc+1}[0]++;
			}
			$color{$current_loc+1}[ $color{$current_loc+1}[0] ] = $channel;
		}
		$strand{$current_loc+1} = "-";
		
		$current_loc = index($seq, $enzyme_rc, $current_loc + 1);
	}
	
	# Remove duplicated values!
	return Uniq(@result);
}

sub Print_cmap_header{
	my ($filename) = @_;
	my $OUT;
	my $tmp = "";
	my $color_uniq = scalar(@color_uniq);
	
	open($OUT, ">$filename") || die ("ERROR: Can't open $filename: $!\n");
	
	#for(my $i = 0; $i < @enzyme_sequences; $i++){
	#	$tmp .= ("# Nickase Recognition Site " . eval($i+1) . ":\t$enzyme_sequences[$i]\n");
	#}
	
	foreach my $key (sort {$a <=> $b} keys %color_to_enzyme_sequence){
		$tmp .= ("# Nickase Recognition Site $key:\t" . $color_to_enzyme_sequence{$key}[0]);
		
		for(my $i = 1; $i < scalar(@{$color_to_enzyme_sequence{$key}}); $i++){
			$tmp .= ",$color_to_enzyme_sequence{$key}[$i]";
		}
		
		$tmp .= "\n";
	}
	
	chomp($tmp);
	
	my $str = << "EOF";
# CMAP File Version:	0.1
# Label Channels:	$color_uniq
$tmp
# Number of Consensus Nanomaps:	N/A
#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence$strand_header
#f int	float	int	int	int	float	float	int	int$strand_header_type
EOF
	print $OUT($str);
	close($OUT);
}

sub Generate_cmap{
	my ($filename, $ID, $loci_ref, $length) = @_;
	my $i;
	my $siteID = 1;
	my $total_loci = 0;
	my $OUT;
	my $length_float = sprintf("%.1f", $length);
	my @sorted_loci = sort {$a <=> $b} @$loci_ref;
	my @positions;  # Two-dimensional array containing label positions of each channel
	
	open($OUT, ">>$filename") || die ("ERROR: Can't open $filename: $!\n");
	
	for($i = 0; $i < @sorted_loci; $i++){
		# "1" in most of the time, ">1" when there are multiple labels at the same position (the array is unique!)
		$total_loci += $color{$sorted_loci[$i]}[0];
	}
	
	for($i = 0; $i < @sorted_loci; $i++){
		my $loci_float = sprintf("%.1f", $sorted_loci[$i]);
		my @sorted_colors = sort {$a <=> $b} @{$color{$sorted_loci[$i]}}[1 .. $color{$sorted_loci[$i]}[0]];
		
		for(my $j = 0; $j < $color{$sorted_loci[$i]}[0]; $j++){
			if($strand){
				print $OUT("$ID\t$length_float\t$total_loci\t$siteID\t", $sorted_colors[$j], "\t$loci_float\t1.0\t1\t1\t", $strand{$sorted_loci[$i]}, "\n");
			}
			else{
				print $OUT("$ID\t$length_float\t$total_loci\t$siteID\t", $sorted_colors[$j], "\t$loci_float\t1.0\t1\t1\n");
			}
			$num_nicks[$sorted_colors[$j]]++;
			$num_nicks_global[$sorted_colors[$j]]++;
			
			$siteID++;
			
			# Storing label positions for calculating the label distances
			push(@{$positions[$sorted_colors[$j]]}, $sorted_loci[$i]);
		}
	}
	
	# Calculating per channel label distances
	foreach my $channel (@color_uniq){
		if(exists($positions[$channel][1])){
			push(@{$label_distances_by_channel[$channel]}, Get_distance(\@{$positions[$channel]}));
			$total_label_distance[$channel] += ($positions[$channel][-1] - $positions[$channel][0]);
		}
	}
	
	# Calculating global label distances
	if(@sorted_loci >= 2){
		push(@label_distance_filtered, Get_distance(\@sorted_loci));
		$total_label_distance_filtered += ($sorted_loci[-1] - $sorted_loci[0]);
	}
	
	if($strand){
		print $OUT("$ID\t$length_float\t$total_loci\t$siteID\t0\t$length_float\t0.0\t1\t0\t.\n");
	}
	else{
		print $OUT("$ID\t$length_float\t$total_loci\t$siteID\t0\t$length_float\t0.0\t1\t0\n");
	}
	$num_cmaps++;
	
	close($OUT);
}

sub Get_distance{
	my ($loci_ref) = @_;
	my @dist;
	for(my $i = 0; $i < @$loci_ref - 1; $i++){
		push(@dist, $loci_ref->[$i+1] - $loci_ref->[$i]);
	}
	return @dist;
}

sub Get_and_set_enzyme_sequence{
	my ($hash_ref, $str_ref) = @_;
	
	foreach my $item (keys %$hash_ref){
		if(uc(substr($item, 0, 3)) eq uc(substr($$str_ref, 0, 3))){
			$$str_ref = $item;
			return uc($hash_ref -> {$item});
		}
	}
	
	print("ERROR: Invalid parameter(s)!\n");
	Usage();
}

sub is_enzyme_name{
	my ($hash_ref, $str) = @_;
	
	my @array = map { uc(substr($_, 0, 3)) } keys %$hash_ref;
	
	if(uc(substr($str, 0, 3)) ~~ @array){
		return 1;
	}
	else{
		return 0;
	}
}

sub is_nt_sequence{
	my ($str) = @_;
	
	for(my $i = 0; $i < length($str); $i++){
		if("ACGTacgt" !~ substr($str, $i, 1)){
			return 0;
		}
	}
	return 1;
}

sub is_enzyme_name_and_sequence{
	my ($str) = @_;
	
	if($str =~ ":"){
		my @name_seq = split(":", $str);
		if(@name_seq == 2 && is_nt_sequence($name_seq[1])){
			return 1;
		}
	}
	return 0;
}

sub Uniq{
	my %seen;
	grep {!$seen{$_}++} @_;
}

sub Uniq_BAK{
	my %seen;
	my @arr_uniq;
	
	foreach (@_){
		if(!exists($seen{$_})){
			$seen{$_} = undef;
			push(@arr_uniq, $_);
		}
	}
	
	return @arr_uniq;
}

sub Get_N50{
	my ($data, $total_len) = @_;
	my @sorted = sort {$b <=> $a} @$data;
	
	my $total = 0;
	my $i;
	
	for($i = 0; $i < @sorted; $i++){
		$total += $sorted[$i];
		if($total >= $total_len / 2){
			return $sorted[$i];
		}
	}
}

sub Get_N90{
	my ($data, $total_len) = @_;
	my @sorted = sort {$b <=> $a} @$data;
	
	my $total = 0;
	my $i;
	
	for($i = 0; $i < @sorted; $i++){
		$total += $sorted[$i];
		if($total >= $total_len * 0.9){
			return $sorted[$i];
		}
	}
}

sub StdDev{
	my ($array_ref) = @_;
	my ($sum, $sum2) = (0, 0);
	
	if(!@$array_ref){
		return 0;
	}
	
	foreach my $item (@$array_ref){
		$sum += $item;
	}
	my $mean = $sum / @$array_ref;
	
	foreach my $item (@$array_ref){
		$sum2 += ($item-$mean)**2;
	}
	
	return sqrt($sum2/@$array_ref);
}

# Created by Zhanyang Zhu
sub getNGaps{
	my ($seq, $minNGapSize) = @_;
	my @gaps;
	# Assert: $minNGapSize should not be less than 1
	while($seq =~ /((N|n){$minNGapSize,})/g){
		push @gaps, [$-[0] + 1, $+[0]];
	}
	return (\@gaps);
}

sub printNGaps{
	my ($ngap_ref, $ngapFileName, $minNGapSize) = @_;
	open NG, ">$ngapFileName" || die ("ERROR: Can't open $ngapFileName for writing: $!\n");
	print NG "#min N gap size: $minNGapSize\n";
	my @chromosomes = sort {$a <=> $b} (keys %$ngap_ref);
	foreach my $chromosome (@chromosomes){
		my $c_gaps_ref= $ngap_ref->{$chromosome};
		my $num_gaps = scalar @$c_gaps_ref;
		for(my $i = 0; $i < $num_gaps; $i++){
            my $start = $ngap_ref->{$chromosome}->[$i]->[0] - 1;
            my $end = $ngap_ref->{$chromosome}->[$i]->[1] - 1;
			print NG "$chromosome\t$start\t$end\tgap\n";
		}
	}
	close(NG);
}

sub GetNumLabelByDist{
  	 my ($distLabel) = @_;
   	 my $saphyr = 1000;    
  	 my $irys = 1500;

  	 my $sumD = 0; ##sum of low dist
  	 my $llc = 0; ##low dist label counts
  	 my @numLabel = [1, 1];
  	 foreach my $d (@{$distLabel}) {
	   if( $d >= $irys) {
		 $numLabel[1]++;
		 $numLabel[2]++;
 		 $sumD = 0;
	 	 $llc = 0;
 	   }
 	   elsif($d >= $saphyr && $d < $irys) {
 		$numLabel[1]++;
 		$sumD = 0;
		 $llc = 0;
	   }
 	   else{
 		$sumD += $d;
		$llc++;
 	   }
	
	  if( $llc > 0 && ($sumD == 0 || $sumD > $saphyr) ){ ## add 1 for two consecutive close labels or more consecutive labels which combining distance bigger than resolution.
 		$numLabel[1]++;
 		$numLabel[2]++;
	 	$sumD = 0;
		$llc = 0;
 	  }
 	}

	return @numLabel;
}


__END__

/home/users/xzhou/PERL/fa2cmap_multi_color.pl -i XXX.fa -e bspq 2 CACTTAAA 1 -o ./

svn propset svn:keywords "Id" fa2cmap_multi_color.pl

=COMMENT
use File::Spec;
use File::Spec::Functions;
File::Spec -> rel2abs($f)

use File::Basename;
my $dir = dirname($f_in_abs);
my $filename = basename($f_in_abs);
print "dir = $dir\n";
print "filename = $filename\n";

use File::Spec::Functions;
$f_log = catfile($dir, $filename);
=cut

