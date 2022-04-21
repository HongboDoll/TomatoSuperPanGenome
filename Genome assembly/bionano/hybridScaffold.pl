# $Id: hybridScaffold.pl 8926 2019-04-23 19:48:33Z jwang $
#!/usr/bin/perl -w
#################################################################################################################################################################################
# File: hybridScaffold.pl                                                                  					
# Date: 07/24/2014                                                                         				
# Purpose: Merge NGS scaffolds with BioNano CMAP                                           				
#                                                                                          			
#                                                                                          		
# Usage:                                                                                   	
#  hybridScaffold.pl <-h> <-n ngs_file> <-b bng_cmap_file> <-c hybrid_config_xml> <-o output_folder> <-B conflict_filter_level> <-N conflict_filter_level> <-f> 
#     <-m molecules_bnx> <-p de_novo_pipeline> <-q de_novo_xml> <-v><-x> <-y> <-e noise_param><-z tar_zip_file><-S><-w enzyme_name>
#     -h    : This help message                                                            											
#     -n    : Input NGS FASTA [required]                                      												
#     -b    : Input BioNano CMAP [required]                                                											
#     -c    : Merge configuration file [required]                                          											
#     -o    : Output folder [required]                                                     											
#     -r    : RefAligner program [required]
#     -B    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contigs [required if not using -M option]
#     -N    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contigs [required if not using -M option]
#     -f    : Force output and overwrite any existing files                                											
#     -x    : Flag to generate molecules to hybrid scaffold alignment and molecules to genome map alignment [optional]
#     -y    : Flag to generate chimeric quality score for the Input BioNano CMAP [optional]
#     -m    : Input BioNano molecules BNX [optional; only required for either the -x or -y option]
#     -p    : Input de novo assembly pipeline directory [optional; only required for -x option]
#     -q    : Input de novo assembly pipeline optArguments XML script [optional; only required for -x option]
#     -e    : Input de novo assembly noise parameter .errbin or .err file [optional; recommended for -y option but not required]
#     -v    : Print pipeline version information																
#     -M    : Input a conflict resolution file indicating which NGS and BioNano conflicting contigs to be cut [optional]
#     -z    : Name of a zipped file to archive the essential output files [optional]
#     -a    : Do not zipped file to archive the essential output files (archive is by default on now)
#     -S    : Only run hybridScaffold up to before Merge steps [optional]
#     -w    : Name of the status text file needed by IrysView [optional]
#     -t    : Perform pre-pairmerge sequence to pre-pairmerge genome map alignment [optional]
#     -u    : Sequence of enzyme recognition site (overrides what has been specified in config XML file, for IrysView only) [optional]
#
#################################################################################################################################################################################

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

# This adds "${CURRENT_SCRIPT_PATH}/scripts/perl5/" direcory at run time to the @INC array
# This script sould sit two levels above the additional Perl modules directory
BEGIN {
	my $script_path = abs_path(dirname($0));
	my $module_path = abs_path($script_path . "/scripts/perl5");
	unshift @INC, $module_path;
	my $lib3;
	if ($] >= 5.010000 && $] <= 5.011000) {
		$module_path = $module_path."/5.10.1"; 
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.014000 && $] <= 5.015000) {
		$module_path = $module_path."/5.14.4"; 
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.016000 && $] <= 5.017000) {
        	$module_path = $module_path."/5.16.3";
        	$lib3 = $module_path."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.018000 && $] <= 5.019000) {
		$module_path = $module_path."/5.18.2";
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	else {
		print "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n";
		die "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n"; 
		exit 1; }
	unshift @INC, $module_path;
	unshift @INC, $lib3;
	#print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}

use IO::Handle;
use PerlIO::Util;
use Getopt::Std;
use Config::Simple;
use File::Path qw(mkpath rmtree);
use File::Slurp;
use File::Copy qw(copy move);
use File::Copy::Recursive qw(fcopy dircopy);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
use XML::Simple;
use DateTime;
use DateTime::Format::Human::Duration;
use File::Basename;
use threads;
use IO::Select;
use IPC::Open3;
use File::Spec;
use File::Find;
sub Init;
sub Usage;
sub CHECK_Config;
sub Version;
sub find_refaligner;
sub log10;
use List::MoreUtils qw(uniq);
use BNG::Utility;

my %opts;
my @cmd = ();	
my $cmdRef = \@cmd;
my ($outResults, $errResults) = ("", "");
my $plDir=abs_path(dirname($0));


#read comand line args
Init();

if(defined($opts{S})){
    print "Running pre-processing only\n";
}

### cur_results.txt ###
# this file indicates where the final results should be
my $hybridScaffoldFinalResultsDir = "hybrid_scaffolds";	
my $modifyNum = getCurMaxModifyNum($opts{o});	# to be used to distinctly identify the number of times a modified conflict status file was fed (ignored if no such file was inputted)
$hybridScaffoldFinalResultsDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");
$hybridScaffoldFinalResultsDir = abs_path("$opts{o}/$hybridScaffoldFinalResultsDir");
eval {mkpath($hybridScaffoldFinalResultsDir)};
if ($@)	{
	print("ERROR: Could not create final output directory $hybridScaffoldFinalResultsDir: $@");
	print("ERROR: Output folder invalid or does not exist!\n");
	Usage();
}
open(OUT_DIR_POINTER, ">$hybridScaffoldFinalResultsDir/cur_results.txt") or dieLog "ERROR: Cannot write a file to indicate which hybrid scaffold final results dir is: $!\n";
print OUT_DIR_POINTER "$hybridScaffoldFinalResultsDir\n";
if (defined($opts{m}))	{
	print OUT_DIR_POINTER "$opts{o}/alignmol_bionano";	print OUT_DIR_POINTER "_M$modifyNum" if (defined($opts{M}));	print OUT_DIR_POINTER "\n";
	print OUT_DIR_POINTER "$opts{o}/alignmol_hybrid";	print OUT_DIR_POINTER "_M$modifyNum" if (defined($opts{M}));	print OUT_DIR_POINTER "\n";
}
my $statusFile = (defined($opts{w})) ? ($opts{w}) : ("status.txt");	
print OUT_DIR_POINTER "$opts{o}/$statusFile\n";	# print the status file
close OUT_DIR_POINTER;

my @s = split("/", $opts{o});

### log file ###
my @NGSpath = split(/\//, $opts{n});
my $NGSfile = $NGSpath[$#NGSpath];
$NGSfile =~ s/\./_/g;	# replaces all . to underscores

my @BNGpath = split(/\//, $opts{b});
my $BNGfile = $BNGpath[$#BNGpath];	# get the actual file name
$BNGfile =~ s/\./_/g;	# replaces all . to underscores

my $temp = $BNGfile;	$temp =~ s/_cmap$//;	# rid the file extension
my $modBNGfile = $temp."_bppAdjust_cmap";
my $log_file2 = "$hybridScaffoldFinalResultsDir/$modBNGfile"."_$NGSfile"."_HYBRID_SCAFFOLD_log.txt";
for (*STDOUT, *STDERR)	{
	# captures the stdout and stderr 
	$_->autoflush;	$_->push_layer(tee=>$log_file2)
} # for	

### status txt (a sub set of the log file) ###
$statusFile = (defined($opts{w})) ? ($opts{w}) : ("status.txt");	
open(STATUS_FILE, ">$opts{o}/$statusFile") or dieLog "ERROR: Cannot write to status file $opts{o}/$statusFile: $!\n";	# if there is an existing $opts{o}/$statusFile, then that would be over-written

	
#print header
printStatus(\*STATUS_FILE, "initial check", "Started", 0);

print "**************************\n";
print "*****BioNano Genomics*****\n";
print "******BNG-NGS Merge*******\n";
print "**************************\n\n";

Version();
print "\n";

my $hostname = `hostname`;
print "Running on host: $hostname\n";
print "Process ID (PID): $$\n";

my $dtStart = DateTime->now;
print "Start time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n\n";

print qx/ps -o args $$/;

#make sure all input files exist 
print "\n";
if (-e $opts{n} && -e $opts{b} && -e $opts{c}) {
	
	print "NGS file: $opts{n} \n" ;
	print "BNG file: $opts{b} \n";
	print "Configuration file: $opts{c} \n"; 
	print "Output folder: $opts{o}\n"; 
	print "Molecules file: $opts{m}\n" if (defined($opts{m}));
	print "User-defined conflict status file: $opts{M}\n" if (defined($opts{M}));
} else {
	print("One or more input files do not exist. \n");
	Usage(); 
}

#make sure input BNG CMAP file contains at least one contig
my (undef, $numContig, undef) = readCMap($opts{b});
if (!($numContig > 0)) {
	my $errMsg = "ERROR: Input BNG file $opts{b} does not contain any contigs!";
	printStatus(\*STATUS_FILE, "initial check", $errMsg, 1);
	dieLog( "$errMsg\n\n");
}

### RefAligner section
my $refaligner = $opts{r};
if ($refaligner =~ /~/)	{
	$refaligner = glob($refaligner);
} else	{
	$refaligner = abs_path($refaligner) if (-e $refaligner);
} # if refaligner

# make sure RefAligner exists
if (! -e $refaligner || ! -x $refaligner || $refaligner !~ /RefAligner$/i)	{
	print "WARNING: RefAligner binary does not exist at $refaligner. Trying to find it...\n";
	my $script_path = abs_path(dirname($0));
	my @s = split('/',$script_path);
	my $val = pop(@s); $val = pop(@s);
	my $home_path = join('/',@s);
	$refaligner = $home_path."/tools/RefAligner";
	
	if (! -e $refaligner || ! -x $refaligner || $refaligner !~ /RefAligner$/i)	{
		$refaligner = glob("~/tools/RefAligner");
		if (! -e $refaligner || ! -x $refaligner || $refaligner !~ /RefAligner$/i)	{
			my $errMsg = "ERROR: RefAligner binary cannot be found!";
			printStatus(\*STATUS_FILE, "initial check", $errMsg, 1);
			dieLog ("$errMsg\n");
		} # if -e refaligner (final try)
	} # if -e refaligner (second try)
} # if -e refaligner (first try)
print "RefAligner binary: $refaligner\n";
### end of RefAligner section

### config file section
#load config file
my $XML = new XML::Simple(KeyAttr=>[]);
my $configRef = $XML->XMLin($opts{c});
my %stageStack = ();
my @commandStack = ();
my $commandStackRef = \@commandStack;
my $baseStage = "";
copy "$opts{c}", "$hybridScaffoldFinalResultsDir" ;
### end of config file section

$baseStage = "global";	
$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
my $theParamValuesRef;
$theParamValuesRef = extractConfig($commandStackRef, "maxmem");	my $maxmem = $theParamValuesRef->[0];
$theParamValuesRef = extractConfig($commandStackRef, "maxthreads");	my $maxthreads = $theParamValuesRef->[0];
# reset
$baseStage = "";	
@$commandStackRef = ();
printStatus(\*STATUS_FILE, "initial check", "Done", 0);


# define some variables that will be needed regardless of whether a conflict status file is fed in
$baseStage = "fasta2cmap";
$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
my @enzymes = ();
if (defined ($opts{u}))	{
	# the optional enzyme information is passed in via command line
	push(@enzymes, uc($opts{u}) );
} else	{
	# just read enzyme information in the XML file
	$theParamValuesRef = extractConfig($commandStackRef, "enzyme");	
	for (my $i = 0; $i < scalar(@$theParamValuesRef); $i += 1)	{
		$enzymes[$i] = uc($theParamValuesRef->[$i]);	# capitalize the entire enzyme name
	} # for i
} # if defined
$theParamValuesRef = extractConfig($commandStackRef, "channelNum");	my $channelNumber = $theParamValuesRef->[0];
$theParamValuesRef = extractConfig($commandStackRef, "minLabels");	my $minLabels = $theParamValuesRef->[0];
$theParamValuesRef = extractConfig($commandStackRef, "minLength");	my $minLength = $theParamValuesRef->[0];
$baseStage = "";
@$commandStackRef = ();

my @file = split(/\./, $opts{n});	my $filename = $opts{n};	$filename =~ s/(\S+)\.\w+$/$1/;	foreach my $enzyme (@enzymes) {$filename .= "_$enzyme"};
my $filename_key = $filename."_".$minLength."kb_".$minLabels."labels_key.txt";
my $ngs_cmap = $filename."_".$minLength."kb_".$minLabels."labels.cmap";
my $bppAdjust = "";					# adjust the bpp of BioNano genome maps after align to sequence in align0
my $readparameters = "";				# to store the noise parameters in align1
my ($assignAlign_r, $assignAlign_q) = ("", "");		# to store those contigs that passed the contig removal process
my ($cutBNFile, $cutSeqFile) = ("", "");
my $outputDir = "";
my ($dtStartStage, $dtEndStage) = (DateTime->now, DateTime->now);
my $span = DateTime::Format::Human::Duration->new();

if (! defined($opts{M}))	{
	# if the user did not incorporate a conflict status file manually, then the pipeline starts normally from the beginning

	my $chimResultDir = "$opts{o}/chim_qual";
	# first check if chimeric quality scores are needed
	if ($opts{B} == 2 || $opts{N} == 2)	{
		# if cut was indicated to be done, then the presence of chimeric quality score is checked, else reset the option to remove conflict
		my %qScores = ();	my $qScoresRef = \%qScores;	my $noQScoreFlag = 0;
		($qScoresRef, $noQScoreFlag) = getQScores($opts{b}, $qScoresRef);

		if ($noQScoreFlag == 1)	{
			print "WARNING: cut at conflict option was selected, but no chimeric quality score is available in the genome map...\n";
			# the BioNano genome maps do not have chimeric quality score
			if ($opts{y})	{
				printStatus(\*STATUS_FILE, "chimeric quality score generation", "Started", 0);
				# user indicates the generation of chimeric quality score
				print "chimeric quality score will be generated\n";	
				$dtStartStage = DateTime->now;
				print "\nBeginning calculation of chimeric score...\n";
				my $chimResultDir = "$opts{o}/chim_qual";
				eval	{	mkpath $chimResultDir	};
				if ($@)	{
					print "Couldn't create $chimResultDir: $@\n";
				} # if

				my $deNovoOpArgsFile = $opts{q};
				my $deNovoNoiseFile = $opts{e};
				my $moleculesDir = $opts{m};
				my ($inMapDir, $inMapFile) = separateDirFile($opts{b});

				unless (chdir $plDir)	{
					my $errMsg = "ERROR: cannot change directory to $plDir: $!";
					printStatus(\*STATUS_FILE, "chimeric quality score generation", $errMsg, 1);
					dieLog ("$errMsg\n");
				} # chdir
				if (defined($opts{e}))	{
					# if a noise parameter file has been provided
					@$cmdRef = ($^X, "scripts/calc_chim_score.pl", "-refAligner", $refaligner, "-inMapFile", $opts{b}, "-outMapFile", "$chimResultDir/$inMapFile", "-xmlFile", $opts{c}, "-noiseParamFile", $deNovoNoiseFile, "-bnxFile", $moleculesDir, "-maxmem", $maxmem);
				} else	{
					# not provided
					@$cmdRef = ($^X, "scripts/calc_chim_score.pl", "-refAligner", $refaligner, "-inMapFile", $opts{b}, "-outMapFile", "$chimResultDir/$inMapFile", "-xmlFile", $opts{c}, "-bnxFile", $moleculesDir, "-maxmem", $maxmem);
				} # if defined opts{e}

				print "Running command: ".(join(" ", @$cmdRef))."\n";
				($outResults, $errResults) = runCommand($cmdRef);
				@$cmdRef = ();
				if (open(OUT, ">$chimResultDir/chim_score.log"))	{
					print OUT $outResults;	close OUT;
				} else	{
					my $errMsg = "ERROR: Cannot write to $chimResultDir/chim_score.log: $!";
					printStatus(\*STATUS_FILE, "chimeric quality score generation", $errMsg, 1);	dieLog ("$errMsg\n");
				} # if open
				if (open(ERR, ">$chimResultDir/chim_score.errlog"))	{
					print ERR $errResults;	close ERR;
				} else	{
					my $errMsg = "ERROR: Cannot write to $chimResultDir/chim_score.errlog: $!";
					printStatus(\*STATUS_FILE, "chimeric quality score generation", $errMsg, 1);	dieLog("$errMsg\n");
				} # if open
				$dtEndStage = DateTime->now;
				$span = DateTime::Format::Human::Duration->new();
				$span = $span->format_duration_between($dtEndStage, $dtStartStage);
				# error check
				if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
					my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
					printStatus(\*STATUS_FILE, "chimeric quality score generation", $errMsg, 1);
				} else	{
					printStatus(\*STATUS_FILE, "chimeric quality score generation", "Done", 0);
				} # if outResults
				errCheck($outResults, $errResults, "ERROR:", "calculation of chimeric score completed in $span.", "ERROR: calculation of chimeric score cannot be completed.");
			} else	{
				# user has not indicated the generation of chimeric quality score
				print "chimeric quality score option was NOT indicated, so ...\n";

				if ($opts{B} == 2)	{
					$opts{B} = 3;
					print "WARNING: when conflicts occur, BioNano maps are thrown out instead of cut\n";
				} # if opts{B}
				if ($opts{N} == 2)	{
					$opts{N} = 3;
					print "WARNING: when conflicts occur, sequences are thrown out instead of cut\n";
				} # if opts{N}
				print "WARNING: B option is now $opts{B} and N option is now $opts{N}\n";
			} # if opts

		} #if noQScoreFlag  

		

	} # if opts{B} or opts{N}

	#run fasta to cmap conversion 
	printStatus(\*STATUS_FILE, "FASTA to CMAP conversion", "Started", 0);

	$dtStartStage = DateTime->now;
	eval{ mkpath("$opts{o}/fa2cmap") };	print "Couldn't create $opts{o}/fa2cmap: $@" if ($@);
	print "\nBeginning FASTA to CMAP conversion...\n";
	print "Using Enzyme: ".(join(" ", @enzymes))."\n";
	print "Minimum Length: $minLength Kb\nMinimum Labels: $minLabels\n";

	chdir $plDir or dieLog( "ERROR: Cannot change directory to $plDir: $!\n");
	@$cmdRef = ($^X, "scripts/fa2cmap_multi_color.pl", "-i", $opts{n}, "-m", $minLabels, "-M", $minLength, "-o", "$opts{o}/fa2cmap", "-e");
	foreach my $enzyme (@enzymes)	{
		push(@$cmdRef, $enzyme, $channelNumber);
	} # foreach enzyme
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand ($cmdRef);
	@$cmdRef = ();
	open (OUT, ">$opts{o}/fa2cmap/fa2cmap.log"); print OUT $outResults."\n"; close OUT;
	open (ERR, ">$opts{o}/fa2cmap/fa2cmap.errlog"); print ERR $errResults."\n"; close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	
	# error check
	
	#check uniqueness of input NGS sequence accession
	my $fasta_dir = $opts{o}."/fa2cmap";
	opendir(DIR, $fasta_dir);
	my @files = grep{/\_key\.txt$/}readdir DIR;
	my $key_file = $fasta_dir."/".$files[0];
	print "key file $key_file\n";
	open(KEY, $key_file) || dieLog ("cannot open fasta2cmap key file");
	my %key_table = ();
	while(<KEY>){
		my @tokens = split(/\t/, $_);
		if($tokens[1]){
			my $accession = (split(/\s+/, $tokens[1]))[0];
			#print "accession is $accession\n";
			if($accession && $key_table{$accession}){
				my $errMessage = "ERROR: Input NGS sequence contain duplicate accession name: $accession";
				printStatus(\*STATUS_FILE, "FASTA to CMAP conversion", $errMessage, 1);
				dieLog("$errMessage\n")
			}
			$key_table{$accession} = 1;
		}
	}
	#other error check
	if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
		my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
		printStatus(\*STATUS_FILE, "FASTA to CMAP conversion", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "FASTA to CMAP conversion", "Done", 0);
	} # if outResults
	errCheck($outResults, $errResults, "ERROR:", "FASTA to CMAP conversion complete in $span.", "ERROR: FASTA to CMAP conversion cannot be completed.");

	
	# convert FASTA header to cmap id
	printStatus(\*STATUS_FILE, "FASTA header to CMAP id conversion in FASTA file", "Started", 0);
	$dtStartStage = DateTime->now;
	print "\nBeginning FASTA header conversion...\n";
	unless (chdir $plDir)	{
		my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
		printStatus(\*STATUS_FILE, "FASTA header to CMAP id conversion in FASTA file", $errMsg, 1);
		dieLog ("$errMsg\n");
	} # unless chdir

	$filename_key = basename($opts{n});
	$filename_key =~ s/(\S+)\.\w+$/$1/;	foreach my $enzyme (@enzymes) {$filename_key .= "_$enzyme"};
	$filename_key = $filename_key."_".$minLength."kb_".$minLabels."labels_key.txt";
	$filename_key = $opts{o}."/fa2cmap/".$filename_key;
	$filename = basename($opts{n});	$filename =~ s/(\S+)\.\w+$/$1/; foreach my $enzyme (@enzymes) {$filename .= "_$enzyme"};
	$filename = $filename."_".$minLength."kb_".$minLabels."labels.cmap";
	$filename = $opts{o}."/fa2cmap/".$filename;

	@$cmdRef = ($^X, "scripts/fa_key_convert.pl", $opts{n}, $filename_key);
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand ($cmdRef);
	@$cmdRef = ();
	open (OUT, ">$opts{o}/fa2cmap/faHeader_to_cmapId.log"); print OUT $outResults."\n"; close OUT;
	open (ERR, ">$opts{o}/fa2cmap/faHeader_to_cmapId.errlog"); print ERR $errResults."\n"; close ERR;	
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check	
	if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
		my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
		printStatus(\*STATUS_FILE, "FASTA header to CMAP id conversion in FASTA file", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "FASTA header to CMAP id conversion in FASTA file", "Done", 0);
	} # if outResults
	errCheck($outResults, $errResults, "ERROR:", "FASTA header conversion complete in $span.", "ERROR: FASTA header conversion cannot be completed.");
	$filename =~ s/(\S+)\.\w+$/$1/;
	$temp = $filename."_CmapIdHeaders.fa";	
	print "\nNew FASTA with CMAP Ids as headers: $temp\n";

	# copy the fasta file to the designated output directory fa2cmap sub-directory
	printStatus(\*STATUS_FILE, "align0", "Started", 0);
	copy "$opts{n}", "$opts{o}/fa2cmap" ;
	$ngs_cmap = $filename.".cmap";

	$ngs_cmap = abs_path($ngs_cmap);
	print "NGS map path: $ngs_cmap\n";

	#make sure input NGS CMAP file contains at least one contit
	(undef, $numContig, undef) = readCMap($ngs_cmap);
	if (!($numContig > 0)) {
		my $errMsg = "ERROR: Input NGS file $ngs_cmap does not contain any contigs!";
		printStatus(\*STATUS_FILE, "align0", $errMsg, 1);
		dieLog( "$errMsg\n\n");
	}

	#perform initial alignment 
	$dtStartStage = DateTime->now;
	print "\nBeginning initial NGS CMAP to BioNano CMAP alignment...\n";
	$outputDir = $opts{o}."/align0";
	eval { mkpath($outputDir) };
	if ($@) {
		print "Couldn't create $outputDir: $@"; }

	# for BioNano maps, check if chimeric quality score cmap file has been generated
	my $align0BNGInput = "";
	my ($inBNGDir, $inBNGFile) = separateDirFile($opts{b});
	if (-e $chimResultDir && -e "$chimResultDir/$inBNGFile")	{
		$align0BNGInput = "$chimResultDir/$inBNGFile";
	} else	{
		$align0BNGInput = $opts{b};
	} # if chimeric quality score cmap file has been generated
	copy "$align0BNGInput" , "$outputDir";

	unless (chdir $outputDir)	{
		my $errMsg = "ERROR: Cannot change directory to $outputDir: $!";
		printStatus(\*STATUS_FILE, "align0", $errMsg, 1);
	} # unless chdir
	@$cmdRef = ($refaligner, "-ref", $ngs_cmap, "-i", $align0BNGInput, "-o", "align0", "-stdout", "-stderr");
	$baseStage = "global";
	$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	$baseStage = "align0";
	$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	if (scalar(@$commandStackRef) == 0)	{
		# no align0 block available
		$baseStage = "align1";
		$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	} # if scalar
	$cmdRef = copyArray($cmdRef, $commandStackRef);	# the rest of the parameters
	$baseStage = "";
	@$commandStackRef = ();
	
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	if (errCheckRefAligner("align0.stdout", "END of output", "Initial alignment complete in $span.", "ERROR: Initial alignment cannot be completed.", 0) != 1)	{
		my $errMsg = "ERROR: Initial alignment cannot be completed.";
		printStatus(\*STATUS_FILE, "align0", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "align0", "Done", 0);
	} # RefAligner check
	errCheckRefAligner("align0.stdout", "END of output", "Initial alignment complete in $span.", "ERROR: Initial alignment cannot be completed.", 1);

	my $align0_errbin = $outputDir."/align0.errbin";

	#rescale BNG input cmap
	printStatus(\*STATUS_FILE, "Rescale BioNano maps", "Started", 0);
	$dtStartStage = DateTime->now;
	print "\nRescaling BioNano CMAP...\n";
	unless (chdir $outputDir)	{
		my $errMsg = "ERROR: Cannot change directory to $outputDir: $!";
		printStatus(\*STATUS_FILE, "Rescale BioNano maps", $errMsg, 1);
		dieLog("$errMsg\n");
	} # unless chdir
	my @suffixes = {".cmap",".CMAP"};
	$bppAdjust = $modBNGfile;	$bppAdjust =~ s/_cmap$//;	# remove the suffix
	$bppAdjust = $opts{o}."/align0/".$bppAdjust;
	@$cmdRef = ($refaligner, "-merge", "-i", $align0BNGInput, "-o", $bppAdjust, "-readparameters", $align0_errbin, "-stdout", "-stderr");
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	if (errCheckRefAligner($bppAdjust.".stdout", "END of output", "Rescaling complete in $span.", "ERROR: Rescaling cannot be completed.", 0) != 1)	{
		my $errMsg = "ERROR: Rescaling cannot be completed.";
		printStatus(\*STATUS_FILE, "Rescale BioNano maps", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "Rescale BioNano maps", "Done", 0);
	} # if RefAligner check
	errCheckRefAligner($bppAdjust.".stdout", "END of output", "Rescaling complete in $span.", "ERROR: Rescaling cannot be completed.", 1);
	$bppAdjust = $bppAdjust.".cmap";


	#perform initial alignment using rescaled BNG
	printStatus(\*STATUS_FILE, "align1", "Started", 0);
	$dtStartStage = DateTime->now;
	print "\nBeginning initial NGS CMAP to rescaled BioNano CMAP alignment...\n";
	$outputDir = $opts{o}."/align1";
	eval { mkpath($outputDir) };
	if ($@) {
	  print "Couldn't create $outputDir: $@"; }

	unless (chdir $outputDir)	{
		my $errMsg = "ERROR: Cannot change directory to $outputDir: $!";
		printStatus(\*STATUS_FILE, "align1", $errMsg, 1);
		dieLog("$errMsg\n");
	} # unless chdir
	@$cmdRef = ($refaligner, "-ref", $ngs_cmap, "-i", $bppAdjust, "-o", "align1", "-stdout", "-stderr");
	$baseStage = "global";
	$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	$baseStage = "align1";
	$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	$cmdRef = copyArray($cmdRef, $commandStackRef);
	$baseStage = "";
	@$commandStackRef = ();
	
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	if (errCheckRefAligner("align1.stdout", "END of output", "Initial rescaled alignment complete in $span.", "ERROR: Initial rescaled alignment cannot be completed.", 0) != 1)	{
		my $errMsg = "ERROR: Initial rescaled alignment cannot be completed.";
		printStatus(\*STATUS_FILE, "align1", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "align1", "Done", 0);
	} # if RefAligner check
	errCheckRefAligner("align1.stdout", "END of output", "Initial rescaled alignment complete in $span.", "ERROR: Initial rescaled alignment cannot be completed.", 1);

	my $align1_xmap = $outputDir."/align1.xmap";
	my $align1_r_cmap = $outputDir."/align1_r.cmap";
	my $align1_q_cmap = $outputDir."/align1_q.cmap";
	$readparameters = abs_path("$opts{o}/align1/align1.errbin");


	#check to make sure that intial alignment XMAP contains alignments
	printStatus(\*STATUS_FILE, "conflict detection", "Started", 0);
	my @int_xmap = read_file($align1_xmap);
	my $headLines=0;
	my $alignLines=0;
	foreach (@int_xmap) {
		if ($_ =~ /^#/ ) {
			$headLines++; }
		else {
			$alignLines++; } }
	if ($alignLines < 1) {
		my $errMsg = "ERROR: No intial alignments found between $opts{n} and $bppAdjust";
		printStatus(\*STATUS_FILE, "conflict detection", $errMsg, 1);
		dieLog("\n$errMsg\n");
	}	else {
		print "\n$alignLines alignments found between $opts{n} and $bppAdjust\n"; 
	} # if no alignment possible between sequences and maps

	#run AssignAlignType script
	$dtStartStage = DateTime->now;
	print "\nBeginning AssignAlignType...\n";
	$outputDir = $opts{o}."/assignAlignType";
	eval { mkpath($outputDir) };
	if ($@) {
	  print "Couldn't create $outputDir: $@"; }
	my $assignAlignType_xmap = $outputDir."/assignAlignType.xmap";
	my $assignAlignType_r_cmap = $outputDir."/assignAlignType_r.cmap";
	my $assignAlignType_q_cmap = $outputDir."/assignAlignType_q.cmap";
	my $conflict_r_cmap = $outputDir."/conflict_ngs_r.cmap";
	my $conflict_q_cmap = $outputDir."/conflict_bng_q.cmap";
	my $breakPoint_txt = $outputDir."/conflicts.txt";
	$baseStage = "assignAlignType";
	$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	$theParamValuesRef = extractConfig($commandStackRef, "T_cutoff");	my $T_cutoff = $theParamValuesRef->[0];	$T_cutoff = -log10($T_cutoff);
	$theParamValuesRef = extractConfig($commandStackRef, "max_overhang");	my $max_overhang = $theParamValuesRef->[0];
	$baseStage = "";
	@$commandStackRef = ();
	my $breakPointFileShiftAmount = 30;	# the amount to shift at the conflict breakpoints so that the first/last label of the alignment is preserved

	unless (chdir $plDir)	{
		my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
		printStatus(\*STATUS_FILE, "conflict detection", $errMsg, 1);
		dieLog("$errMsg\n");
	} # unless chdir
	@$cmdRef = ("$^X", "./scripts/AssignAlignType.pl", $align1_xmap, $align1_r_cmap, $align1_q_cmap, $assignAlignType_xmap, $assignAlignType_r_cmap, $assignAlignType_q_cmap, $T_cutoff, $max_overhang, $ngs_cmap, $bppAdjust, $breakPoint_txt, $breakPointFileShiftAmount);

	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	open(OUT, ">$outputDir/assignAlignType.log"); 	print OUT "$outResults";	close OUT;
	open(ERR, ">$outputDir/assignAlignType.errlog");	print ERR "$errResults";	close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
		my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
		printStatus(\*STATUS_FILE, "conflict detection", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "conflict detection", "Done", 0);
	}
	errCheck($outResults, $errResults, "ERROR:", "AssignAlignType complete in $span.", "ERROR: AssignAlignType cannot be completed.");

	#print number of BioNano and NGS contigs that have been flagged as conflicting
	my @bn; my @ngs;
	my @sticky = read_file($assignAlignType_xmap);
	for (my $i = 0; $i < scalar(@sticky); $i += 1)	{
		next if ($sticky[$i] =~ /^#/);
		my @s = split("\t", $sticky[$i]); 
		push (@bn, $s[1]);
		push (@ngs, $s[2]);
	}
	@bn = uniq @bn;
	@ngs = uniq @ngs;
	print scalar(@bn)." BNG contigs have been flagged as conflicting\n";
	print scalar(@ngs)." NGS contigs have been flagged as conflicting\n";

	#check AssignAlign outputs and Merge inputs contain at least one contig
	$assignAlign_r = abs_path("$opts{o}/assignAlignType/assignAlignType_r.cmap");
	(undef, $numContig, undef) = readCMap($assignAlign_r);
	if (!($numContig > 0)) {;
		print "WARNING: AssignAlign output file $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting.\n";
		warn "WARNING: AssignAlign output file $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting.\n";
	}
	$assignAlign_q = abs_path("$opts{o}/assignAlignType/assignAlignType_q.cmap");
	(undef, $numContig, undef) = readCMap($assignAlign_q);
	if (!($numContig > 0)) {
		print "WARNING: AssignAlign output file $assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting.\n";
		warn "WARNING: AssignAlign output file $assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting.\n";
	}

	# now perform cut at conflicts; run cut_conflicts.pl
	printStatus(\*STATUS_FILE, "conflict resolution", "Started", 0);
	$dtStartStage = DateTime->now;
	print "\nBeginning cut_conflicts...\n";
	$baseStage = "cut_conflicts";
	$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	$theParamValuesRef = extractConfig($commandStackRef, "window_size");	my $windowSize = $theParamValuesRef->[0];
	$theParamValuesRef = extractConfig($commandStackRef, "min_quality_score_threshold");	my $qScoreThreshold = $theParamValuesRef->[0];
	$theParamValuesRef = extractConfig($commandStackRef, "min_coverage_threshold");	my $covThreshold = $theParamValuesRef->[0];
	$baseStage = "";
	@$commandStackRef = ();

	my $oriBNFile = $bppAdjust;	
	my @theOriBNFileContent = split(/\//, $oriBNFile);	my $theOriBNFileName = $theOriBNFileContent[$#theOriBNFileContent];	my $theOutBNFileName = $theOriBNFileName;	$theOutBNFileName =~ s/\.cmap$/_cut.cmap/;
	my $oriSeqFile = $ngs_cmap;	
	my @theOriSeqFileContent = split(/\//, $oriSeqFile);	my $theOriSeqFileName = $theOriSeqFileContent[$#theOriSeqFileContent];	my $theOutSeqFileName = $theOriSeqFileName;	$theOutSeqFileName =~ s/\.cmap$/_cut.cmap/;
	my $conflictFile = $breakPoint_txt;
	my $theOutputDir = "$outputDir/cut_conflicts";	mkdir $theOutputDir if (! -e $theOutputDir);
	$cutSeqFile = abs_path("$opts{o}/assignAlignType/cut_conflicts/$theOutSeqFileName");
	$cutBNFile = abs_path("$opts{o}/assignAlignType/cut_conflicts/$theOutBNFileName");

	chdir $plDir or dieLog( "ERROR: Cannot change directory to $plDir: $!\n");
	# the rescue must be at least 2 times the max_overhang used in AssignAlignType
	@$cmdRef = ("$^X", "./scripts/cut_conflicts.pl", "-align1XmapFile", $align1_xmap, "-align1GMFile", $align1_q_cmap, "-align1SeqFile", $align1_r_cmap, "-maxOverhang", $max_overhang * 2, "-breakPointFileShiftAmount", $breakPointFileShiftAmount, "-oriGMFile", $oriBNFile, "-oriSeqFile", $oriSeqFile, "-conflictFile", $conflictFile, "-outDir", $theOutputDir, "-outGMFile", $cutBNFile, "-outSeqFile", $cutSeqFile, "-windowSize", $windowSize, "-qScoreThreshold", $qScoreThreshold, "-covThreshold", $covThreshold, "-refAligner", $refaligner);

	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	open(OUT, ">$theOutputDir/cut_conflicts.log"); 	print OUT $outResults;	close OUT;
	open(ERR, ">$theOutputDir/cut_conflicts.errlog");	print ERR $errResults;	close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
		my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
		printStatus(\*STATUS_FILE, "conflict resolution", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "conflict resolution", "Done", 0);
	} # if outResults
        errCheck($outResults, $errResults, "ERROR:", "cut_conflicts complete in $span.", "ERROR: cut_conflicts cannot be completed.");

        #copy conflict cut status file to one that has enzyme name
        my $new_cut_file = "conflicts_cut_status";
        foreach my $enzyme (@enzymes) {$new_cut_file .= "_$enzyme"};

        my $copy_cmd = "cp $theOutputDir/conflicts_cut_status.txt $theOutputDir/$new_cut_file\.txt";
        `$copy_cmd`;
       
  	# print the number of cut contigs
	(undef, $numContig, undef) = readCMap($cutBNFile);
	print "$numContig BNG contigs are found after the cut-conflict step\n";
	(undef, $numContig, undef) = readCMap($cutSeqFile);
	print "$numContig NGS contigs are found after the cut-conflict step\n";
} else	{
	# if the user defined a conflict status file manually, then the pipeline starts from this point on

	printStatus(\*STATUS_FILE, "manual conflict resolution", "Started", 0);	
	# check for the existence of some important files
	$bppAdjust = $modBNGfile;	$bppAdjust =~ s/_cmap$/.cmap/;	$bppAdjust = "$opts{o}/align0/$bppAdjust"; dieLog ("ERROR: BioNano bpp adjusted CMAP file does not exist in align0/: $bppAdjust\n") if (! -e $bppAdjust);
	$filename_key = basename($opts{n});	$filename_key =~ s/(\S+)\.\w+$/$1/; foreach my $enzyme (@enzymes) {$filename_key .= "_$enzyme"}; $filename_key = $filename_key."_".$minLength."kb_".$minLabels."labels_key.txt"; $filename_key = $opts{o}."/fa2cmap/".$filename_key; 
	if (! -e $filename_key)	{
		my $errMsg = "ERROR: fasta to cmap key file does not exist: $filename_key";
		printStatus(\*STATUS_FILE, "manual conflict resolution", $errMsg, 1);	dieLog("$errMsg\n");
	} # if no filename_key
	$filename = basename($opts{n});	$filename =~ s/(\S+)\.\w+$/$1/;	foreach my $enzyme (@enzymes) {$filename .= "_$enzyme"}; $filename = $filename."_".$minLength."kb_".$minLabels."labels.cmap"; $filename = $opts{o}."/fa2cmap/".$filename;
	$ngs_cmap = $filename;
	if (! -e $ngs_cmap)	{
		my $errMsg = "ERROR: sequence cmap file does not exist: $ngs_cmap";
		printStatus(\*STATUS_FILE, "manual conflict resolution", $errMsg, 1);	dieLog("$errMsg\n");
	} # if no ngs_cmap
	my $fastaName = $opts{n};	my @fastaNameContent = split(/\//, $fastaName);	$fastaName = $fastaNameContent[$#fastaNameContent];	$fastaName = "$opts{o}/fa2cmap/$fastaName";
	if (! -e $fastaName)	{
		my $errMsg = "ERROR: the sequence fasta file has not been copied to the output fa2cmap directory: $fastaName";
		printStatus(\*STATUS_FILE, "manual conflict resolution", $errMsg, 1);	dieLog("$errMsg\n");
	} # if no fastaName
		
	# check for the presence of align0 error bin file
	$readparameters = abs_path("$opts{o}/align1/align1.errbin"); 
	if (! -e $readparameters)	{
		my $errMsg = "ERROR: file does not exist: $readparameters";
		printStatus(\*STATUS_FILE, "manual conflict resolution", $errMsg, 1);	dieLog("$errMsg\n");
	} # if no readparameters
	
	# now perform cut at conflicts; run cut_conflicts.pl
	print "\nBeginning cut_conflicts based on input a conflict file...\n";

	my $oriBNFile = $bppAdjust;	
	my @theOriBNFileContent = split(/\//, $oriBNFile);	my $theOriBNFileName = $theOriBNFileContent[$#theOriBNFileContent];	my $theOutBNFileName = $theOriBNFileName;	$theOutBNFileName =~ s/\.cmap$/_cut.cmap/;
	my $oriSeqFile = $ngs_cmap;	
	my @theOriSeqFileContent = split(/\//, $oriSeqFile);	my $theOriSeqFileName = $theOriSeqFileContent[$#theOriSeqFileContent];	my $theOutSeqFileName = $theOriSeqFileName;	$theOutSeqFileName =~ s/\.cmap$/_cut.cmap/;
	my $theOutputDir = "$opts{o}/assignAlignType/cut_conflicts_M$modifyNum";	mkdir $theOutputDir if (! -e $theOutputDir);
	$cutSeqFile = abs_path("$theOutputDir/$theOutSeqFileName");
	$cutBNFile = abs_path("$theOutputDir/$theOutBNFileName");

	unless (chdir $plDir)	{
		my $errMsg = "ERROR: Cannot change directory to $plDir: $!\n";
		printStatus(\*STATUS_FILE, "manual conflict resolution", $errMsg, 1);	dieLog("$errMsg\n");
	} # unless chdir	
	@$cmdRef = ("$^X", "./scripts/cut_conflicts.pl", "-oriGMFile", $oriBNFile, "-oriSeqFile", $oriSeqFile, "-outDir", $theOutputDir, "-outGMFile", $cutBNFile, "-outSeqFile", $cutSeqFile, "-modBkptStatusFile", $opts{M}, "-refAligner", $refaligner);
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	open(OUT, ">$theOutputDir/cut_conflicts.log");	print OUT $outResults;	close OUT;
	open(ERR, ">$theOutputDir/cut_conflicts.errlog");	print ERR $errResults;	close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
		my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
		printStatus(\*STATUS_FILE, "manual conflict resolution", $errMsg, 1);	
	} else	{
		printStatus(\*STATUS_FILE, "manual conflict resolution", "Done", 0);
	} # if outResults
	errCheck($outResults, $errResults, "ERROR:", "cut_conflicts complete in $span.", "ERROR: cut_conflicts cannot be completed.");	
        
	#copy conflict cut status file to one that has enzyme name
        my $new_cut_file = "conflicts_cut_status";
        foreach my $enzyme (@enzymes) {$new_cut_file .= "_$enzyme"};

        my $copy_cmd = "cp $theOutputDir/conflicts_cut_status.txt $theOutputDir/$new_cut_file\.txt";
`$copy_cmd`;

	# print the number of cut contigs
	(undef, $numContig, undef) = readCMap($cutBNFile);
	print "$numContig BNG contigs are found after the cut-conflict step\n";
	(undef, $numContig, undef) = readCMap($cutSeqFile);
	print "$numContig NGS contigs are found after the cut-conflict step\n";
} # if defined a conflict status file 


#we only run the pre-processing stage, do not merge contigs yet

if(defined($opts{S})){
    exit 0;
}



# run MergeNGS_BN
printStatus(\*STATUS_FILE, "pairwise merge", "Started", 0);
$dtStartStage = DateTime->now;
print "\nBeginning MergeNGS_BN...\n";
$outputDir = $opts{o}."/mergeNGS_BN";	$outputDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");
eval { mkpath($outputDir) };
if ($@) {
	print "Couldn't create $outputDir: $@"; }

$baseStage = "mergeNGS_BN";
$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
$theParamValuesRef = extractConfig($commandStackRef, "merge_Tvalue");	my $merge_Tvalue = $theParamValuesRef->[0];
$theParamValuesRef = extractConfig($commandStackRef, "id_shift");	my $id_shift = $theParamValuesRef->[0];
$baseStage = "";
@$commandStackRef = ();

### but first, we need to define a variable that stores the input files to the merge step
my ($mergeBN, $mergeNgs) = ("", "");
if (defined($opts{M}))	{
	# user defined a conflict file
	($mergeBN, $mergeNgs) = ($cutBNFile, $cutSeqFile);
	print "*Using user-defined conflict-cut BioNano and sequence CMAP files*\n";
} else	{
	if ($opts{B} == 1 && $opts{N} == 1)	{
		($mergeBN, $mergeNgs) = ($bppAdjust, $ngs_cmap);
		print "*Using all BioNano and all sequence CMAP*\n";
	} elsif ($opts{B} == 1 && $opts{N} == 2)	{
		($mergeBN, $mergeNgs) = ($bppAdjust, $cutSeqFile);
		print "*Using all BioNano and conflict-cut sequence CMAP*\n";
	} elsif ($opts{B} == 1 && $opts{N} == 3)	{
		($mergeBN, $mergeNgs) = ($bppAdjust, $assignAlign_r);
		print "*Using all BioNano and non-conflicting-only sequence CMAP*\n";
	} elsif ($opts{B} == 2 && $opts{N} == 1)	{
		($mergeBN, $mergeNgs) = ($cutBNFile, $ngs_cmap);
		print "*Using conflict-cut BioNano and all sequence CMAP*\n";
	} elsif ($opts{B} == 2 && $opts{N} == 2)	{
		($mergeBN, $mergeNgs) = ($cutBNFile, $cutSeqFile);
		print "*Using conflict-cut BioNano and conflict-cut sequence CMAP*\n";
	} elsif ($opts{B} == 2 && $opts{N} == 3)	{
		($mergeBN, $mergeNgs) = ($cutBNFile, $assignAlign_r);
		print "*Using conflict-cut BioNano and non-conflicting-only sequence CMAP*\n";
	} elsif ($opts{B} == 3 && $opts{N} == 1)	{
		($mergeBN, $mergeNgs) = ($assignAlign_q, $ngs_cmap);
		print "*Using non-conflicting-only BioNano and all sequence CMAP*\n";
	} elsif ($opts{B} == 3 && $opts{N} == 2)	{
		($mergeBN, $mergeNgs) = ($assignAlign_q, $cutSeqFile);
		print "*Using non-conflicting-only BioNano and conflict-cut sequence CMAP*\n";
	} else	{
		($mergeBN, $mergeNgs) = ($assignAlign_q, $assignAlign_r);
		print "*Using non-conflicting-only BioNano and non-conflicting-only sequence CMAP*\n";
	} # if cut
} # if a modified conflict file is specified

# make sure that the input files contain at least one contig
(undef, $numContig, undef) = readCMap($mergeNgs);
if ($numContig < 1)	{
	my $errMsg = "ERROR: before merge, sequence file $mergeNgs does not contain any contigs! You can try using all (no filter) sequence option.";
	printStatus(\*STATUS_FILE, "pairwise merge", $errMsg, 1);
	dieLog("$errMsg\n\n");
} # if numContig
(undef, $numContig, undef) = readCMap($mergeBN);
if ($numContig < 1)	{
	my $errMsg = "ERROR: before merge, BioNano file $mergeBN does not contain any contigs! You can try using all (no filter) BioNano option.";
	printStatus(\*STATUS_FILE, "pairwise merge", $errMsg, 1);
	dieLog("$errMsg\n\n");
} # if numContig

# we need to have the id shift value to be bigger than the largest ngs cmap id value
my ($maxSeqId, undef) = readCmapIdLength($mergeNgs);
$id_shift = int($id_shift);
$id_shift = ($id_shift <= 0) ? (100000) : ($id_shift);	# properly set the value
while ($id_shift < $maxSeqId)	{
	$id_shift = $id_shift * 10;
} # while id_shift

unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "pairwise merge", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/MergeNGS_BN.pl", "-outputDir", $outputDir, "-refaligner", $refaligner, "-logFile", "mergeNGS_BN_script.log", "-ngs_cmap_fn", $mergeNgs, "-bng_cmap_fn", $mergeBN, "-id_shift", $id_shift, "-readparameters", $readparameters, "-xmlFile", $opts{c});

print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
open(OUT, ">$outputDir/MergeNGS_BN.log");	print OUT "$outResults";	close OUT;
open(ERR, ">$outputDir/MergeNGS_BN.errlog");	print ERR "$errResults";	close ERR;

$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
	printStatus(\*STATUS_FILE, "pairwise merge", $errMsg, 1);
} else	{
	printStatus(\*STATUS_FILE, "pairwise merge", "Done", 0);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "MergeNGS_BN complete in $span.", "ERROR: MergeNGS_BN cannot be completed");


#find NGS contigs that went into hybrid assembly and make new NGS cmap
# $outputDir is pointing to $opts{o}/mergeNGS_BN(_M$modifyNum)
printStatus(\*STATUS_FILE, "extract sequences in hybrid scaffold pairmerge step", "Started", 0);
my $final_outputDir = "$opts{o}/align_final";	$final_outputDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");	
eval { mkpath($final_outputDir) };
if ($@) {
	print "Couldn't create $final_outputDir: $@"; }
my $hybridScaffoldFileName = "step2.hybrid.cmap";
copy "$outputDir/$hybridScaffoldFileName", "$final_outputDir/$hybridScaffoldFileName" ;

$dtStartStage = DateTime->now;
print "\nBeginning extraction of used and not used NGS contigs in hybrid scaffold pairmerge step...\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "extract used sequences in hybrid scaffold pairmerge step", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
if (defined($opts{M}) || $opts{N} == 2)	{
	# if the user indicated to cut
	@$cmdRef = ($^X, "scripts/find_used_not_used_ngs.pl", $outputDir, $id_shift, $cutSeqFile, $refaligner);
} else	{
	# if the user did not indicate cut, take the original input NGS 
	@$cmdRef = ($^X, "scripts/find_used_not_used_ngs.pl", $outputDir, $id_shift, $ngs_cmap, $refaligner);
} # if optsN
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
open(OUT, ">$outputDir/find_used_not_used_ngs.log"); print OUT $outResults;  close OUT;
open(ERR, ">$outputDir/find_used_not_used_ngs.errlog"); print ERR $errResults;	close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
	printStatus(\*STATUS_FILE, "extract used sequences in hybrid scaffold pairmerge step", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "Extraction of NGS contigs complete in $span.", "ERROR: Extraction of NGS contigs cannot be completed.");

copy "$outputDir/all_used_ngs.cmap", "$final_outputDir/filtered_NGS.cmap";
copy "$outputDir/all_not_used_ngs.cmap", "$final_outputDir/filtered_not_used_NGS.cmap";


#align the NGS contigs to hybrid assembly
$dtStartStage = DateTime->now;
print "\nBeginning alignment of NGS cmap to Hybrid CMAP...\n";

my $prefixOrig = $BNGfile."_".$NGSfile."_HYBRID_SCAFFOLD";
my $prefix = $modBNGfile."_".$NGSfile."_NGScontigs_HYBRID_SCAFFOLD";

# $final_outputDir is pointing to $opts{o}/align_final(_M$modifyNum)
my $filteredNGS_file = $mergeNgs;
my $hybrid_cmap = "$final_outputDir/$hybridScaffoldFileName";

unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change to $plDir: $!";
	printStatus(\*STATUS_FILE, "align sequences to hybrid scaffolds", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ($^X, "./scripts/align_final_seq_two_passes.pl", "-outDir", $final_outputDir, "-outFilePrefix", $prefix, "-refAligner", $refaligner, "-xmlFile", $opts{c}, "-seqFile", $mergeNgs, "-scaffoldFile", $hybrid_cmap, "-logFile", "$final_outputDir/align_final_seq_script_log.txt");

print "Running command: ", (join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
open(OUT, ">$final_outputDir/align_final_seq.log"); print OUT $outResults; close OUT;
open(ERR, ">$final_outputDir/align_final_seq.errlog"); print ERR $errResults; close ERR;
# error check
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
if ($outResults =~/ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
	printStatus(\*STATUS_FILE, "align sequences to hybrid scaffolds", $errMsg, 1);
} else	{
	printStatus(\*STATUS_FILE, "align sequences to hybrid scaffolds", "Done", 0);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "align sequences to hybrid scaffolds completed  in $span.", "ERROR: align sequences to hybrid scaffolds cannot be completed");


#Align original BNG contigs to Hybrid map 
#find BNG contigs that went into hybrid assembly and make new BNG cmap
printStatus(\*STATUS_FILE, "align used BioNano maps to hybrid scaffolds pairmerge step", "Started", 0);
$dtStartStage = DateTime->now;
# $outputDir is pointing to $opts{o}/mergeNGS_BN(_M$modifyNum)
print "\nBeginning extraction of used and not used BNG contigs in hybrid scaffold pairmerge step...\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
if (defined($opts{M}) || $opts{B} == 2)	{
	# the user indicated to cut
	@$cmdRef = ($^X, "scripts/find_used_not_used_bn.pl", $outputDir, $id_shift, $cutBNFile, $refaligner);
} else	{
	# the user indicated not to cut (put in the original bpp-adjusted BioNano genome maps)
	@$cmdRef = ($^X, "scripts/find_used_not_used_bn.pl", $outputDir, $id_shift, $bppAdjust, $refaligner);
} # if optsB
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
open(OUT, ">$outputDir/find_used_not_used_bn.log"); print OUT $outResults;	close OUT;
open(ERR, ">$outputDir/find_used_not_used_bn.errlog");	print ERR $errResults;	close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
	printStatus(\*STATUS_FILE, "extract used BioNano maps in hybrid scaffold pairmerge step", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "Extraction of BNG contigs complete in $span.", "ERROR: Extraction of BNG contigs cannot be completed.");

copy "$outputDir/all_used_bng.cmap", "$final_outputDir/filtered_BNG.cmap";
copy "$outputDir/all_not_used_bng.cmap", "$final_outputDir/filtered_not_used_BNG.cmap";

#align the BNG genome maps in hybrid back to hybrid assembly
$dtStartStage = DateTime->now;
print "\nBeginning alignment of BNG cmap to Hybrid CMAP...\n";

my $prefixBNG = $modBNGfile."_".$NGSfile."_BNGcontigs_HYBRID_SCAFFOLD";

my $filteredBNG_file = "$final_outputDir/filtered_BNG.cmap";
$hybrid_cmap = "$final_outputDir/$hybridScaffoldFileName";

unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change to $plDir: $!";
	printStatus(\*STATUS_FILE, "extract BioNano maps in hybrid scaffolds", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ($^X, "./scripts/align_final_bng.pl", "-outDir", $final_outputDir, "-outFilePrefix", $prefixBNG, "-refAligner", $refaligner, "-xmlFile", $opts{c}, "-bngFile", $mergeBN, "-scaffoldFile", $hybrid_cmap, "-logFile", "$final_outputDir/align_final_bng_log.txt");

print "Running command: ", (join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
open(OUT, ">$final_outputDir/align_final_bng.log");	print OUT $outResults;	close OUT;
open(ERR, ">$final_outputDir/align_final_bng.errlog");	print ERR $errResults;	close ERR;
# error check
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
if ($outResults =~/ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
	printStatus(\*STATUS_FILE, "align Bionano genome maps to hybrid scaffolds", $errMsg, 1);
} else	{
	printStatus(\*STATUS_FILE, "align Bionano genome maps to hybrid scaffolds", "Done", 0);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "align Bionano genome maps to hybrid scaffolds completed  in $span.", "ERROR: align Bionano genome maps to hybrid scaffolds cannot be completed");

#merge hybrid cmap with all NGS not participated in scaffolding process
# $outputDir is pointing to $opts{o}/mergeNGS_BN(_M$modifyNum)
# $final_outputDir is pointing to $opts{o}/align_final(_M$modifyNum)
printStatus(\*STATUS_FILE, "combine hybrid scaffolds with left over sequences", "Started", 0);
$dtStartStage = DateTime->now;
print "\nMerging Hybrid CMAP with NGS not participated in the hybrid scaffold...\n";
my $notUsedNGSCmap = "$outputDir/all_not_used_ngs.cmap";
my $hybridNotUsedNGS = "HYBRID_SCAFFOLD_notUsedNGS_merged.cmap";
unless (chdir $final_outputDir)	{
	my $errMsg = "ERROR: Cannot change to $final_outputDir: $!";
	printStatus(\*STATUS_FILE, "combine hybrid scaffolds with left over sequences", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ($refaligner, "-f", "-i", $hybrid_cmap, "-i", $notUsedNGSCmap, "-o", "HYBRID_SCAFFOLD_notUsedNGS_merged", "-merge", "-minsites", 0, "-stdout", "-stderr");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
if (errCheckRefAligner("HYBRID_SCAFFOLD_notUsedNGS_merged.stdout", "END of output", "Merging Hybrid CMAP with naive NGS CMAP complete in $span.", "ERROR: Merging Hybrid CMAP with naive NGS CMAP cannot be completed.", 0) != 1)	{
	my $errMsg = "ERROR: Merging Hybrid CMAP with naive NGS CMAP cannot be completed.";
	printStatus(\*STATUS_FILE, "combine hybrid scaffolds with left over sequences", $errMsg, 1);
} else	{
	printStatus(\*STATUS_FILE, "combine hybrid scaffolds with left over sequences", "Done", 0);
} # if refaligner check
errCheckRefAligner("HYBRID_SCAFFOLD_notUsedNGS_merged.stdout", "END of output", "Merging Hybrid CMAP with naive NGS CMAP complete in $span.", "ERROR: Merging Hybrid CMAP with naive NGS CMAP cannot be completed.", 1);

#merge hybrid cmap with all BNG not participated in scaffolding process
printStatus(\*STATUS_FILE, "combine hybrid scaffolds with left over BioNano maps", "Started", 0);
$dtStartStage = DateTime->now;
print "\nMerging Hybrid CMAP with BNG CMAP not participated in the hybrid scaffold...\n";
my $notUsedBNGCmap = "$outputDir/all_not_used_bng.cmap";
my $hybridNotUsedBNG = "HYBRID_SCAFFOLD_notUsedBNG_merged.cmap";
unless (chdir $final_outputDir)	{
	my $errMsg = "ERROR: Cannot change to $final_outputDir: $!";
	printStatus(\*STATUS_FILE, "combine hybrid scaffolds with left over BioNano maps", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ($refaligner, "-f", "-i", $hybrid_cmap, "-i", $notUsedBNGCmap, "-o", "HYBRID_SCAFFOLD_notUsedBNG_merged", "-merge", "-stdout", "-stderr");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
if (errCheckRefAligner("HYBRID_SCAFFOLD_notUsedBNG_merged.stdout", "END of output", "Merging Hybrid CMAP with naive BNG CMAP complete in $span.", "ERROR: Merging Hybrid CMAP with naive BNG CMAP cannot be completed.", 0) != 1)	{
	my $errMsg = "ERROR: Merging Hybrid CMAP with naive BNG CMAP cannot be completed.";
	printStatus(\*STATUS_FILE, "combine hybrid scaffolds with left over BioNano maps", $errMsg, 1);
} else	{
	printStatus(\*STATUS_FILE, "combine hybrid scaffolds with left over BioNano maps", "Done", 0);
} # if refaligner check
errCheckRefAligner("HYBRID_SCAFFOLD_notUsedBNG_merged.stdout", "END of output", "Merging Hybrid CMAP with naive BNG CMAP complete in $span.", "ERROR: Merging Hybrid CMAP with naive BNG CMAP cannot be completed.", 1);

# for TGH, we need a pre-pairmerge sequence (reference) to pre-pairmerge BNG alignment (query)
if (defined($opts{t}))	{
	printStatus(\*STATUS_FILE, "pre-pairmerge sequence alignment to pre-pairmerge BioNano maps", "Started", 0);
	$dtStartStage = DateTime->now;
	print "\nBeginning alignment between pre-pairmerge sequence and pre-pairmerge BioNano maps...\n";

	unless (chdir $final_outputDir)	{
		my $errMsg = "ERROR: Cannot change to $final_outputDir: $!\n";
		printStatus(\*STATUS_FILE, "pre-pairmerge sequence alignment and pre-pairmerge BioNano maps", $errMsg, 1);
		dieLog("$errMsg\n");
	} # unless chdir
	my $theTempOutputPrefix = "$final_outputDir/".(join("_", "E", @enzymes, "Q", "prePairmerge_BNG", "A", "prePairmerge_NGS")); 
	
	@$cmdRef = ($refaligner, "-ref", $mergeNgs, "-i", $mergeBN, "-o", $theTempOutputPrefix, "-stdout", "-stderr", "-XmapStatWrite", "$theTempOutputPrefix.stats");
	$baseStage = "global";
	$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	$baseStage = "align1";
	$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	$cmdRef = copyArray($cmdRef, $commandStackRef); # the rest of the parameters
	#not using BestRef 1
	my $cmdRefStr = join(" ", @$cmdRef);
	$cmdRefStr =~ s/-BestRef 1/-BestRef 0/;
	my @ary = split(/\s/, $cmdRefStr);
	$cmdRef = \@ary;

	$baseStage = "";
	@$commandStackRef = ();

	print "Running command: ". (join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	# error check
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	if (errCheckRefAligner("$theTempOutputPrefix.stdout", "END of output", "Alignment between pre-pairmerge sequence and pre-pairmerge BioNano maps complete in $span.", "ERROR: alignment between pre-pairmerge sequence and pre-pairmerge BioNano maps cannot be completed.", 0) != 1)	{
		my $errMsg = "ERROR: alignment between pre-pairmerge sequence and pre-pairmerge BioNano maps.";
		printStatus(\*STATUS_FILE, "pre-pairmerge sequence alignment to pre-pairmerge BioNano maps", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "pre-pairmerge sequence alignment to pre-pairmerge BioNano maps", "Done", 0);
	} # if refaligner check
	errCheckRefAligner("$theTempOutputPrefix.stdout", "END of output", "Alignment between pre-pairmerge sequence and pre-pairmerge BioNano maps complete in $span.", "ERROR: alignment between pre-pairmerge sequence and pre-pairmerge BioNano maps cannot be completed.", 1);
} # if defined opts t

#generate AGP and FASTA file
# $outputDir is pointing to $opts{o}/mergeNGS_BN(_M$modifyNum)
# $final_outputDir is pointing to $opts{o}/align_final(_M$modifyNum)
printStatus(\*STATUS_FILE, "AGP and FASTA generation", "Started", 0);
$dtStartStage = DateTime->now;
my $agpFastaDir = "$opts{o}/agp_fasta";	$agpFastaDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");	
mkdir $agpFastaDir if (! -e $agpFastaDir);
my $fastaName = $opts{n};       my @fastaNameContent = split(/\//, $fastaName); $fastaName = $fastaNameContent[$#fastaNameContent];     $fastaName = "$opts{o}/fa2cmap/$fastaName";
print "\nBeginnning construction of AGP and FASTA file of the scaffolded and unscaffolded sequences...\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change to $plDir: $!";
	printStatus(\*STATUS_FILE, "AGP and FASTA generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
my $enzymeString = "";
foreach my $enzyme (@enzymes)	{
	$enzymeString .= "$enzyme $channelNumber ";	
} # foreach enzyme
$enzymeString =~ s/\s$//;
if (defined($opts{M}))	{
	# the user indicated to use manually cut results
	# enzyme argument here is somewhat weird
	my $cutNGSCoordFile = "$opts{o}/assignAlignType/cut_conflicts_M$modifyNum/auto_cut_NGS_coord_translation.txt"; 
	if (! -e $cutNGSCoordFile)	{
		my $errMsg = "ERROR: file indicating NGS cut coordinates does not exists: $cutNGSCoordFile";
		printStatus(\*STATUS_FILE, "AGP and FASTA generation", $errMsg, 1);
		dieLog("$errMsg\n");
	} # if no cutNGSCoordFile
	@$cmdRef = ($^X, "scripts/ExportAGP.pl", "-i", "$final_outputDir/$prefix.xmap", "-c", "$final_outputDir/$prefix"."_r.cmap", "-o", $agpFastaDir, "-m", $filename_key, "-s", $fastaName, "-e", $enzymeString, "-r", $cutNGSCoordFile);
} elsif ($opts{N} == 2)	{
	# no manual cut, auto cut was done instead
	my $cutNGSCoordFile = "$opts{o}/assignAlignType/cut_conflicts/auto_cut_NGS_coord_translation.txt";	
	if (! -e $cutNGSCoordFile)	{
		my $errMsg = "ERROR: file indicating NGS cut coordinates does not exists: $cutNGSCoordFile";
		printStatus(\*STATUS_FILE, "AGP and FASTA generation", $errMsg, 1);
		dieLog("$errMsg\n");
	} # if no cutNGSCoordFile
	@$cmdRef = ($^X, "scripts/ExportAGP.pl", "-i", "$final_outputDir/$prefix.xmap", "-c", "$final_outputDir/$prefix"."_r.cmap", "-o", $agpFastaDir, "-m", $filename_key, "-s", $fastaName, "-e", $enzymeString, "-r", $cutNGSCoordFile);
} else	{
	# the user indicated no cut to the results
	# enzyme argument here is somewhat weird
	@$cmdRef = ($^X, "scripts/ExportAGP.pl", "-i", "$final_outputDir/$prefix.xmap", "-c", "$final_outputDir/$prefix"."_r.cmap", "-o", $agpFastaDir, "-m", $filename_key, "-s", $fastaName, "-e", $enzymeString);	
} # if opts
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
open(OUT, ">$agpFastaDir/xmap2agp.log") or dieLog ("ERROR: Cannot write to $agpFastaDir/xmap2agp.log: $!\n");	print OUT $outResults;	close OUT;
open(ERR, ">$agpFastaDir/xmap2agp.errlog") or dieLog ("ERROR: Cannot write to $agpFastaDir/xmap2agp.errlog: $!\n"); print ERR $errResults; close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
	printStatus(\*STATUS_FILE, "AGP and FASTA generation", $errMsg, 1);
} else	{
	printStatus(\*STATUS_FILE, "AGP and FASTA generation", "Done", 0);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "AGP and FASTA generation complete in $span.", "ERROR: AGP and FASTA generation cannot be completed.");


#Calculate some stats about CMAPs
printStatus(\*STATUS_FILE, "statistical report generation", "Started", 0);
$dtStartStage = DateTime->now;
print "\nCalculating statistics...\n\n";

### hybrid_scaffold_informatics_report.txt (a sub set of the log file) ###
open(INFO_REPORT, ">$hybridScaffoldFinalResultsDir/hybrid_scaffold_informatics_report.txt") or dieLog "ERROR: Cannot write to informatics report $hybridScaffoldFinalResultsDir/hybrid_scaffold_informatics_report.txt: $!\n";

print "Original BioNano Genome Map statistics:\n";
print INFO_REPORT "Original BioNano Genome Map statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $opts{b});
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# errorCheck
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for the original BioNano CMAP file";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for the original BioNano CMAP file");
# print statistics
print "$outResults\n";	print INFO_REPORT "$outResults\n";

print "Bpp-adjusted BioNano Genome Map statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $bppAdjust);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# errorCheck
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for the bpp-adjusted BioNano CMAP file";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for the bpp-adjusted BioNano CMAP file");
# print statistics
print "$outResults\n";

print "Original NGS sequences statistics:\n";
print INFO_REPORT "Original NGS sequences statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $ngs_cmap);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# errorCheck
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for original NGS sequences";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for original NGS sequences");
# print statistics
print "$outResults\n";	print INFO_REPORT "$outResults\n";

if (defined ($opts{M}) || $opts{N} == 2)	{
	# conflict resolution was done
	print "Conflict Resolution from BNG-NGS alignment:\n";
	print INFO_REPORT "Conflict Resolution from BNG-NGS alignment:\n";
	unless (chdir $plDir)	{
		my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
		printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
		dieLog("$errMsg\n");
	} # unless chdir
	my $conflictsFile = (defined($opts{M})) ? ("$opts{o}/assignAlignType/cut_conflicts_M$modifyNum/conflicts_cut_status.txt") : ("$opts{o}/assignAlignType/cut_conflicts/conflicts_cut_status.txt");
	@$cmdRef = ("$^X", "./scripts/calc_conflicts_cut_stats.pl", $conflictsFile);
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	# errorCheck
	if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
		my $errMsg = "ERROR: cannot calculate statistics for conflict cut status file";
		printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	} # if outResults
	errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for conflict cut status file");
	print "$outResults\n";	print INFO_REPORT "$outResults\n";
} # if defined

print "Before merge: BioNano Genome Map statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $mergeBN);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# errorCheck
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for BioNano Genome Map right before the Merge";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for BioNano Genome Map right before the Merge");
# print statistics
print "$outResults\n";

print "Before merge: NGS sequences statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $mergeNgs);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# errorCheck
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for NGS sequences right before the Merge";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for NGS sequences right before the Merge");
# print statistics
print "$outResults\n";	

print "BNG Genome Map in hybrid scaffold statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", "$final_outputDir/$prefixBNG"."_q.cmap");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# errorCheck
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for BNG Genome Map aligned to hybrid scaffold";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for BNG Genome Map aligned to hybrid scaffold");
# print statistics
print "$outResults\n";

print "NGS sequences in hybrid scaffold (CMAP) statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", "$final_outputDir/$prefix"."_q.cmap");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for NGS sequences aligned to hybrid scaffold (CMAP)";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);	
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics NGS sequences aligned to hybrid scaffold (CMAP)");
# print statistics
print "$outResults\n";

print "Hybrid scaffold Map statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", "$final_outputDir/$hybridScaffoldFileName");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for hybrid scaffold Map";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for hybrid scaffold Map");
# print statistics
print "$outResults\n";

print "NGS FASTA sequence in hybrid scaffold statistics:\n";
print INFO_REPORT "NGS FASTA sequence in hybrid scaffold statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir

if (defined($opts{M}) || $opts{N} == 2)	{
	my ($postCutSeqDir, $postCutSeqFile) = ("$opts{o}/agp_fasta", $opts{n});
	$postCutSeqDir .= "_M$modifyNum" if (defined($opts{M}));
	my @postCutSeqFileContent = split(/\//, $postCutSeqFile);
	$postCutSeqFile = $postCutSeqFileContent[$#postCutSeqFileContent];
	$postCutSeqFile .= ".cut.fasta";

	@$cmdRef = ($^X, "scripts/calc_scaffolded_seq_fasta_stats.pl", "$postCutSeqDir/$postCutSeqFile", "$agpFastaDir/$prefix"."_NOT_SCAFFOLDED.fasta");
} else	{
	@$cmdRef = ($^X, "scripts/calc_scaffolded_seq_fasta_stats.pl", $opts{n}, "$agpFastaDir/$prefix"."_NOT_SCAFFOLDED.fasta");
}
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for NGS FASTA sequence in hybrid scaffold";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for NGS FASTA sequence in hybrid scaffold");
# print statistics
print "$outResults\n";	print INFO_REPORT "$outResults\n";

print "Hybrid scaffold FASTA statistics:\n";
print INFO_REPORT "Hybrid scaffold FASTA statistics:\n";
unless (chdir $plDir) {
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!\n";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ($^X, "scripts/calc_fasta_stats.pl", "$agpFastaDir/$prefix.fasta");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/) {
	my $errMsg = "ERROR: cannot calculate statistics for hybrid scaffold FASTA statistics";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for hybrid scaffold FASTA");
# print statistics
print "$outResults\n";	print INFO_REPORT "$outResults\n";

print "Hybrid scaffold FASTA plus not scaffolded NGS FASTA statistics:\n";
print INFO_REPORT "Hybrid scaffold FASTA plus not scaffolded NGS FASTA statistics:\n";
unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!\n";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ($^X, "scripts/calc_hybrid_scaffold_not_scaffolded_seq_fasta_stats.pl", "$agpFastaDir/$prefix.fasta", "$agpFastaDir/$prefix"."_NOT_SCAFFOLDED.fasta");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
	my $errMsg = "ERROR: cannot calculate statistics for hybrid scaffold FASTA plus not scaffolded NGS FASTA";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for hybrid scaffold FASTA plus not scaffolded NGS FASTA");
# print statistics
print "$outResults\n";	print INFO_REPORT "$outResults\n";

#calculate XMAP stats and output to .stats file

$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
print "\nCalculating CMAP statistics complete in $span.\n";

$dtStartStage = DateTime->now;
print "\nCalculating XMAP statistics...\n";

unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $final_outputDir, "$prefix.xmap", $hybridScaffoldFileName, "$prefix.xmap.stats");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
open(OUT, ">$final_outputDir/calc_xmap_stats.log");	print OUT $outResults; close OUT;
open(ERR, ">$final_outputDir/calc_xmap_stats.errlog");	print ERR $errResults;	close ERR;
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/) {
	my $errMsg = "ERROR: cannot calculate statistics for xmap $prefix.xmap";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefix.xmap");
# print statistics
print "$outResults\n";

unless (chdir $plDir)	{
	my $errMsg = "ERROR: Cannot change directory to $plDir: $!";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
	dieLog("$errMsg\n");
} # unless chdir
@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $final_outputDir, "$prefixBNG.xmap", $hybridScaffoldFileName, "$prefixBNG.xmap.stats");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
@$cmdRef = ();
open(OUT, ">$final_outputDir/calc_xmap_stats.log");	print OUT $outResults; close OUT;
open(ERR, ">$final_outputDir/calc_xmap_stats.errlog");	print ERR $errResults;	close ERR;
# error check
if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/) {
	my $errMsg = "ERROR: cannot calculate statistics for xmap $prefixBNG.xmap";
	printStatus(\*STATUS_FILE, "statistical report generation", $errMsg, 1);
} else	{
	printStatus(\*STATUS_FILE, "statistical report generation", "Done", ((defined $opts{x}) ? (0) : (1)));	# also close the status file, if molecule-alignment is not to be executed
} # if outResults
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefixBNG.xmap");
# print statistics
print "$outResults\n";

$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
print "XMAP statistics calculation complete in $span.\n";

### now copy the important results files to the final output directory ###
$hybridScaffoldFileName = $modBNGfile."_$NGSfile"."_HYBRID_SCAFFOLD.cmap";
copy "$final_outputDir/step2.hybrid.cmap", "$hybridScaffoldFinalResultsDir/$hybridScaffoldFileName";	# hybrid scaffold CMAP file

# prefix is  $prefix = $modBNGfile."_".$NGSfile."_NGScontigs_HYBRID_SCAFFOLD";
copy "$final_outputDir/$prefix.xmap", "$hybridScaffoldFinalResultsDir";
copy "$final_outputDir/$prefix"."_q.cmap", "$hybridScaffoldFinalResultsDir";
copy "$final_outputDir/$prefix"."_r.cmap", "$hybridScaffoldFinalResultsDir";

# prefixBNG is $prefixBNG = $modBNGfile."_".$NGSfile."_BNGcontigs_HYBRID_SCAFFOLD";
copy "$final_outputDir/$prefixBNG.xmap", "$hybridScaffoldFinalResultsDir";
copy "$final_outputDir/$prefixBNG"."_q.cmap", "$hybridScaffoldFinalResultsDir";
copy "$final_outputDir/$prefixBNG"."_r.cmap", "$hybridScaffoldFinalResultsDir";

# copy the cut status text file, only if users specified cut was to be done (i.e. N 2 or B 2 or M)
if (defined($opts{M}) || $opts{N} == 2 || $opts{B} == 2)	{
	my $sourceDir = "$opts{o}/assignAlignType/cut_conflicts";	$sourceDir .= "_M$modifyNum" if (defined($opts{M}));
	copy "$sourceDir/conflicts_cut_status.txt", "$hybridScaffoldFinalResultsDir";
	copy "$sourceDir/auto_cut_NGS_coord_translation.txt", "$hybridScaffoldFinalResultsDir";
	copy "$sourceDir/auto_cut_BN_coord_translation.txt", "$hybridScaffoldFinalResultsDir";

	my $extraCutFile = "conflicts_cut_status";
	foreach my $enzyme (@enzymes)	{$extraCutFile .= "_$enzyme"};
	$extraCutFile .= ".txt";
	copy "$sourceDir/$extraCutFile", "$hybridScaffoldFinalResultsDir";
} # if to copy

# copy the bed files indicating whether the new fragments were cut or not
if (defined($opts{M}) || $opts{N} == 2)	{
	my $sourceDir = "$opts{o}/assignAlignType/cut_conflicts";	$sourceDir .= "_M$modifyNum" if (defined($opts{M}));
	copy "$sourceDir/ngs_pre_cut_annotations.bed", "$hybridScaffoldFinalResultsDir";
} # if to copy
if (defined($opts{M}) || $opts{B} == 2)	{
	my $sourceDir = "$opts{o}/assignAlignType/cut_conflicts";       $sourceDir .= "_M$modifyNum" if (defined($opts{M}));
	copy "$sourceDir/bn_pre_cut_projected_ngs_coord_annotations.bed", "$hybridScaffoldFinalResultsDir";
} # if to copy

# copy the AGP/FASTA/trim files 
# output a manifest file for NCBI export
open(OUT_MANIFEST, ">$hybridScaffoldFinalResultsDir/ncbi_manifest.txt") or dieLog "ERROR: Cannot write the NCBI manifest file $hybridScaffoldFinalResultsDir/ncbi_manifest.txt: $!\n";
my $theSourceDir = "$opts{o}/agp_fasta";	$theSourceDir .= "_M$modifyNum" if (defined($opts{M}));
copy "$theSourceDir/$prefix.agp", "$hybridScaffoldFinalResultsDir";			# AGP file
copy "$theSourceDir/$prefix.gap", "$hybridScaffoldFinalResultsDir";			# gap file
copy "$theSourceDir/$prefix.fasta", "$hybridScaffoldFinalResultsDir";			# output FASTA with inter-contig nick sites
copy "$theSourceDir/$prefix"."_NCBI.fasta", "$hybridScaffoldFinalResultsDir";		# output FASTA without inter-contig nick sites
copy "$theSourceDir/$prefix"."_NOT_SCAFFOLDED.fasta", "$hybridScaffoldFinalResultsDir";	# output FASTA of non-scaffolded sequence contigs
copy "$theSourceDir/$prefix"."_trimHeadTailGap.coord", "$hybridScaffoldFinalResultsDir"; # .coord file

# pre-scaffold/post-conflict resolution FASTA	
if (defined($opts{M}) || $opts{N} == 2) {
                copy "$theSourceDir/$fastaNameContent[$#fastaNameContent].cut.fasta", "$hybridScaffoldFinalResultsDir";
		print OUT_MANIFEST "$fastaNameContent[$#fastaNameContent].cut.fasta\n";
} else     {
                copy "$opts{n}", "$hybridScaffoldFinalResultsDir";
		my @seqFileNameContent = split(/\//, $opts{n});
		print OUT_MANIFEST "$seqFileNameContent[$#seqFileNameContent]\n";
}
print OUT_MANIFEST "$prefix.agp\n";
print OUT_MANIFEST "$prefix"."_trimHeadTailGap.coord\n";
print OUT_MANIFEST "$hybridScaffoldFileName\n";
close OUT_MANIFEST;


# copy breakpoint file (delinearing conflicts between the sequence and genome maps), and align1 alignment files
$theSourceDir = "$opts{o}/assignAlignType";
copy "$theSourceDir/conflicts.txt", "$hybridScaffoldFinalResultsDir";	
$theSourceDir = "$opts{o}/align1"; 
my $prefixNGS_BNG = join("_", @enzymes)."_$modBNGfile"."_$NGSfile"."_BNGcontigs_NGScontigs";
copy "$theSourceDir/align1.xmap", "$hybridScaffoldFinalResultsDir/$prefixNGS_BNG.xmap";
editReferenceQueryMapsHeader("$hybridScaffoldFinalResultsDir/$prefixNGS_BNG.xmap", "$prefixNGS_BNG");
copy "$theSourceDir/align1_r.cmap", "$hybridScaffoldFinalResultsDir/$prefixNGS_BNG"."_r.cmap";
copy "$theSourceDir/align1_q.cmap", "$hybridScaffoldFinalResultsDir/$prefixNGS_BNG"."_q.cmap";

# copy the key files (pre- and post-cut key files)
my @filenameKeyContent = split(/\//, $filename_key);
my $filenameKeyFile = $filenameKeyContent[$#filenameKeyContent];	# take the actual file name of the key file
copy "$filename_key", "$hybridScaffoldFinalResultsDir";
# new key file would exist if conflict resolution for sequence was indicated by opts{N}
copy "$opts{o}/agp_fasta/$filenameKeyFile.cut.txt", "$hybridScaffoldFinalResultsDir" if (-e "$opts{o}/agp_fasta/$filenameKeyFile.cut.txt"); 

###

# align molecules BNX to BioNano and hyrbid scaffold [optional]
if (defined $opts{x})	{
	my $theRefAligner = $refaligner;
	my $globalMaxThreads = $maxthreads;
	my $deNovoPipelineDir = $opts{p};	# opts{p} should be absolute path
	my $deNovoOpArgsFile = $opts{q};	# opts{q} should be absolute path

	my ($moleculesDir, $moleculesFile) = separateDirFile(abs_path($opts{m}));	# opts{m} should be absolute path
	my ($bioNanoMapDir, $bioNanoMapFile) = separateDirFile(abs_path($mergeBN));	# the genome maps right before the merge step
	my ($hybridMapDir, $hybridMapFile) = ($hybridScaffoldFinalResultsDir, $hybridScaffoldFileName);

	my $autoNoiseOutDir = "$opts{o}/auto_noise";	$autoNoiseOutDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");
	my $alignMolBioNanoOutDir = "$opts{o}/alignmol_bionano";	$alignMolBioNanoOutDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");
	my $alignMolHybridOutDir  = "$opts{o}/alignmol_hybrid";	$alignMolHybridOutDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");

	printStatus(\*STATUS_FILE, "molecules alignment to BioNano maps and hybrid scaffolds", "Started", 0);
	@$cmdRef = ();
	$dtStartStage = DateTime->now;
	print "\nBeginning molecules alignment to BioNano genome maps and hybrid scaffolds...\n";
	unless (chdir $plDir)	{
		my $errMsg = "ERROR: Cannot change to $plDir: $!";
		printStatus(\*STATUS_FILE, "molecules alignment to BioNano maps and hybrid scaffolds", $errMsg, 1);
		dieLog("$errMsg\n");
	} # unless chdir
	@$cmdRef = ($^X, "scripts/align_molecules.pl", "-refAligner", $refaligner, "-globalMaxThreads", $globalMaxThreads, "-deNovoPipelineDir", $deNovoPipelineDir, "-deNovoOpArgsFile", $deNovoOpArgsFile, "-moleculesDir", $moleculesDir, "-moleculesFile", $moleculesFile, "-bioNanoMapDir", $bioNanoMapDir, "-bioNanoMapFile", $bioNanoMapFile, "-hybridMapDir", $hybridMapDir, "-hybridMapFile", $hybridMapFile, "-autoNoiseOutDir", $autoNoiseOutDir, "-alignMolBioNanoOutDir", $alignMolBioNanoOutDir, "-alignMolHybridOutDir", $alignMolHybridOutDir);
	
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	open(OUT, ">$autoNoiseOutDir/align_molecules.log");	print OUT $outResults;	close OUT;
	open(ERR, ">$autoNoiseOutDir/align_molecules.errlog");	print ERR $errResults;	close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	if ($outResults =~ /ERROR:/ || $errResults =~ /ERROR:/)	{
		my $errMsg = ($outResults =~ /ERROR:/) ? ($outResults) : ($errResults);
		printStatus(\*STATUS_FILE, "molecules alignment to BioNano maps and hybrid scaffolds", $errMsg, 1);
	} else	{
		printStatus(\*STATUS_FILE, "molecules alignment to BioNano maps and hybrid scaffolds", "Done", 1);	# close the status file
	} # if outResults
	errCheck($outResults, $errResults, "ERROR:", "Alignment of molecules to BioNano genome maps and hybrid scaffolds complete in $span.", "ERROR: Alignment of molecules to BioNano genome maps and hybrid scaffolds cannot be completed.");
} # if defined x 

#if (defined($opts{z}))	{
 if (!defined($opts{a})){
	### zipping the necessary results for IrysView
	@$cmdRef = ();
	$dtStartStage = DateTime->now;
	print "\nBeginning archive and zip of result files...\n";
	chdir $opts{o} or dieLog "ERROR: Cannot change directory to $opts{o}: $!\n";
	if (!defined($opts{z})){
		$opts{z} = "hybridScaffold_archive.tar.gz"
	}
	
	my $tarZipFile = "$opts{o}/$opts{z}";
	my $tarName = "$opts{o}/$opts{z}";	$tarName =~ s/\.tar\.gz$//;	# remove the extension
	# create the destination directory
	mkdir $tarName if (! -e $tarName);	# there should not be another tarName directory already (as it is unique)

	# if alignment of molecules are done
	if (defined($opts{x}))	{
		my $copySource = "alignmol_bionano";	$copySource .= "_M$modifyNum" if (defined($opts{M}));	$copySource .= "/merge";
		my $copyDest = "$tarName/alignmol_bionano";	$copyDest .= "_M$modifyNum" if (defined($opts{M}));	
		dircopy("$copySource/", "$copyDest/") or dieLog "ERROR: cannot copy directory $copySource to directory $copyDest\n";
		chdir $tarName or dieLog("ERROR: cannot change directory to $tarName: $!\n");
		my $tarZipSource = "alignmol_bionano";	$tarZipSource .= "_M$modifyNum" if (defined($opts{M}));	# relative path of copyDest
		my $tarZipDest = $tarZipSource;	$tarZipDest .= ".tar.gz";
		@$cmdRef = ("tar", "-czf", $tarZipDest, $tarZipSource);
		($outResults, $errResults) = runCommand($cmdRef);
		@$cmdRef = ();
		# error check
		errCheck($outResults, $errResults, "fail", "The tar zip command completed.", "ERROR: the tar zip command cannot be completed.");
		# remove copyDest directory 
		rmtree($copyDest) or dieLog("ERROR: cannot remove directory $copyDest for tar zip of $tarZipSource: $!\n");
		chdir $opts{o} or dieLog("ERROR: Cannot change directory to $opts{o}: $!\n");
		
		$copySource = "alignmol_hybrid";	$copySource .= "_M$modifyNum" if (defined($opts{M}));	$copySource .= "/merge";
		$copyDest = "$tarName/alignmol_hybrid";	$copyDest .= "_M$modifyNum" if (defined($opts{M}));	
		dircopy("$copySource/", "$copyDest/") or dieLog("ERROR: cannot copy directory $copySource to directory $copyDest\n");
		chdir $tarName or dieLog("ERROR: cannot change directory to $tarName: $!\n");
		$tarZipSource = "alignmol_hybrid";	$tarZipSource .= "_M$modifyNum" if (defined($opts{M}));	# relative path of copyDest
		$tarZipDest = $tarZipSource;	$tarZipDest .= ".tar.gz";
		@$cmdRef = ("tar", "-czf", $tarZipDest, $tarZipSource);
		($outResults, $errResults) = runCommand($cmdRef);
		@$cmdRef = ();
		# error check
		errCheck($outResults, $errResults, "fail", "The tar zip command completed.", "ERROR: the tar zip command cannot be completed.");
		# remove copyDest directory
		rmtree($copyDest) or dieLog("ERROR: cannot remove directory $copyDest for tar zip of $tarZipSource: $!\n");
		chdir $opts{o} or dieLog("ERROR: Cannot change directory to $opts{o}: $!\n");
	} # if defined molecules

	# make another dummy fa2cmap directory
	mkdir "$tarName/fa2cmap" if (! -e "$tarName/fa2cmap");	
	# copy necessary hybrid scaffolds _M directory
	my $copySource = "hybrid_scaffolds";	$copySource .= "_M$modifyNum" if (defined($opts{M}));
	dircopy("$copySource/", "$tarName/$copySource/") or dieLog("ERROR: cannot copy directory $copySource to directory $tarName\n");
	# copy cur_results.txt
	$copySource .= "/cur_results.txt";
	copy($copySource, $tarName) or dieLog("ERROR: cannot copy file $copySource to directory $tarName\n");
	# the status file 
	$copySource = "$opts{o}/$statusFile";
	copy($copySource, $tarName) or dieLog("ERROR: cannot copy file $copySource to directory $tarName\n");
		
	
	# the actual tar zip command
	my $tarNameRelativePath = $opts{z};	$tarNameRelativePath =~ s/\.tar\.gz$//;
	print "tar file: $tarZipFile\n";
	@$cmdRef = ( "tar", "-czf", $tarZipFile, $tarNameRelativePath);	# tar zip the directory tarName into tarZipFile			 
	($outResults, $errResults) = runCommand($cmdRef);
	@$cmdRef = ();
	# error check
	errCheck($outResults, $errResults, "fail", "The tar zip command completed.", "ERROR: the tar zip command cannot be completed.");
	# remove the destination directory as the zip file is already complete
	rmtree($tarName) or dieLog("ERROR: cannot remove directory $tarName used for tar zip: $!\n");
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheck($outResults, $errResults, "fail", "Archive and zip of result files completed in $span.", "ERROR: archive and zip of result files cannot be completed.");
} # if not defined opts{a}


print "\n";
my $dtEnd = DateTime->now;
print "End time: "; print join ' ', $dtEnd->ymd, $dtEnd->hms; print "\n\n";

$span = DateTime::Format::Human::Duration->new();
print 'Total elapsed time: ', $span->format_duration_between($dtEnd, $dtStart); print "\n";

print "\n\nMerging of $opts{b} with $opts{n} is complete.\n\n";


print "END of output\n";





######################################################################
#                           Subroutines                              #
######################################################################

# this subroutine just prints message to the file handle
# if closeFlag is 1, then also close the file handle
sub printStatus	{
	my ($fileHandle, $stageName, $stageStatus, $closeFlag) = @_;
	print $fileHandle "{Stage:\"$stageName\",Status:\"$stageStatus\"}\n";
	if ($closeFlag == 1)	{
		close $fileHandle;
	} # if closeFlag
} # printStatus

# this subroutine is used to check RefAligner processes to see if it finishes successfully
# it assumes that -stderr and -stdout was used to run RefAligner
sub errCheckRefAligner	{
	my ($file, $completeString, $completeMsg, $dieMsg, $exitFlag) = @_;
	open(IN, "$file") or dieLog ("ERROR: Cannot open $file: $!\n");
	my $presence = 0;
	while (my $line = <IN>)	{
		if ($line =~ /$completeString/)	{
			# if the line contains the string that indicates successful completion
			$presence = 1;
		} # if line
	} # while line
	close IN;

	return $presence if ($exitFlag != 1);
	
	if ($presence == 0)	{
		dieLog ("ERROR: $dieMsg\n");
	} else	{
		print "$completeMsg\n";
	} # if presence
	return $presence;
} # errCheckRefAligner

# this subrountine is used to call non-RefAligner processes to see if there is error in the child
# it assumes that the child process prints out "ERROR:" tag to all messages
sub errCheck	{
	my ($outResults, $errResults, $errorString, $completeMsg, $dieMsg) = @_;
	if ($outResults =~ /$errorString/)	{
		dieLog ("$outResults\nERROR: $dieMsg\n"); } 
	elsif ($errResults =~ /$errorString/)	{
		dieLog ("$errResults\nERROR: $dieMsg\n"); } 
	else {
		print "$completeMsg\n"; } 
}

sub runCommand  {
        my ($argsRef) = @_;
#print "runCommand: before open3\n";
        my $pid = open3(my $CMD_IN, my $out, my $err, @$argsRef);
#print "runCommand: after open3\n";

        close($CMD_IN);

        my $outResults = "";
        my $errResults = "";
        my $sel = new IO::Select;
        $sel->add($out, $err);
        while(my @fhs = $sel->can_read) {
                foreach my $fh (@fhs) {
                        my $line = <$fh>;
                        unless(defined $line) {
                                $sel->remove($fh);
                                next;
                        } # unless line
                        if($fh == $out) {
                                $outResults .= "$line";
                                #print "$line";
                        }elsif($fh == $err) {
                                $errResults .= "$line";
                                #print "$line";
                        }else{
                                dieLog ("ERROR: This should never execute!");
                        } # if fh
                } # foreach fh
        } # while
#print "runCommand: before waitpid\n";
        my $ret=waitpid ($pid, 0); # reap the exit code
#print "runCommand: after waitpid\n";
        return ($outResults, "$errResults");
} # runCommand

sub Init{
	my $opt_string = 'hn:b:c:o:r:B:N:fxym:vM:p:q:e:z:Sw:tu:a';
	if(!getopts("$opt_string", \%opts)){
		print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	Usage() if $opts{h};
	if ($opts{v}) {
		Version();
		exit 0; 
	} # if opts v
	
	if (! defined $opts{b} || $opts{b} !~ /\w+/ || ! -e $opts{b} )	{
		 print ("\nERROR: Please provide a valid BioNano genome map file, option -b\n");
		Usage();
	}
	if (! defined $opts{n} || $opts{n} !~ /\w+/ || ! -e $opts{n})	{
		print ("\nERROR: Please provide a valid sequence fasta file, option -n\n");
		Usage();
	} # if 
	if (! defined $opts{c} || $opts{c} !~ /\w+/ || ! -e $opts{c})	{
		print ("\nERROR: Please provide a valid XML configuration file for the hybrid scaffold pipeline, option -c\n");
		Usage();
	} # if 
	if (! defined $opts{r} || $opts{r} !~ /\w+/)	{
		print ("\nERROR: Please enter the RefAligner binary file, option -r\n");
		Usage();
	} # if 
	
	$opts{b} = abs_path($opts{b});
	$opts{n} = abs_path($opts{n}); 
	$opts{c} = abs_path($opts{c});
	$opts{r} = abs_path($opts{r});
	
	
	# output directory	
	if(defined($opts{o})) {
		$opts{o} = abs_path($opts{o}); 
	} # if 
	
	# user defined conflict resolution file, which can affect output directory, B and N options, and f overwrite option
	if (! defined($opts{M}))	{
		# if the user annotated conflict resolution file is not provided, then, must specify how to resolve conflict
		if ((! defined($opts{N}) || ! defined($opts{B})) || ($opts{N} !~ /^\d$/ || $opts{B} !~ /^\d$/ || !(1 <= $opts{B} && $opts{B} <= 3) || !(1 <= $opts{N} && $opts{N} <= 3)))	{
			print ("\nERROR: Please input in a conflict-resolution value between 1 and 3 for -N and -B to handle conflicts.\n");
			Usage();
		} # if optsB

		# a conflict status file was not given, then start anew
		if ($opts{f})	{
			rmtree($opts{o}) if (-d $opts{o});	# if exists an output directory
		} else	{
			dieLog( "\nERROR: Output folder: $opts{o} already exists! Use -f option to overwrite.\n") if (-d $opts{o});
		} # if force overwrite
		eval { mkpath($opts{o}) };
		if ($@) {
			print("\nERROR: Couldn't create $opts{o}: $@");
			print("\nERROR: Output folder invalid or does not exist! \n");
                	Usage();
        	} # if
	} else	{
		# check if the user annotated conflict resolution file exists
		if (! -e $opts{M})	{
			dieLog ("\nERROR: the conflict status file does not exist: $opts{M}\n");
		} # if not exists
		
		# a conflict status file was given, then there must already be data from a previous hybrid scaffold run (ignore -f)
		if (! -d $opts{o})	{
			dieLog("\nERROR: Output folder: $opts{o} does not exist! Please point to the output of a previous hybrid scaffold run.\n");
		} # if exists an output direcotry
                $opts{M} = abs_path($opts{M});
	} # defined M

	# user wants to perform alignment of molecules to BioNano genome maps and Hybrid scaffolds
	if ($opts{x})	{
		if (! defined($opts{m}) || $opts{m} !~ /\w+/ || ! -e $opts{m})	{
			print ("\nERROR:  Please provide a valid molecules BNX file for alignment to BioNano genome maps and hybrid scaffolds, option -m\n");
			Usage();
		} # check molecule file
		if (! defined($opts{p}) || $opts{p} !~ /\w+/ || ! -d $opts{p})	{
			print ("\nERROR: alignment of molecules option has been indicated, but a valid de novo assembly pipeline script directory is not provided, option -p\n");
			Usage();
		} # check pipeline script directory
		if (! defined($opts{q}) || $opts{q} !~ /\w+/ || ! -e $opts{q})	{
			print ("\nERROR: alignment of molecules option has been indicated, but a valid de novo assembly pipeline parameter file (optArguments.xml) is not provided, option -q\n");
			Usage();
		} # check pipeline xml file
		$opts{m} = abs_path($opts{m});
		$opts{p} = abs_path($opts{p});
		$opts{q} = abs_path($opts{q});
	} # if opts{x}

	# user wants to perform chimeric quality score generation for the input BioNano genome maps
	if ($opts{y})	{
		if (! defined($opts{m}) || $opts{m} !~ /\w+/ || ! -e $opts{m})	{
			print ("\nERROR:  Please provide a valid molecules BNX file for the generation of chimeric quality score for the input BioNano genome maps, option -m\n");
			Usage();
		} # check molecule file
		if (defined($opts{e}) && ! -e $opts{e} )	{
			print ("\nERROR: Generation of chimeric quality score option has been indicated, and a noise parameter file was provided, but it is not valid, option -e\n");
			Usage();
		} # check the .err file or .errbin file
		$opts{m} = abs_path($opts{m});
		$opts{e} = abs_path($opts{e}) if (defined($opts{e}));	
	} # if opts{y}

	# for IrysView to archive file
	if (defined($opts{z}) && $opts{z} !~ /tar\.gz$/)	{
		print ("\nERROR: Please provide a valid name for a tar gzip file (e.g. hybrid_results.tar.gz), option -z\n");
		Usage();
	} # if opts{z}
} # Init

sub Version{
	my $dtStartStage = DateTime->now;

	my @revs;
	###my @pipelineFiles = process_files($plDir);
	# only process a shorter list of files to find the version number
	my @pipelineFiles = process_files_shortList($plDir);
	foreach (@pipelineFiles) {
		if (!-d $_) {
			#if ($_ !~ /svn/ && $_ !~ /perl5/) {
				open (FILE,$_);
				my @file = <FILE>; 
				close FILE;
				my $start = '$Id:'; my $end = '$';
				#print "$_\n";
				foreach (@file) {
					if ($_ =~ m/\$Id/) { 
						chomp($_);
						my ($wanted) = $_ =~ /\$Id:(.*)\$/;
						$wanted =~ s/^\s+|\s+$//g;
						push @revs, $wanted;
						#print $wanted."\n"; 
						last; 
						} } } } 
						#}
	
	my $revLine = ""; my $revNum=0;
	@revs = uniq @revs;
	foreach my $rev (@revs)	{
		my @content = split(/\s+/, $rev);
		next if (scalar(@content) < 2);
		if ($content[1] >= $revNum)	{
			$revNum = $content[1];
			$revLine = $rev;
		} # if content
	} # foreach rev

	
	print "$revLine\n";  

	my $dtEndStage = DateTime->now;
	my $span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	#print "Total time: $span\n";
} # Version

sub process_files_shortList	{
	my ($thePath) = @_;
	my @theFiles = ();
	opendir(DIR, "$thePath") or dieLog("process_files_shortList: cannot open dir $thePath: $!\n");
	my @topFiles = grep {$_ =~ /\.pl$/ || $_ =~ /\.pm$/} readdir DIR;	# add hybridScaffold.pl
	closedir DIR;
	foreach my $topFile (@topFiles)	{
		push(@theFiles, "$thePath/$topFile");
	} # foreach topFile

	# scripts sub-directory
	opendir(DIR, "$thePath/scripts") or dieLog("process_files_shortList: cannot open dir $thePath/scripts: $!\n");
	my @secLevelFiles = grep {$_ =~ /\.pl$/ || $_ =~ /\.pm$/} readdir DIR;
	closedir DIR;	
	foreach my $secLevelFile (@secLevelFiles)	{
		push(@theFiles, "$thePath/scripts/$secLevelFile");
	} # foreach secLevelFile

	# BNG sub directory
	opendir(DIR, "$thePath/scripts/perl5/BNG") or dieLog("process_files_shortList: cannot open dir $thePath/scripts/perl5/BNG: $!\n");
	my @thirdLevelFiles = grep {$_ =~ /\.pl$/ || $_ =~ /\.pm$/} readdir DIR;
	closedir DIR;
	foreach my $thirdLevelFile (@thirdLevelFiles)	{
		push(@theFiles, "$thePath/scripts/perl5/BNG/$thirdLevelFile");
	} # foreach thirdLevelFile

	return @theFiles;
} # process_files_shortList

sub process_files {
	my $path = shift;
    opendir (DIR, $path)
        or dieLog ("ERROR: Unable to open $path: $!");

    # We are just chaining the grep and map
    # This is the same as: LIST = map(EXP, grep(EXP, readdir()))
    my @files =
        # Third: Prepend the full path
        map { $path . '/' . $_ }
        # Second: take out '.' and '..'
        grep { !/^\.{1,2}$/ }
        # First: get all files
        readdir (DIR);

    closedir (DIR);

    for (@files) {
        if (-d $_) {
            # Add all of the new files from this directory
            # (and its subdirectories, and so on... if any)
            if (abs_path($_) !~ /5.10.1/ && abs_path($_) !~ /5.14.4/ && abs_path($_) !~ /svn/ && abs_path($_) =~ /scripts/) {
				push @files, process_files ($_); } }

        else {
            # Do whatever you want here =) .. if anything.
        }
    }
    # NOTE: we're returning the list of files
    return @files;
}

sub Usage{
	print << "EOF";
	
Usage: perl hybridScaffold.pl <-h> <-n ngs_file> <-b bng_cmap_file> <-c hybrid_config_xml> <-o output_folder> <-B conflict_filter_level> <-N conflict_filter_level> <-f> 
      <-m molecules_bnx> <-p de_novo_pipeline> <-q de_novo_xml> <-v> <-x> <-y> <-e noise_param><-z tar_zip_file><-S>
      -h    : This help message         
      -n    : Input NGS FASTA [required]
      -b    : Input BioNano CMAP  [required]
      -c    : Merge configuration file [required]
      -o    : Output folder [required]
      -r    : RefAligner program [required]
      -B    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contig [required if not using -M option]
      -N    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contig [required if not using -M option]
      -f    : Force output and overwrite any existing files
      -x    : Flag to generate molecules to hybrid scaffold alignment and molecules to genome map alignment [optional]
      -y    : Flag to generate chimeric quality score for the Input BioNano CMAP [optional]
      -m    : Input BioNano molecules BNX [optional; only required for either the -x or -y option]
      -p    : Input de novo assembly pipeline directory [optional; only required for -x option]
      -q    : Input de novo assembly pipeline optArguments XML script [optional; only required for -x option]
      -e    : Input de novo assembly noise parameter .errbin or .err file [optional; recommended for -y option but not required]
      -v    : Print pipeline version information
      -M    : Input a conflict resolution file indicating which NGS and BioNano conflicting contigs to be cut [optional] 
      -z    : Name of a zipped file to archive the essential output files [optional]
      -S    : Only run hybridScaffold up to before Merge steps [optional]
      -w    : Name of the status text file needed by IrysView [optional]
      -t    : Perform pre-pairmerge sequence to pre-pairmerge genome map alignment [optional]
      -u    : Sequence of enzyme recognition site (overrides what has been specified in config XML file, for IrysView only) [optional]
      
EOF
	exit 1;
} # Usage

sub numeric	{	$a	<=>	$b	}

sub parseConfig	{
	my ($configRef, $theStage, $commandStackRef) = @_;
	if(!exists $configRef->{$theStage}){
	    return $commandStackRef;
	}
	foreach my $flagIncludeKey (keys %{$configRef->{$theStage}})	{
		if ($flagIncludeKey =~ /flag/)	{
			# inside flag section
			# check to see if there is just one flag, in which case under flag is a hash table; if more than one entry under flag, then an array
			if (ref($configRef->{$theStage}{$flagIncludeKey}) eq "ARRAY")	{
			# iterate the array of attributes
			for (my $i = 0; $i < scalar(@{$configRef->{$theStage}{$flagIncludeKey}}); $i += 1)	{
				my $tempAttr = "";
				my %tempVals = ();
				foreach my $flagKey (keys %{$configRef->{$theStage}{$flagIncludeKey}[$i]})	{
					# capture the attribute and values
					if ($flagKey =~ /attr/)	{
						# check if attr has "-" in front
						$tempAttr = ($configRef->{$theStage}{$flagIncludeKey}[$i]{$flagKey} !~ /^-/) ? ("-$configRef->{$theStage}{$flagIncludeKey}[$i]{$flagKey}") : ($configRef->{$theStage}{$flagIncludeKey}[$i]{$flagKey});
					} elsif ($flagKey =~ /val(\d+)/)	{
						$tempVals{$1} = $configRef->{$theStage}{$flagIncludeKey}[$i]{$flagKey};
					} else	{
						# no op, only want attr and val\d+
					} # if flagKeys
				} # foreach flagKeys
				# at the end of the XML line, now assign the sorted values to attr
				if ($tempAttr ne "")	{
					# make sure there is an attribute to store
					push(@$commandStackRef, $tempAttr);
					foreach my $vals(sort numeric keys %tempVals)	{
						push(@$commandStackRef, $tempVals{$vals});
					} # foreach vals
				} # if tempAttr
			} # for i
			} elsif (ref($configRef->{$theStage}{$flagIncludeKey}) eq "HASH")	{
				my $tempAttr = "";
				my %tempVals = ();
				foreach my $flagKey (keys %{$configRef->{$theStage}{$flagIncludeKey}})	{
					# capture the attribute and values
					if ($flagKey =~ /attr/)	{
						# check if attr has "-" in front
						$tempAttr = ($configRef->{$theStage}{$flagIncludeKey}{$flagKey} !~ /^-/) ? ("-$configRef->{$theStage}{$flagIncludeKey}{$flagKey}") : ($configRef->{$theStage}{$flagIncludeKey}{$flagKey});
					} elsif ($flagKey =~ /val(\d+)/)	{
						$tempVals{$1} = $configRef->{$theStage}{$flagIncludeKey}{$flagKey};
					} else	{
						# no op, only want attr and val\d+
					} # if flagKey
				} # foreach flagKey
				if ($tempAttr ne "")	{
					# make sure there is an attribute to store
					push(@$commandStackRef, $tempAttr);
					foreach my $vals (sort numeric keys %tempVals)	{
						push(@$commandStackRef, $tempVals{$vals});
					} # foreach vals
				} # if tempAttr
			} else	{
				# op operation, as it must be an array or hash under flag
			} # if under flag is of type array or hash
		} # if flag
		if ($flagIncludeKey =~ /include/)	{
			# inside include section
			# check if there are multiple includes, in which case it will be an array under include; single include, then a hash under include	
			if (ref ($configRef->{$theStage}{$flagIncludeKey}) eq "ARRAY")	{
				# iterate each include
				for (my $i = 0; $i < scalar(@{$configRef->{$theStage}{$flagIncludeKey}}); $i += 1)	{
					foreach my $includeKey (keys %{$configRef->{$theStage}{$flagIncludeKey}[$i]})	{
						if ($includeKey =~ /val(\d+)/)	{
							# find out which stage to include, and recursively parse that stage
							my $childStage = $configRef->{$theStage}{$flagIncludeKey}[$i]{$includeKey};
							$commandStackRef = parseConfig($configRef, $childStage, $commandStackRef);
						} # if includeKey
					} # foreach includeKey
				} # for i
			} elsif (ref ($configRef->{$theStage}{$flagIncludeKey}) eq "HASH")	{
				# only 1 stage to include
				foreach my $includeKey (keys %{$configRef->{$theStage}{$flagIncludeKey}})	{
					if ($includeKey =~ /val(\d+)/)	{
						# find out which stage to include, and recursively parse that stage
						my $childStage = $configRef->{$theStage}{$flagIncludeKey}{$includeKey};
						$commandStackRef = parseConfig($configRef, $childStage, $commandStackRef);
					} # if includeKey
				} # foreach includeKey
			} else	{
				# no op, as it must an array or hash under include
			} # if under include is of type array or hash type
		} # if include 
	} # foreach flagIncludKeys
	return $commandStackRef
} # parseConfig

sub extractConfig	{
	my ($commandStackRef, $theKey) = @_;
	# this subrountine iterates the commandStack array, and extracts the values of a given key
	# it returns an array containing all the values of that key
	my @theValues = ();
	# remember that a key must have "-" in front, and all subsequent entries  without "-" are the values
	$theKey = "-$theKey";
	my $curKey = "";	
	my $found = -1;
	for (my $i = 0; $i < scalar(@$commandStackRef) && $found != 1; $i += 1)	{
		my $curEntry = $commandStackRef->[$i];
		if ($curEntry =~ /^-/)	{
			# entry is a key
			if ($curEntry =~ /^$theKey$/)	{
				# this is the key of interest
				$found = 0;	# flip found to 0
			} # if seeing the key of interest 
			$curKey = $curEntry;

			# now check if we already have seen the key of interest
			if ($curKey !~ /^$theKey$/ && $found == 0)	{
				$found = 1;
			} # already seen
		} else	{
			# entry is a value
			if ($curKey =~ /^$theKey$/ && $found == 0)	{
				push(@theValues, $curEntry);
			} # if curKey
		} # if curEntr
	} # for i
	return \@theValues;
} # extractConfig

sub copyArray 	{
	my ($toArrayRef, $fromArrayRef) = @_;
	for (my $i = 0; $i < scalar(@$fromArrayRef); $i += 1)	{
		push(@$toArrayRef, $fromArrayRef->[$i]);
	} # for i
	return $toArrayRef;
} # copyArray

sub find_refaligner {
	if (! -d $_) {
		if ($_ =~ /RefAligner/i) {
			$refaligner = $File::Find::name; }}
}

sub log10 {
	my $n = shift;
    return log($n)/log(10);
}

# subroutine to determine the modify number for this current run (by looking at the largest modify number, if exists, in the "hybrid_scaffolds" output directory)
sub getCurMaxModifyNum	{
	my ($aDir) =@_;
	my $modifyNum = 1;
	opendir(DIR, $aDir) or dieLog ("ERROR: cannot access output directory $aDir: $!\n");
	my @alignFinalDirs = grep {$_ =~ /^hybrid_scaffolds/i} readdir DIR;
	closedir DIR;

	my @alignFinalModifyNum = ();
	for (my $i = 0; $i < scalar(@alignFinalDirs); $i += 1)	{
		if ($alignFinalDirs[$i] =~ /.+_M(\d+)$/)	{
			push(@alignFinalModifyNum, $1);
		} # if alignFinalDirs
	} # for i

	if (scalar(@alignFinalModifyNum) > 0)	{
		@alignFinalModifyNum = sort{$a <=> $b}(@alignFinalModifyNum);
		$modifyNum = $alignFinalModifyNum[$#alignFinalModifyNum] + 1;
	} # if scalar

	return $modifyNum;
} # getCurMaxModifyNum

# this subroutine determines if chimeric quality scores exist in a given file
sub getQScores	{
	my ($file, $qScoresRef) = @_;
	# read in the genome map cmap file, ASSUMING that the chimQuality is on column 10
	my $noQScoreFlag = 0;
	open(IN, "$file") or dieLog "getQScores: cannot open $file: $!\n";
	my $foundChimQuality = -1;
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# assumes that data line comes after header lines
			if ($line =~ /^#h\s+/)	{
				# find that column storing the ChimQuality
				$line =~ s/^#h\s+//;
				my @headerContent = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@headerContent); $i += 1)	{
					$foundChimQuality = $i if ($headerContent[$i] =~ /ChimQuality/i);
				} # for i
			} # if line
			next;
		} # if line
		$line =~ s/^\s+/\t/;    $line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		if ($foundChimQuality == -1 || $content[$foundChimQuality] !~ /\d+/)	{
			warn "WARNING: getQScores cmap file=$file does not have chimeric quality score\n";
			$noQScoreFlag = 1;
			%$qScoresRef = ();	# empty any qScores currently stored, and just flag this file as having no qScores
			return ($qScoresRef, $noQScoreFlag);
		} # if scalar
		my ($id, $labelChannel, $start, $coverage, $chimQuality) = ($content[0], $content[4], $content[5], $content[7], $content[$foundChimQuality]);
		# skip if label channel is 0    (line indicating contig length)
		next if ($labelChannel == 0);
		push(@{$qScoresRef->{$id}}, {start => $start, coverage => $coverage, chimQuality => $chimQuality});
	} # while line
	close IN;
	return ($qScoresRef, $noQScoreFlag);
} # getQScores

sub editReferenceQueryMapsHeader	{
	my ($file, $headerSubString) = @_;
	# the sole purpose of this subroutine is to edit the # Reference Maps From: and # Query Maps From: header lines (so that IrysView will not complaining)
	my @lines = ();
	open(IN, $file) or dieLog("ERROR: editReferenceQueryMapsHeader: cannot read in $file: $!\n");
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^(#\s+Reference\s+Maps\s+From:)/)	{
			# _r.cmap
			$line = "$1\t$headerSubString"."_r.cmap";
		} elsif ($line =~ /^(#\s+Query\s+Maps\s+From:)/)	{
			# _q.cmap
			$line = "$1\t$headerSubString"."_q.cmap";
		} else	{
			# no op
		} # if line
		push(@lines, $line);
	} # while line
	close IN;

	open(OUT, ">$file") or dielog("ERROR: editReferenceQueryMapsHeader: cannot write to $file: $!\n");
	foreach my $line (@lines)	{
		print OUT "$line\n";
	} # foreach line
	close OUT;
} # editReferenceQueryMapsHeader

# this subrountine reads in a cmap file, and determines its largest cmap id
sub readCmapIdLength	{
	my ($file) = @_;
	# read in a cmap file, and determine the largest sequence id
	my $maxSeqId = -1;
	my %theLengths = ();
	open(IN, "$file") or dieLog "ERROR: readCmapId: reading in $file: $!\n";
	while (my $line = <IN>) {
		chomp $line;    $line =~ s/\r//g;
		next if ($line =~ /^#/);        # skip header lines
		my @content = split(/\t/, $line);
		my ($id, $aLength) = @content[0..1];
		$aLength = int($aLength);       # again change to integer
		$theLengths{$id} = $aLength;
		$maxSeqId = ($maxSeqId < $id) ? ($id) : ($maxSeqId);
	} # while line
	close IN;
	return ($maxSeqId, \%theLengths);
} # readCmapIdLength

sub separateDirFile     {
        my ($fullFile) = @_;
        my ($dir, $file) = ("" ,"");
        if ($fullFile =~ /^$/)  {
                return ($dir, $file);
        } # if fullFile
        my $tempFile = $fullFile;
        my @tempContent = split(/\//, $tempFile);
        $file = $tempContent[$#tempContent];

        if ($fullFile =~ /^\// && scalar(@tempContent) == 2)    {
                # the file is in the "/" the very top directory
                # tempContent[0] is "", tempContent[1] is "$file"
                return ("/", $file);
        } # if fullFile
        $dir = $tempFile;
        $dir =~ s/$file$//;
        $dir =~ s/\/$//;        # remove trailing /
        # if only a file with no path has been specified (meaning the file is in the current directory)
        # fullFile did not start with "/" either
        $dir = ($dir =~ /^$/) ? (".") : ($dir); 
        return ($dir, $file);
} # separateDirFile

