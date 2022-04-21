#20170728 EL - run molecule-to-reference alignment
#Run SNR filtering for all label channels, auto-noise, and alignmolvref
#The script assumes that you have access to RefAligner and python assembly pipeline

description = """Wrapper for running molecule-to-reference alignment using Pipeline scripts.
Final alignment results are in <output>/contigs/alignmolvref/merge/"""

import argparse
import os
import subprocess
import sys
import tarfile
#import shutil

#Run SNR filtering on all label channels (histogram-based method).
#This script should do SNR filtering by histogram method for BOTH color channels (if present)
# and if that fails, should apply default SNR threshold: currently uses only one snrThresh argument
# for both channels (first command line arg to this script), but should be extended to two in the future
def run_snr_filtering(sample_prefix, sample_bnx, pipeline, ra, snrThresh):

    out_bnx = sample_prefix + "_snr_filtered.bnx" #create output filename

    path_export_command = 'export PATH=$PATH:' + ra + '; ' #export path

    args = ['/usr/bin/perl']
    args += [pipeline + '/filter_SNR_dynamic.pl']
    args += ['-i', sample_bnx]
    args += ['-o', out_bnx]
    args += ['-d', str(snrThresh)]

    args = ' '.join(args)
    args = path_export_command + args

    print args

    os.system(args)

    return(out_bnx)

#call RefAligner with -minSNR X Y (for each color)
# here snrThresh is a vector of size 1 OR 2
def run_snr_filteringRA(sample_prefix, sample_bnx, pipeline, ra, snrThresh):

    out_bnx = sample_prefix + "_snr_filtered" #create output filename: prefix only for RefAligner

    args = [os.path.join(ra,"RefAligner"), "-o", out_bnx, "-stdout", "-merge", "-bnx", "-f",
            "-i", sample_bnx]
    args += ["-minSNR", str(snrThresh[0])]
    if len(snrThresh) > 1 : #should always be true, but RefAligner is ok either way
        args += [str(snrThresh[1])]

    print "Creating SNR filtered bnx", out_bnx+".bnx"

    stdout = joberr = ""
    retc = 1
    try:
        retc, stdout, joberr = util.runJob(args, returnstdout=True) #bc returnstdout is True, return tuple
    except OSError, e: 
        print "Error: SNR filter RefAligner job failed (%s), check output: %s" % (e,out_bnx+".stdout")
        return None

    return out_bnx+".bnx"

#Run autonoise to estimate error parameters
def run_autonoise(outdir, sample_prefix, sample_bnx, ref, pipeline, ra, nthreads, optArgs, color):

    args = ['python2.7']
    args += [pipeline + '/pipelineCL.py']
    args += ['-y', '-x', '-z']
    args += ['-T', nthreads]
    args += ['-j', nthreads]
    args += ['-b', sample_bnx]
    args += ['-l', outdir]
    args += ['-e', sample_prefix]
    args += ['-a', optArgs]
    args += ['-t', ra]
    args += ['-r', ref]
    args += ['-F', str(color)]

    print ' '.join(args)

    s1 = subprocess.Popen(args)
    s1.wait()

    out_bnx =  os.path.join(outdir, 'contigs/auto_noise/autoNoise1_rescaled.bnx')

    return(out_bnx)

#Run molecule-to-reference alignment based on rescaled molecules
def run_alignmol(outdir, sample_bnx, ref, pipeline, ra, nthreads, optArgs, color):

    errbin = os.path.join(outdir, 'contigs/auto_noise/autoNoise1.errbin')
    if not util.checkFile(errbin) :
        print "Error: .errbin file not found, check autonoise job: " + errbin
        sys.exit(1)

    args = ['python2.7']
    args += [pipeline + '/runAlignMol.py']
    args += ['-q', ref]
    args += ['-r']
    args += ['-b', sample_bnx]
    args += ['-a', optArgs ]
    args += ['-o', os.path.join(outdir, 'contigs/alignmolvref')]
    args += ['-t', ra]
    args += ['-T', nthreads]
    args += ['-j', nthreads]
    args += ['-E', errbin]
    args += ['-p', pipeline]
    args += ['-F', str(color)]

    print ' '.join(args)

    s1 = subprocess.Popen(args)
    s1.wait()

    #return(" ".join(args)) #not used

#Tar and compress alignment output
def compress_output(outdir):

    filelist = []
    filepath = outdir + "/contigs/alignmolvref/merge/"
    
    output_tar = outdir + "/alignments.tar.gz"

    #traverse output directory and identify files of interest
    for root, dirs, files in os.walk(filepath):
        for file in files:
            if ".xmap" in file or ".cmap" in file:
                filelist.append(root + file)

    #create tar ball
    tar = tarfile.open(output_tar, "w:gz")
    for file in filelist:
        tar.add(file, arcname=file.split("/")[-1])
    tar.close()

#Get input parameters
def get_args():

    parser = argparse.ArgumentParser(description=description)

    snrDefault = 3.5 #default for both channels
    optArgsDefault = "optArguments_haplotype_saphyr.xml"

    parser.add_argument('--prefix', dest='prefix', help='sample name prefix (required)', type=str, required=True)
    parser.add_argument('--mol', dest='mol', help='input molecule bnx (required)', type=str, required=True)
    parser.add_argument('--ref', dest='ref', help='input reference cmap (required)', type=str, required=True)
    parser.add_argument('--ra', dest='ra', help='RefAligner directory (required)', type=str, required=True)
    parser.add_argument('--nthreads', dest='nthreads', help='Number of threads (required)', type=str, required=True)
    parser.add_argument('--output', dest='output', help='output dir (optional, default pwd)', type=str, default="")
    optArgsHelp = 'optArguments.xml file (optional, default "%s" in RefAligner dir)' % optArgsDefault
    parser.add_argument('--optArgs', dest='optArgs', help=optArgsHelp, type=str, default=optArgsDefault)
    parser.add_argument('--snrFilter', help='Label SNR filter method: 0 for histogram (default), 1 for fixed threshold (see below), 2 for NO filter', type=int, default=0)
    #nargs="+" means that at least one arg must be given (but this is not required, so can also omit the option)
    snrThreshHelp = 'Label SNR threshold: if 1 specified for previous argument (or 0 specified, but insufficient data), use argument as fixed threshold (ignored if histogram method used or no SNR filter); if 2-color bnx, can specify threshold for second color as optional second argument [default %.1f, %.1f]' % (snrDefault, snrDefault)
    parser.add_argument('--snrThreshold', help=snrThreshHelp, type=float, default=[snrDefault, snrDefault], nargs="+")
    parser.add_argument('--dosnrThresh', help='determines whether to create threshold labels based on SNR', action='store_true')
    parser.add_argument('--color', help='Color channel for alignment: replace -usecolor X in optArgs with this, must be either 1 or 2 [default OFF]', default=0, type=int)
    parser.add_argument('--pipeline', dest='pipeline', help='Pipeline directory (optional, defaults to script dir)', type=str)

    args = parser.parse_args()

    if args.pipeline: #check all Pipeline dependencies
        cwd = args.pipeline
    else:
        cwd = os.path.split(os.path.realpath(__file__))[0] #this is path of this script
        if not os.path.isfile(os.path.join(cwd,"utilities.py")): #if still not here, last try is actual cwd
            cwd = os.getcwd() #still check this below
    args.pipeline = cwd

    #this is the only one imported here and in runCharacterize
    if not os.path.isfile(os.path.join(cwd,"utilities.py")):
        print "ERROR: utilities.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)

    sys.path.append(cwd)
    global util
    import utilities as util

    if not util.checkFile(args.mol, ".bnx"):
        print "Error: bnx file not found or doesn't end with .bnx: " + args.mol
        sys.exit(1)

    if not util.checkFile(args.ref, ".cmap"):
        print "Error: reference file not found or doesn't end with .cmap: " + args.ref
        sys.exit(1)

    if not util.checkDir(args.ra, checkWritable=False, makeIfNotExist=False):
        print "Error: RefAligner dir not found: " + args.ra
        sys.exit(1)

    if not util.checkFile(os.path.join(args.ra, "RefAligner")):
        print "Error: RefAligner binary not found in dir: " + args.ra
        sys.exit(1)

    if args.optArgs == optArgsDefault:
        args.optArgs = os.path.join(args.ra, optArgsDefault)
    if not util.checkFile(args.optArgs, ".xml"):
        print "Error: optArgs file not found or doesn't end with .xml: " + args.optArgs
        sys.exit(1)

    if args.output:
        if not util.checkDir(args.output): #this will make args.output if doesn't exist
            print "Error: output dir cannot be created: " + args.output
            sys.exit(1)
    else :
        args.output = os.getcwd()

    if args.snrFilter < 0 or args.snrFilter > 2:
        print "Error: argument --snrFilter must be 0, 1, or 2. (%i)" % args.snrFilter
        sys.exit(1)

    snrThresh = args.snrThreshold
    if len(snrThresh) < 2: #at least one arg is requried (len >= 1)
        snrThresh.append(snrDefault)
    elif len(snrThresh) > 2: #discard extra args
        snrThresh = snrThresh[:2]
    if snrThresh[0] < 0 or snrThresh[0] > 30:
        print "Error: first argument to --snrThreshold must be >= 0 AND <= 30. (%f)" % snrThresh[0]
        sys.exit(1)
    if snrThresh[1] < 0 or snrThresh[1] > 30:
        print "Error: second argument to --snrThreshold must be >= 0 AND <= 30. (%f)" % snrThresh[1]
        sys.exit(1)
    args.snrThreshold = snrThresh

    if args.color < 0 or args.color > 2:
        print "Error: argument --color must be >= 0 AND <= 2. (%i)" % args.color
        sys.exit(1)

    return args

def main():

    args = get_args()
    
    bnx = args.mol
    
    #Step 1: threshold labels based on SNR (optional)
    if args.dosnrThresh:

        if args.snrFilter == 0:
            #need to add second snrThreshold: see comments at fn def
            bnx = run_snr_filtering(os.path.join(args.output,args.prefix), args.mol, args.pipeline, args.ra, args.snrThreshold[0])

        elif args.snrFilter == 1:
            bnx = run_snr_filteringRA(os.path.join(args.output,args.prefix), args.mol, args.pipeline, args.ra, args.snrThreshold)
        
        if not bnx or not util.checkFile(bnx, ".bnx"): #should be path to filtered bnx, if not present, exit
            print "Error: SNR filter failed"
            return

    #Step 2: estimate error parameters using autoNoise in Pipeline
    scaled_bnx = run_autonoise(args.output, args.prefix, bnx, args.ref, args.pipeline, args.ra, args.nthreads,
                               args.optArgs, args.color)

    if not util.checkFile(scaled_bnx, ".bnx"):
        print "Error: autonoise result bnx file not found, check autonoise job: " + scaled_bnx
        sys.exit(1)

    #Step 3: align molecule to reference using runAlignMol.py
    run_alignmol(args.output, scaled_bnx, args.ref, args.pipeline, args.ra, args.nthreads, args.optArgs, args.color)
    
    #Step 4: compress output
    compress_output(args.output)

if __name__ == '__main__':
    main()
