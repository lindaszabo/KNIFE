# given an existing comprehensive junction fasta file, limit it based on specified criteria.
# Current options include limiting to only reg, dup, or rev and specifying only within gene,
# all, or only when genes are different. Typically want to call for an entire directory, but
# it is implemented so you can call on just a specific file too.

# usage 1: python limitFasta.py -d /srv/gsfs0/projects/salzman/Linda/junctionIndex/fastas/fasta1Mb
#                 -o /srv/gsfs0/projects/salzman/Linda/junctionIndex/fastas/fasta1Mb_scrambled
#                 -t rev -p _scrambled -v 
# usage 2: python limitFasta.py -s /srv/gsfs0/projects/salzman/Linda/junctionIndex/fasta1Mb/chr10_junctions.fa
#                 -o /srv/gsfs0/projects/salzman/Linda/junctionIndex/fasta1Mb_reg_inGene
#                 -t reg -b within -p _reg_inGene -v

from Bio import SeqIO
import argparse
import re
import os
import utils_os

# we are basing whether we include the junction on
# whether the 2 sides of the junctions are within the same gene
# and whether it is a dup, reg, or rev junction
# return True if the junction meets criteria and should be included in the output
def inBounds(match):
    if args.juncType != match.group(3):
        return False
    if args.bounds == 'all':
        return True
    if args.bounds == 'within':
        return match.group(1) == match.group(2)
    else:
        return match.group(1) != match.group(2)

    
    
# param fileName: full path to template fasta file  
def writeLimitedFasta(fileName):
        
    fileBase = os.path.splitext(os.path.basename(fileName))[0]  # name of file without extension
    
    if args.verbose:
        print fileName
        print fileBase
        print args.outDir + "/" + fileBase + args.postpend
    
    in_handle = open(fileName, "rU")
    out_handle = open(args.outDir + "/" + fileBase + args.postpend + ".fa", "wb")
    
    fasta_sequences = SeqIO.parse(in_handle,"fasta")
    for fasta in fasta_sequences:
        match = id_patt.search(fasta.id)
        if match and inBounds(match):
            SeqIO.write(fasta, out_handle, "fasta")
    
    out_handle.close()
    in_handle.close()

if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--fastaDir', help='directory containing fasta files to limit. Either -d or -s required.')
    parser.add_argument('-s', '--singleFile', help='fasta file to limit. Either -d or -s required.')
    parser.add_argument('-o', '--outDir', help='directory to output files. Will be created if it does not exist.', required=True)
    parser.add_argument('-t', '--juncType', help='type of junction to limit to', choices=['reg', 'rev', 'dup'], required=True)
    parser.add_argument('-b', '--bounds', help='type of junction to limit to', choices=['all', 'within', 'between'], default='all')
    parser.add_argument('-p', '--postpend', help='text to add to end of output file names', required=True)
    parser.add_argument('-v', '--verbose', help='print info about data obtained', action='store_true')
    
    args = parser.parse_args()
    
    if args.fastaDir and args.singleFile:
        sys.exit("Only 1 of -d and -s can be specified")
    
    # create output directory if it does not exist    
    utils_os.createDirectory(args.outDir)
    
    # fasta junction ids look like chr10|TTC40:134751179|MYC1:134722640|reg|-
    # where we are interested in TTC40, MYC1, and reg
    id_patt = re.compile(".+?\|(.+?):.+?\|(.+?):.+?\|(.+?)\|.*")
    
    if args.fastaDir:
        # loop through files directory
        for f in os.listdir(args.fastaDir):
            if f.endswith(".fa"): # only parse if this is a fasta file
                writeLimitedFasta('/'.join([args.fastaDir,f]))
    elif args.singleFile:
        writeLimitedFasta(args.singleFile)
        