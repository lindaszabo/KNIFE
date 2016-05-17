# Creates a fasta and fastq file containing all reads for the sample that did not align to any of
# the Bowtie2 indices in the directory "orig/unaligned". It also creates a fastq file under
# "orig/unaligned/forDenovoIndex" that contains the subset of unaligned reads that are long
# enough to be used for the de novo index (NTRIM + 8). This forDenovoIndex file had to be
# created because Bowtie does not gracefully handle when a read gets trimmed to a negative length.
#
# This script is called by qualityStatsSingleSample.sh when the sampleStats reports are being generated,
# so it returns the number of unaligned reads to be reported. This is the total number of unaligned reads,
# not just those that were long enough to be used for the de novo index.

import argparse
import os
import gzip
from ParseFastQ import ParseFastQ

def addAlignedIds(samFile):
    handle = open(samFile, "rU")

    for line in handle:
        if not line.startswith("@"): # ignore header lines
            aligned[line.strip().split()[0]] = None  #
    
    handle.close()
    
if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--readFile', help='path to fastq file containing reads you attempted to align', required=True)
    parser.add_argument('-n', '--fileBaseName', help='base name used in alignment output files', required=True)
    parser.add_argument('-a', '--alignmentDir', help='path to directory containing alignments output for this dataset (will end with /orig)', required=True)
    parser.add_argument('-t', '--ntrim', help='ntrim used in denovo pipeline', type=int, required=True)
    parser.add_argument('-v', '--verbose', help='print info about data obtained', action='store_true')
    args = parser.parse_args()
    
    aligned = {}
    unaligned = {}
    
    minReadLen = args.ntrim + 8  # we want to have at least 8 nt left after read is trimmed for denovo split reads
    
    try:
        # put each id from sam files in a dict
        addAlignedIds("".join(["/".join([args.alignmentDir, "genome", args.fileBaseName]), "_genome_output.sam"]))
        addAlignedIds("".join(["/".join([args.alignmentDir, "junction", args.fileBaseName]), "_junction_output.sam"]))
        addAlignedIds("".join(["/".join([args.alignmentDir, "reg", args.fileBaseName]), "_reg_output.sam"]))
        addAlignedIds("".join(["/".join([args.alignmentDir, "ribo", args.fileBaseName]), "_ribo_output.sam"]))
        
        # get each id from fastq file, add it to a dict if not in the aligned dict
        parser = ParseFastQ(args.readFile)
        for (seqHeader, seqStr, qualHeader, qualStr) in parser:
            title = seqHeader.split()[0][1:]
            if title not in aligned:
                unaligned[title] = (seqHeader, seqStr, qualHeader, qualStr)

        # write list of unaligned ids to fasta and fastq files
        out_handle = open(args.alignmentDir + "/unaligned/unaligned_" + args.fileBaseName + ".fasta", "wb")
        out_fq_handle = open(args.alignmentDir + "/unaligned/unaligned_" + args.fileBaseName + ".fq", "wb")
        out_denovo_handle = open(args.alignmentDir + "/unaligned/forDenovoIndex/unaligned_" + args.fileBaseName + ".fq", "wb") 
        for id in unaligned.keys():
            # fasta file
            out_handle.write(">" + id + "\n")
            out_handle.write(unaligned[id][1] + "\n")
            # fastq file
            out_fq_handle.write("%s\n%s\n%s\n%s\n" % (unaligned[id][0], unaligned[id][1], unaligned[id][2], unaligned[id][3]))
            # fastq file for denovo index
            if len(unaligned[id][1]) >= minReadLen:
                out_denovo_handle.write("%s\n%s\n%s\n%s\n" % (unaligned[id][0], unaligned[id][1], unaligned[id][2], unaligned[id][3]))
            
            
        out_handle.close()
        out_fq_handle.close()
        out_denovo_handle.close()
        
        # return size of unaligned dict
        print len(unaligned)
    except Exception as e:
        print "Exception:", e