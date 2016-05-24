# Takes a pickled exon dictionary and creates a fasta file containing all of the 
# possible exon pair junctions (x-x, x-y, y-x) within a sliding window. Default
# window size is 100Kb.
#
# Since I was running into memory issues when creating junctions for 1Mb window
# and storing the junctions as a pkl, I reduced the memory footprint
# by just keeping track of keys and writing out the junction object to file if it
# does not exist in list of keys already. Also sped it up significantly by using
# a dictionary instead of a list for keeping track of previously seen junctions.
#
# By default it loops through all files in a directory, but if you pass -s filename
# then it will just run for that single file. In general, I recommend using it in
# single mode, then you can get the files for the entire genome in just a few hours.

# usage 1 (whole directory): python makeJunctionsAndWriteFasta.py -w 2000 -e output/exons -r output/records -f output/fasta -v
# usage 2 (single file): python makeJunctionsAndWriteFasta.py -w 2000 -s exon.pkl -r output/records -f output/fasta -v

#        all arguments are optional

### TODO:
### 1. add error checking and exception handling
### 4. streamline code
### 5. add logging instead of printing to stdout?

import argparse
from Bio import SeqIO
from itertools import product
import cPickle as pickle
import os
import utils_os
import re
import sys
from utils_junction import junction

EMPTY_SEQUENCE = "NOSEQUENCE" # to indicate junction sequence contains all Ns so do not include in fasta output
        
# pattern like exonsByStrand_chr1.pkl
# where we are interested in chr1 as the file id
patt_exonfilename = re.compile(".*exonsByStrand_(.+?)\.pkl")

# file name like exonsByStrand_chr1.pkl
patt_exonfile = re.compile("exonsByStrand_.+\.pkl")

# Helper function called to actually do the writing of the junction sequence to
# the fasta file. First checks to make sure this junction has not already been written.
#
# param curSeqRec: sequenceRecord object that contains the sequence for this chromosome
# param id: something like chr# used to identify the current chromosome and all related files
# param allExons: master dictionary of all exons for this chromosome 
# param curExons: list of exon keys for the exons within the current window
# param junctions: running dict of all junction ids already observed
# param handle: open file handle to write fasta entries to
#
# return: junctions dict with the new junction ids written in this round appended
def writeAllPairs(curSeqRec, id, allExons, curExons, junctions, handle):
    allPairs = product(curExons, repeat=2)  # create all exon pairs including x-x, x-y, y-x
    for a,b in allPairs:
        try:
            exonA = allExons[a] 
            exonB = allExons[b]
            if exonA.strand == 1:
                curJuncId = (exonA.location.end, exonB.location.start)
            else:
                curJuncId = (exonA.location.start, exonB.location.end)
            
            # if this junction is not already accounted for, write it to the file
            if not curJuncId in junctions:
                curJunc = junction(curSeqRec, id, exonA, exonB)
                junctions[curJuncId] = None  # don't need to store anything in the dictionary, just need the keys hashed
                if str(curJunc) != EMPTY_SEQUENCE:
                    handle.write(str(curJunc.printHeader(args.name1, args.name2)))  # junction fasta header
                    handle.write(str(curJunc))  # prints out the sequence
        except Exception as e:
            print "Exception"
            print e
            print "error:", sys.exc_info()[0]
            print "parsing features for", a, b
            
    return junctions

# Move along the sliding window, get all exons within the window, write out all of those
# junctions, then move along the chromosome by 1 exon and repeat until you make it
# to the end of the chromosome. To be considered within a window, both exons must be
# entirely included in the window. In each window, some of the junctions have already
# been written because both exons were in a previous window. In this case the junction
# is ignored so it will only ever be written to the file once.
#
# This function goes through the + strand and - strand exons separately, all are written to
# the same fasta file though.
# 
# param fileId: something like chr# which is contained in names of all related files for this chromosome
# param exonFile: name of the exon file for this chromosome
def createJunctions(fileId, exonFile):
    if args.verbose:
        print fileId
    
    # get the sequence record for this chromosome (exons just contain annotation data)
    exonSeqRec = pickle.load(open(args.recordDir + '/rec_' + fileId + '.pkl', 'rb'))
    
    # and get the exons that belong to this sequence
    chrExonsByDirection = pickle.load(open(args.exonDir + '/' + exonFile, 'rb'))
    
    allJunctions = {} # (endExon1, startExon2) identifier to keep track of which junctions we have done already 
    
    if args.verbose:
        print len(chrExonsByDirection)
    
    outf = open(args.fastaDir + '/' + fileId + '_junctions.fa', 'wb')
    
    # there should be 2 elems in each file, 1 for + strand exons and 1 for - strand exons
    for elem in chrExonsByDirection:
        chrId, strand, exons = elem
        
        exonsInRange = exons.keys()
        
        if len(exonsInRange) > 0:
            if strand == 1:
                # start from begining and work forward for + strand
                startPos = min(x[0] for x in exonsInRange)
                # when we are looking at a set of exons including this last one, we are done with the sliding window
                endPos = max(x[1] for x in exonsInRange)
                
                while True:
                    # the starting point was already established, now limit to those within window
                    exonsInRange = [x for x in exonsInRange if x[1] <= startPos + args.window]
                    
                    if len(exonsInRange) > 0:
                        # write these junctions to file if they haven't already been written
                        allJunctions = writeAllPairs(exonSeqRec, chrId, exons, exonsInRange, allJunctions, outf)
                        
                        # if max in exonsInRange is endPos then we're done
                        if max(x[1] for x in exonsInRange) == endPos:
                            break
                        # otherwise move over 1 exon and repeat
                        exonsInRange = [x for x in exons.keys() if x[0] > startPos]
                        startPos = min(x[0] for x in exonsInRange)   
            else:
                # start from end and work backwards for - strand
                startPos = max(x[1] for x in exonsInRange)
                endPos = min(x[0] for x in exonsInRange)
                
                while True:
                    # the starting point was already established, now limit to those within window
                    exonsInRange = [x for x in exonsInRange if x[0] >= startPos - args.window]
                    if len(exonsInRange) > 0:
                        allJunctions = writeAllPairs(exonSeqRec, chrId, exons, exonsInRange, allJunctions, outf)
                        
                        # if min in exonsInRange is endPos then we're done
                        if min(x[0] for x in exonsInRange) == endPos:
                            break
                        # otherwise move over 1 exon and repeat
                        exonsInRange = [x for x in exons.keys() if x[1] < startPos]
                        startPos = max(x[1] for x in exonsInRange)
                
    outf.close()
    
    if args.verbose:
        print len(allJunctions)


if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--window', help='sliding window size to create junctions', default=1000000, type=int)
    parser.add_argument('-e', '--exonDir', help='directory containing exon pickle files', default='output/exons')
    parser.add_argument('-s', '--singleFile', help='path to single exon file that should be parsed for junctions')
    parser.add_argument('-r', '--recordDir', help='directory containing seq record pickle files', default='output/records')
    parser.add_argument('-f', '--fastaDir', help='directory to output junction fasta files, will be created if does not exist', default='output/fasta')
    parser.add_argument('-n1', '--name1', help='name of field in gtf to use for gene names', default='gene_name')
    parser.add_argument('-n2', '--name2', help='name of field in gtf to use for gene names if n1 does not exist', default='gene_id')
    parser.add_argument('-v', '--verbose', help='print extra debugging info', action='store_true')
    args = parser.parse_args()

    # create fastaDir if it doesn't exist
    utils_os.createDirectory(args.fastaDir)
    
    # only run for a single file
    if args.singleFile:
        if args.verbose:
            print "running for single file", args.singleFile
        exonObj = os.path.basename(args.singleFile)
        fileId = utils_os.getFileId(patt_exonfilename, 1, args.singleFile)  # usually something like chr# which is in the name of each created file
        createJunctions(fileId, exonObj)
    else:
        if args.verbose:
            print "running for directory", args.exonDir
        # or loop through files in exonDir to create junction file for each
        for exonObj in os.listdir(args.exonDir):
            if patt_exonfile.search(exonObj): # only parse if this is an exon pickled file
                fileId = utils_os.getFileId(patt_exonfilename, 1, exonObj)  # usually something like chr# which is in the name of each created file
                createJunctions(fileId, exonObj)
    

    
    
    

    
    
 

