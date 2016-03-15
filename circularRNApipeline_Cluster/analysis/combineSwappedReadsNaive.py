# Before running this script, findCircularRNA.sh must have been run twice. The first
# time is as described in the main project README. For the second run, create a new directory
# for the reads and re-label the read1 file to read2 and vice versa. Run the pipeline again
# on these swapped reads. 

# After circularRNApipeline has been run on original files and swapped files (where Rd1 was renamed Rd2),
# this script consolidates the naive report output from the 2 files into a single report. If a read was 
# assigned to the same junction based on R1 and R2, it will not be double-counted. If unaligned mode was run,
# consolidated denovo reports will be created as well.

# This is aimed at being able to directly compare our results to others who count junctional reads from
# both R1 and R2. We report the p-value reported by each run, and the sum of all reads
# from both runs, as long as they are not the same id.

# parameters:
# -a is path to original run output (directory that contains reports, ids, etc subdirectories)
# -b is path to swapped run output (directory that contains reports, ids, etc subdirectories)
# -q is read id style (complete or appended, same as used in pipeline run)
# -v is verbose flag to output processing status to terminal

# output: a directory called combinedReports is created as a subdirectory of the path specified as -a (original run output directory).
#         separate files are created for circular and linear junctions for each of the sample files

# usage: python combineSwappedReadsGLM.py -a /home/linda/alignments/cancerData/circReads
#                                         -b /home/linda/alignments/cancerDataSwapped/circReads
#                                         -q complete
#                                         -v
import argparse
import os
import utils_os
import re

patt_filename = re.compile(".+?_report.txt")


class junctionAlignment:
    
    def __init__(self, id, numGreads, numBreads, numUreads, pval):
        self.id = id  # junction id
        self.numGreads = int(numGreads) # number of good reads (linear or circular)
        self.numBreads = int(numBreads) # number of bad reads (anomaly or decoy)
        self.numUreads = int(numUreads) # number of unmapped reads
        self.pval = pval # posterior probability of
        self.numGreads2 = self.numBreads2 = self.numUreads2 = self.pval2 = 0  # will be updated if we get more in the 2nd run
        self.numSamples = 1 # how many runs it was detected in
        self.reads = {}  # read ids that supported this junction
    
    def __str__(self):
        return "\t".join([self.id, str(self.numGreads), str(self.numBreads), str(self.numUreads), str(self.pval),
                          str(self.numGreads2), str(self.numBreads2), str(self.numUreads2),
                          str(self.pval2), str(len(self.reads))])   

def populateJunctions(fileName, fileType):
    handle = open(fileName, "rU")
    
    if fileType == "A":
        if args.verbose:
            print "processing file A", fileName
        for line in handle:
            # first time populating junctions so just add all we see
            if not line.startswith("@"):
                vals = line.strip().split()
                if "|reg|" in vals[0]:
                    # linear junction
                    junctions[vals[0]] = junctionAlignment(vals[0], vals[1], vals[2], vals[3], vals[7])  # junction name, linear, anomaly, unmapped, pval
                else:
                    # circular junction
                    junctions[vals[0]] = junctionAlignment(vals[0], vals[5], vals[6], vals[3], vals[7])  # junction name, circ, decoy, unmapped, pval
        if args.verbose:
            print "processed junctions", len(junctions)
    elif fileType == "B":
        if args.verbose:
            print "processing file B", fileName
        for line in handle:
            if not line.startswith("@"):
                vals = line.strip().split()
                # update or create new and then update
                if vals[0] not in junctions:
                    junctions[vals[0]] = junctionAlignment(vals[0], 0, 0, 0, 0)  # create empty junction since it wasn't in first report
                if "|reg|" in vals[0]:
                    # linear junction
                    junctions[vals[0]].numGreads2 = int(vals[1])
                    junctions[vals[0]].numBreads2 = int(vals[2])
                else:
                    # circular junction
                    junctions[vals[0]].numGreads2 = int(vals[5])
                    junctions[vals[0]].numBreads2 = int(vals[6])
                junctions[vals[0]].numUreads2 = int(vals[3])
                junctions[vals[0]].pval2 = vals[7]
        if args.verbose:
            print "processed junctions", len(junctions)
    
    handle.close()

# read in id files and assign reads to junctions    
def updateReadCounts(fileName):
    handle = open(fileName, "rU")
    
    for line in handle:
        id, category, posR1, qualR1, aScoreR1, numNR1, readLenR1, junctionR1, strandR1, posR2, qualR2, aScoreR2, numNR2, readLenR2, junctionR2, strandR2 = line.strip().split()

        if category.startswith("linear") or category.startswith("circ"):
            if args.fastqIdStyle == "appended":
                id = id[:-1]

            if junctionR1 in junctions:
                # update reads for read 1 aligning to junction
                junctions[junctionR1].reads[id] = None
            if junctionR2 in junctions:
                # update reads for read 2 aligning to junction
                junctions[junctionR2].reads[id] = None    
    handle.close()
                
# param fileType: a or b, matching with whether the file came from aDir or bDir
def combineResults(rFileName, iFileName):
    
    populateJunctions("/".join([args.dirA, "reports", rFileName]), "A")  # read in results from 1st reports file to populate objects
    populateJunctions("/".join([args.dirB, "reports", rFileName]), "B")  # update with results from 2nd reports file
    
    updateReadCounts("/".join([args.dirA, "ids", iFileName]))  # update with read ids from 1st report file             
    updateReadCounts("/".join([args.dirB, "ids", iFileName]))  # update with read ids from 2nd report file 
    
    
    # print out
    handle = open("/".join([args.dirA, "combinedReports", "naive"+rFileName]), "wb")
    handle.write("junction\torig_circOrLinear\torig_decoyOrAnom\torig_unmapped\torig_pval\t")
    handle.write("swapped_circOrLinear\tswapped_decoyOrAnom\tswapped_unmapped\tswapped_pval\ttotal_reads\n")
    for j in junctions:
        handle.write(str(junctions[j]) + "\n")
    handle.close()
        
        
    
if __name__  == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--dirA', help='directory containing original pipeline reports (parent of reports, ids, etc)', required=True)
    parser.add_argument('-b', '--dirB', help='directory containing swapped pipeline reports (parent of reports, ids, etc directories)', required=True)
    parser.add_argument('-q', '--fastqIdStyle', help='type of read ids used', required=True, choices=['appended', 'complete'])
    parser.add_argument('-v', '--verbose', help='print extra debugging info', action='store_true')
    args = parser.parse_args()
    
    utils_os.createDirectory("/".join([args.dirA, "combinedReports"]))
    
    for f in os.listdir("/".join([args.dirA, "reports"])):
        reportFileName = os.path.basename(f)
        if args.verbose:
            print reportFileName

        if patt_filename.search(reportFileName):
            idFileName = reportFileName.replace("report", "_output")  # could be denovo or annotated report file, all handled the same
            if args.verbose:
                print reportFileName
                print idFileName
            
            junctions = {}
            combineResults(reportFileName, idFileName)
    