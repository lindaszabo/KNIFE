# Before running this script, findCircularRNA.sh must have been run twice. The first
# time is as described in the main project README. For the second run, create a new directory
# for the reads and re-label the read1 file to read2 and vice versa. Run the pipeline again
# on these swapped reads. 

# After circularRNApipeline has been run on original files and swapped files (where Rd1 was renamed Rd2),
# this script consolidates the GLM output from the 2 files into a single report. If a read was 
# assigned to the same junction based on R1 and R2, it will not be double-counted.

# This is aimed at being able to directly compare our results to others who count junctional reads from
# both R1 and R2. We report the posterior probability reported by each run, and the sum of all reads
# from both runs, as long as they are not the same id.

# parameters:
# -a is path to original run output (directory that contains glmReports, ids, etc subdirectories)
# -b is path to swapped run output (directory that contains glmReports, ids, etc subdirectories)
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

patt_GLMfilename = re.compile(".+?__.*JuncProbs.txt")


class junctionAlignment:
    
    def __init__(self, id, numReads, posterior):
        self.id = id  # junction id
        self.numReads = int(numReads) # number of reads
        self.numReads2 = 0  # will be updated if we get more in the 2nd run
        self.posterior = float(posterior) # posterior probability of 
        self.posterior2 = 0 # will be updated if we get more in the 2nd run
        self.numSamples = 1 # how many runs it was detected in
        self.reads = {}  # read ids that supported this junction
    
    def __str__(self):
        return "\t".join([self.id, str(self.numReads), str(self.posterior), str(self.numReads2),
                          str(self.posterior2), str(len(self.reads))])   

def populateJunctions(fileName, fileType):
    handle = open(fileName, "rU")
    
    if fileType == "A":
        if args.verbose:
            print "processing file A", fileName
        for line in handle:
            if not line.startswith("junction"):
                vals = line.strip().split()
                # first time populating junctions
                junctions[vals[0]] = junctionAlignment(vals[0], vals[1], vals[2])  # junction name, count, posterior
        if args.verbose:
            print "processed junctions", len(junctions)
    elif fileType == "B":
        if args.verbose:
            print "processing file B", fileName
        for line in handle:
            if not line.startswith("junction"):
                vals = line.strip().split()
                # update or create new and then update
                if vals[0] not in junctions:
                    junctions[vals[0]] = junctionAlignment(vals[0], 0, 0)  # create empty junction since it wasn't in first report
                junctions[vals[0]].numReads2 = int(vals[1])
                junctions[vals[0]].posterior2 = float(vals[2])
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
    
    populateJunctions("/".join([args.dirA, "glmReports", rFileName]), "A")  # read in results from 1st reports file to populate objects
    populateJunctions("/".join([args.dirB, "glmReports", rFileName]), "B")  # update with results from 2nd reports file
    
    updateReadCounts("/".join([args.dirA, "ids", iFileName]))  # update with read ids from 1st report file             
    updateReadCounts("/".join([args.dirB, "ids", iFileName]))  # update with read ids from 2nd report file 
    
    
    # print out
    handle = open("/".join([args.dirA, "combinedReports", rFileName]), "wb")
    handle.write("junction\torig_count\torig_posterior\tswapped_count\tswapped_posterior\ttotal_reads\n")
    for j in junctions:
        handle.write(str(junctions[j]) + "\n")
    handle.close()
        
        
    
if __name__  == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--dirA', help='directory containing original pipeline reports (parent of glmReports, ids, etc)', required=True)
    parser.add_argument('-b', '--dirB', help='directory containing swapped pipeline reports (parent of glmReports, ids, etc directories)', required=True)
    parser.add_argument('-q', '--fastqIdStyle', help='type of read ids used', required=True, choices=['appended', 'complete'])
    parser.add_argument('-v', '--verbose', help='print extra debugging info', action='store_true')
    args = parser.parse_args()
    
    utils_os.createDirectory("/".join([args.dirA, "combinedReports"]))
    
    for f in os.listdir("/".join([args.dirA, "glmReports"])):
        reportFileName = os.path.basename(f)
        if args.verbose:
            print reportFileName

        if patt_GLMfilename.search(reportFileName):
            idFileName = reportFileName.replace("circJuncProbs", "output").replace("linearJuncProbs", "output")  # could be linear or circular file
            if args.verbose:
                print reportFileName
                print idFileName
            
            junctions = {}
            combineResults(reportFileName, idFileName)
    