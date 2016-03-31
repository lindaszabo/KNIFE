
# Assumes preprocessAlignedReads.sh has been called and files with read ids for genome,
# junction, reg, ribo, and transcriptome alignments have been output.

# Loops through all id files to find those Rd1 ids that are aligned to a junction and not to the ribosome
# or genome. Then goes back to original sam files to get data for reads that aligned to each of these junctions
# in order to create JuncObjs with all data for the reads that aligned to a junction.

# Then assigns each read to a bucket based on where Rd1 and Rd2 aligned. Outputs reports based on analysis using the
# Naive method described in the paper into the reports directory. A very early method using a hard cutoff on alignment
# scores is still implemented here (passing alignment score thresholds for read1 as -a1 and read2 as -a2),
# but should not be used as the naive method performs better. 

import argparse
import os
import utils_os
from utils_juncReads_minimal import *
from scipy.stats import scoreatpercentile
from scipy.stats import poisson
from math import ceil
import sys

JUNC_MIDPOINT = 150  

# look in the ribo id file for this sample, and in the associated genome id file for this
# sample to find the read ids that aligned to the genome or the ribosome.
# Then look in the associated regular junction overlapped read id file and select
# only those ids that are junction overlapped and did not align to ribosome or genome.
# Populate the global regIds with these ids. Then look in the associated all-junction
# overlapped id file and populate the global nonRegIds variable with these read ids
# that did not align to the ribosome, genome, or regular junctions
def selectCandidateIds():
    
    # add the suffix if we used an alternate overlap specification and output the overlapped ids to an alternate location
    regIdDir = "reg" + args.junctionIdDirSuffix
    juncIdDir = useJuncStr + args.junctionIdDirSuffix
        
    idDir = "/".join([args.parentDir, "orig", "ids"])
    
    # only need to look for ids to exclude if we are not looking at previously unaligned reads
    if not args.unalignedMode:
        # get ribo ids
        try:
            handle = open("".join(["/".join([idDir, "ribo", args.sampleId]), "_ribo_output.txt"]), "rU")
            if args.fastqIdStyle == "appended":
                for line in handle:
                    ignoreIds[line.strip().split()[0][:-1]] = None
            else:
                for line in handle:
                    ignoreIds[line.strip().split()[0]] = None
            handle.close()
        except Exception as e:
            print "Exception"
            print e
            print "parsing ribo ids for", line
        
        # load genome aligned id file for the same sample and add to ignoreIds
        handle = open("".join(["/".join([idDir, "genome", args.sampleId]), "_genome_output.txt"]), "rU")
            
        if args.fastqIdStyle == "appended":
            for line in handle:
                try:
                    ignoreIds[line.strip().split()[0][:-1]] = None
                except TypeError as e:
                    print "error parsing genome ids for", line
                    print "Type error({0}): {1}".format(e.errno, e.strerror)
                    print "error:", sys.exc_info()[0]
                except:
                    print "error parsing genome ids for", line
                    print "error:", sys.exc_info()[0] 
        else:
            for line in handle:
                ignoreIds[line.strip().split()[0]] = None       
        handle.close()
        
        
        # load reg-junction aligned id file for same sample and load ids not in ignoreIds
        try:
            handle = open("".join(["/".join([idDir, regIdDir, args.sampleId]), "_reg_output.txt"]), "rU")
            for line in handle:
                if args.fastqIdStyle == "appended":
                    testId = line.strip().split()[0][:-1]
                else:
                    testId = line.strip().split()[0]
                    
                if testId not in ignoreIds:
                    regIds[testId] = None
            handle.close()
        except Exception as e:
            print "Exception"
            print e
            print "parsing reg ids for", line
            
    
    # load junction aligned id file (or de novo aligned id file) for same sample and load ids not in ignoreIds or regJuncIds
    try:
        handle = open("_".join(["/".join([idDir, juncIdDir, args.sampleId]), useJuncStr, "output.txt"]), "rU")
        for line in handle:
            if args.fastqIdStyle == "appended":
                testId = line.strip().split()[0][:-1]
            else:
                testId = line.strip().split()[0]
                
            if testId not in ignoreIds and testId not in regIds:
                nonRegIds[testId] = None   
        handle.close()
    except Exception as e:
        print "Exception"
        print e
        print "error parsing junction ids for", line 
    
    # print out these ids for later debugging use
    try:
        if args.unalignedMode:
            out_handle = open("".join(["/".join([idDir, "denovoNonGR", args.sampleId]), "_output.txt"]), "wb")
        else:
            out_handle = open("".join(["/".join([idDir, "juncNonGR", args.sampleId]), "_output.txt"]), "wb")
        
        for i in nonRegIds:
            out_handle.write(i)
            out_handle.write("\n")
        for i in regIds:
            out_handle.write(i)
            out_handle.write("\n")
        out_handle.close()
    except Exception as e:
        print "Exception"
        print e
        print "writing ids for", i

      
# Get mismatch rate for the decoys in this dataset. This is total_mismatches / total_bases
def getDecoyMismatchRate():
    
    aScore = 0  # sum of all alignment scores
    numBases = 0
    
    for j in junctions:
        for r in junctions[j].decoyReads:
            aScore += r.readStat
            numBases += r.juncRead.readLen
            # this is average aScore, so it is aScore for read 1 or avg for read1 and read2 if mate is not None
            if r.useMate != None:
                aScore += r.readStat  # need to double the contribution of aScore since that was the avg of the 2 mates
                numBases += r.useMate.readLen  # and need to add the length of the mate
        
    if numBases == 0:
        return None
    else:
        return (aScore / -6.0) / numBases  # this is the mismatch rate per base observed in the decoys

# param alignedReads: array of juncReadObj that aligned to the junction 
def getPval(alignedReads):
    if len(alignedReads) > 0:
        useMMrate = globalDecoyMMrate
        
        # total number of mismatches observed for all reads aligning to this junction, rounded to get integer which is required for poisson.cdf
        num_mm = int(ceil(sum([x.readStat for x in alignedReads]) / -6.0))  
        num_bases = sum([x.juncRead.readLen for x in alignedReads]) + sum([x.useMate.readLen for x in alignedReads if x.useMate != None])
        
        return 1 - poisson.cdf(num_mm, useMMrate*num_bases)
    else:
        return "-"

def getReadScores(alignedReads):
    scores = []
    if len(alignedReads) > 0:
        a1Scores=[int(r.juncRead.aScore) for r in alignedReads]
        if args.singleEnd:
            scores = a1Scores
        else:
            a2Scores=[int(r.useMate.aScore) for r in alignedReads]
            scores = zip(a1Scores, a2Scores)
    
    return scores

# per-junction p-value is calculated assuming number of mismatches is Poisson(0.01*avg_read_length)
# .01 is the high end of the Illumina sequencing error rate
# prints all junctions, even those that have no circular reads
def reportCircularReads(cutoff):
    
    if args.unalignedMode:
        report_handle = open("_".join(["/".join([args.parentDir, args.outDirName, "reports", args.sampleId]), "denovo_report.txt"]), "wb")
    else:
        report_handle = open("_".join(["/".join([args.parentDir, args.outDirName, "reports", args.sampleId]), "report.txt"]), "wb")
        
    if cutoff:
        report_handle.write("".join(["@Global alignment score cutoff: ", str(cutoff), "\n"]))
    else:
        report_handle.write("".join(["@GUsing Poisson distribution with p=: ", str(globalDecoyMMrate), "\n"]))
    report_handle.write("@junction\tlinear\tanomaly\tunmapped\tmultimapped\tcirc\tdecoy\t")
    
    # if we specified a read cutoff, the value returned is an FDR, otherwise it is a p-value using naive model
    if cutoff:
        report_handle.write("FDR")
    else:
        report_handle.write("pvalue")
    
    report_handle.write("\tscores\n")
    
    for j in junctions:
        # print out the junction and read id for all good reads
        # either both mates need to pass threshold, or single-end read needs to pass Rd1 threshold
        if cutoff:
            numGood = 0  # track any linear or circular that passed the score threshold
            for elem in junctions[j].circularReads:
                if int(elem.juncRead.aScore) >= int(cutoff[0]) and (args.singleEnd or int(elem.useMate.aScore) >= int(cutoff[1])):
                    numGood += 1
                
            for elem in junctions[j].linearReads:
                if int(elem.juncRead.aScore) >= int(cutoff[0]) and (args.singleEnd or int(elem.useMate.aScore) >= int(cutoff[1])):
                    numGood += 1 
                
            # print out global junction reads stats
            numCandidates = len(junctions[j].circularReads) + len(junctions[j].linearReads)
            numBad = numCandidates - numGood
            if numCandidates > 0:
                sig_stat = float(numBad) / numCandidates
            else:
                sig_stat = "-"
        elif "|reg|" not in junctions[j].id:  # if it's a circular junction
            sig_stat = getPval(junctions[j].circularReads)  # p-value
        elif "|reg|" in junctions[j].id:  # or it is a linear junction
            sig_stat = getPval(junctions[j].linearReads)  # p-value
        else:
            print "cound not get sig stat for", junctions[j].id
            sig_stat = "NA"
        
        report_handle.write(str(j) + "\t")
        report_handle.write(str(len(junctions[j].linearReads)) + "\t")
        report_handle.write(str(len(junctions[j].anomalyReads)) + "\t")
        report_handle.write(str(len(junctions[j].unmappedReads)) + "\t")
        report_handle.write(str(len(junctions[j].multimappedReads)) + "\t")
        report_handle.write(str(len(junctions[j].circularReads)) + "\t")
        report_handle.write(str(len(junctions[j].decoyReads)) + "\t")
        report_handle.write(str(sig_stat) + "\t")
        
        if "|reg|" not in junctions[j].id:  # if it's a circular junction
            scores = getReadScores(junctions[j].circularReads)  
        elif "|reg|" in junctions[j].id:  # or it is a linear junction
            scores = getReadScores(junctions[j].linearReads)  # scores to print in report
        else:
            scores = []
            
        if len(scores) > 10:
            # if we have lots of scores for this junction, just write out the quantiles
            for i in xrange(0,101,10):    
                if args.singleEnd:
                    report_handle.write(str(scoreatpercentile(scores, i, interpolation_method='lower')) + ",")
                else:
                    report_handle.write(str(scoreatpercentile([x[0] for x in scores], i, interpolation_method='lower')) + ":" +
                                        str(scoreatpercentile([x[1] for x in scores], i, interpolation_method='lower')) + ",")
        else:
            # otherwise we have just a few so let's just print them all out
            for s in scores:
                report_handle.write(str(s) + ",")
        
        report_handle.write("\n")
    
        # store FDR for this junction so we can decide if we think it is an artifact or not
        if sig_stat != "-" and sig_stat != "NA":
            junctions[j].fdr = sig_stat  # will need to update this variable name in juncObj if we are using p-values in the future instead of FDRs
        
    report_handle.close()
            
def reportAllReadIds2(cutoff):
    if args.unalignedMode:
        id_handle = open("_".join(["/".join([args.parentDir, args.outDirName, "ids", args.sampleId]), "denovo__output.txt"]), "wb")
    else:
        id_handle = open("_".join(["/".join([args.parentDir, args.outDirName, "ids", args.sampleId]), "_output.txt"]), "wb")
        
    # use the column names used in R
    id_handle.write("\t".join(["id", "class", "pos", "qual", "aScore", "numN", "readLen", "junction", "strand",
                               "posR2", "qualR2", "aScoreR2", "numNR2", "readLenR2", "junctionR2", "strandR2"]) + "\n")
    
    for j in junctions:
        # print circular read ids
        for elem in junctions[j].circularReads:
            # if we pass a cutoff for read1 and read2 scores, then we are calculating an FDR which we want to keep low
            # otherwise we are calculating a probability of observing these reads under the null that this really is a circle and then we want to get rid of those with low p-values
            if junctions[j].fdr and ((cutoff and float(junctions[j].fdr) > args.reportFDR)
                or (not cutoff and float(junctions[j].fdr) < args.reportFDR)):
                readClass = "circArtifact"
            elif cutoff:
                if int(elem.juncRead.aScore) >= int(cutoff[0]) and (args.singleEnd or int(elem.useMate.aScore) >= int(cutoff[1])):
                    readClass = "circStrong"
                else:
                    readClass = "circFailed"
            else:
                readClass = "circStrong"
                
            myInfo = "\t".join([str(elem.juncRead.offset), str(elem.juncRead.mapQual), str(elem.juncRead.aScore),
                                       str(elem.juncRead.numN), str(elem.juncRead.readLen), str(j), str(elem.juncRead.flag)])
            if args.singleEnd:
                pairInfo = "\t".join(["NA", "NA", "NA", "NA", "NA", "NA", "NA"])
            else:
                pairInfo = "\t".join([str(elem.useMate.offset),
                                       str(elem.useMate.mapQual), str(elem.useMate.aScore), str(elem.useMate.numN), str(elem.useMate.readLen),
                                       str(elem.useMate.refName), str(elem.useMate.flag)])
            
            id_handle.write("\t".join([str(elem.juncRead.name), readClass, myInfo, pairInfo]) + "\n")
        
        # print linear read ids
        for elem in junctions[j].linearReads:
            if junctions[j].fdr and ((cutoff and float(junctions[j].fdr) > args.reportFDR)
                or (not cutoff and float(junctions[j].fdr) < args.reportFDR)):
                readClass = "linearArtifact"
            elif cutoff:
                if int(elem.juncRead.aScore) >= int(globalCutOff[0]) and (args.singleEnd or int(elem.useMate.aScore) >= int(globalCutOff[1])):
                    readClass = "linearStrong"
                else:
                    readClass = "linearFailed"
            else:
                readClass = "linearStrong"
            
            myInfo = "\t".join([str(elem.juncRead.offset), str(elem.juncRead.mapQual), str(elem.juncRead.aScore),
                                       str(elem.juncRead.numN), str(elem.juncRead.readLen), str(j), str(elem.juncRead.flag)])
            if args.singleEnd:
                pairInfo = "\t".join(["NA", "NA", "NA", "NA", "NA", "NA", "NA"])
            else:
                pairInfo = "\t".join([str(elem.useMate.offset),
                                       str(elem.useMate.mapQual), str(elem.useMate.aScore), str(elem.useMate.numN), str(elem.useMate.readLen),
                                       str(elem.useMate.refName), str(elem.useMate.flag)])
            
            id_handle.write("\t".join([str(elem.juncRead.name), readClass, myInfo, pairInfo]) + "\n")
            
        # print multimapped read ids
        for elem in junctions[j].multimappedReads:
            myInfo = "\t".join([str(elem.juncRead.offset), str(elem.juncRead.mapQual), str(elem.juncRead.aScore),
                                       str(elem.juncRead.numN), str(elem.juncRead.readLen), str(j), str(elem.juncRead.flag)])
            if args.singleEnd or not elem.useMate:
                pairInfo = "\t".join(["NA", "NA", "NA", "NA", "NA", "NA", "NA"])
            else:
                pairInfo = "\t".join([str(elem.useMate.offset),
                                       str(elem.useMate.mapQual), str(elem.useMate.aScore), str(elem.useMate.numN), str(elem.useMate.readLen),
                                       str(elem.useMate.refName), str(elem.useMate.flag)])
            
            id_handle.write("\t".join([str(elem.juncRead.name), "multimapped", myInfo, pairInfo]) + "\n")
            
        # print decoy read ids
        if not args.singleEnd:
            for elem in junctions[j].decoyReads:
                myInfo = "\t".join([str(elem.juncRead.offset), str(elem.juncRead.mapQual), str(elem.juncRead.aScore),
                                       str(elem.juncRead.numN), str(elem.juncRead.readLen), str(j), str(elem.juncRead.flag)])

                pairInfo = "\t".join([str(elem.useMate.offset),
                                       str(elem.useMate.mapQual), str(elem.useMate.aScore), str(elem.useMate.numN), str(elem.useMate.readLen),
                                       str(elem.useMate.refName), str(elem.useMate.flag)])
        
                id_handle.write("\t".join([str(elem.juncRead.name), "decoy", myInfo, pairInfo]) + "\n")
                
            # print anomaly read ids
            for elem in junctions[j].anomalyReads:
                myInfo = "\t".join([str(elem.juncRead.offset), str(elem.juncRead.mapQual), str(elem.juncRead.aScore),
                                       str(elem.juncRead.numN), str(elem.juncRead.readLen), str(j), str(elem.juncRead.flag)])

                pairInfo = "\t".join([str(elem.useMate.offset),
                                       str(elem.useMate.mapQual), str(elem.useMate.aScore), str(elem.useMate.numN), str(elem.useMate.readLen),
                                       str(elem.useMate.refName), str(elem.useMate.flag)])
        
                id_handle.write("\t".join([str(elem.juncRead.name), "anomaly", myInfo, pairInfo]) + "\n")
                
            # print unmapped read ids
            for elem in junctions[j].unmappedReads:
                myInfo = "\t".join([str(elem.juncRead.offset), str(elem.juncRead.mapQual), str(elem.juncRead.aScore),
                                       str(elem.juncRead.numN), str(elem.juncRead.readLen), str(j), str(elem.juncRead.flag)])
                pairInfo = "\t".join(["NA", "NA", "NA", "NA", "NA", "NA", "NA"])
                id_handle.write("\t".join([str(elem.juncRead.name), "unmapped", myInfo, pairInfo]) + "\n")
            
    id_handle.close()

# loop through sam file, either creating a juncReadObj & juncObj (if this is read1s)
# or updating mateGenome, mateJunction, or mateRegJunction if this is an alignment file from read2
# to genome, all junctions, or regular junctions respectively.
# The sam file may contain both aligned and unaligned reads. Only aligned reads will actually be stored.
# param readType: "r1rj" is Rd1 to regular junction index, "r1j" is Rd1 to all junction,
#                 "gMate" is Rd2 to genome, "jMate" is Rd2 to all junction, "rjMate" is Rd2 to regular junction
def parseSam(samFile, readType):
    if args.verbose:
        print "samFile:" + samFile
        print "readType:" + readType
        
    handle = open(samFile, "rU")

    for line in handle:
        if not line.startswith("@"): # ignore header lines
            try:
                read = newReadObj(line.strip().split(), args.fastqIdStyle)
                # only need to store info if the read actually aligned 
                if read.aScore:
                    readBase = read.baseName  # part of the read that is the same between Rd1 and Rd2
                    # read1 that didn't map to genome, or ribo and then mapped to reg junctions or to all junctions and not to reg junction
                    # since in unaligned mode we report all reads, we want to only count the first time in the file we see it (output order is highest alignment score first)
                    if ((readType == "r1j" and readBase in nonRegIds and readBase not in juncReads) or
                        (readType == "r1rj" and readBase in regIds and readBase not in juncReads)):
                        curJuncReadObj = juncReadObj(read)  # create object for this read (mate info empty for now) to put in junction and juncReads dicts
                        juncReads[readBase] = curJuncReadObj  # will look up this obj later to populate mates
                            
                        # this is first read aligning to this junction, need to create juncObj first
                        if not read.refName in junctions:
                            curJuncObj = juncObj(read.refName)
                            junctions[read.refName] = curJuncObj
                        # initially just append all reads to unknownReads, we will later remove some and move to circularReads, decoyReads, or unmappedReads
                        junctions[read.refName].unknownReads.append(curJuncReadObj)
                    elif readBase in juncReads: # this is a mate so we only care about it if it's Rd1 was stored
                        if readType == "gMate" and not juncReads[readBase].mateGenomic:  # only if we haven't already found the primary genomic alignment
                            juncReads[readBase].mateGenomic = read
                        else:
                            # this is a junction mate, need to check offset 
                            if (int(read.offset) >= (JUNC_MIDPOINT - int(read.readLen) + args.overhang + 1) and
                                int(read.offset) <= (JUNC_MIDPOINT - args.overhang + 1)):
                                # if it overlaps the junction, add it to the appropriate mate field
                                if readType == "jMate" and not juncReads[readBase].mateJunction:  # only if we haven't already found the primary junction alignment
                                    juncReads[readBase].mateJunction = read
                                elif readType == "rjMate" and not juncReads[readBase].mateRegJunction:  # only if we haven't already found the primary reg junction alignment
                                    juncReads[readBase].mateRegJunction = read
                                elif readType == "dMate" and not juncReads[readBase].mateDenovoJunction:  # only if we haven't already found the primary denovo junction alignment
                                    juncReads[readBase].mateDenovoJunction = read
            except Exception as e:
                print "Exception"
                print e
                print "error:", sys.exc_info()[0]
                print "parsing sam output for", line
                
    handle.close()
    
def updateReads():
    
    # for each junction
    for j in junctions:
        # all reads were initially assigned to unknownReads, need to see which are circular,etc
        while len(junctions[j].unknownReads) > 0:
            r = junctions[j].unknownReads.pop()
            r.updateInfo(junctions[j], args.singleEnd, args.unalignedMode) # figure out which mate is best, whether this means it is decoy or not
            
            # assign read to correct bucket
            if r.readType == "c":
                junctions[j].circularReads.append(r)
            elif r.readType == "d":
                junctions[j].decoyReads.append(r)
            elif r.readType == "l":
                junctions[j].linearReads.append(r)
            elif r.readType == "a":
                junctions[j].anomalyReads.append(r)
            elif r.readType == "i":
                # we are ignoring this read because Rd1 aligned to reg and Rd2 aligned to scrambled
                pass
            else:
                if r.mateGenomic or r.mateJunction or r.mateRegJunction or r.mateDenovoJunction:
                    # the mate did map, but had very poor mapping quality so we are discarding it
                    junctions[j].multimappedReads.append(r)
                else:
                    # either mate was not in the file or did not align to either genome or junction index
                    junctions[j].unmappedReads.append(r)
                    

if __name__  == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--parentDir', help='path to alignment parent directory for this dataset', required=True)
    parser.add_argument('-s', '--sampleId', help='string used to identify this sample', required=True)
    parser.add_argument('-o', '--outDirName', help='name of directory to output files under parentDir, created if does not exist', required=True)
    parser.add_argument('-q', '--fastqIdStyle', help='type of read ids used', required=True, choices=['appended', 'complete'])
    parser.add_argument('-f', '--reportFDR', help='FDR to select circular artifacts in report', type=float, default=0.9)
    parser.add_argument('-e', '--seqErrorRate', help='per base sequencing error rate', type=float, default=0.01)
    parser.add_argument('-a1', '--aScore1', help='specify bowtie min AS for read1 instead of using global FDR cutoff', type=int)
    parser.add_argument('-a2', '--aScore2', help='specify bowtie min AS for read1 instead of using global FDR cutoff', type=int)
    parser.add_argument('-oh', '--overhang', help='how much you want the read to overlap a junction to be considered aligned to junction', type=int, required=True)
    parser.add_argument('-j', '--junctionIdDirSuffix', help='suffix appended to junc and reg to find overlapped reads for this run', default='')
    parser.add_argument('-se', '--singleEnd', help='is this single end read data', action='store_true')
    parser.add_argument('-u', '--unalignedMode', help='is this an unaligned mode run', action='store_true')
    parser.add_argument('-v', '--verbose', help='print extra debugging info', action='store_true')
    args = parser.parse_args()
    
    if args.verbose:
        print "parentDir:", args.parentDir
        print "sampleId:", args.sampleId
        print "outDirName:", args.outDirName
        print "style:", args.fastqIdStyle
        print "a1:", args.aScore1
        print "a2:", args.aScore2
        print "overhang:", args.overhang
        print "id dir suffix:", args.junctionIdDirSuffix
        print "unaligned:", args.unalignedMode
    
    # make output dirs if they don't exist
    utils_os.createDirectory("/".join([args.parentDir, args.outDirName]))
    utils_os.createDirectory("/".join([args.parentDir, args.outDirName, "reports"]))  # these are the reports using the naive method
    utils_os.createDirectory("/".join([args.parentDir, args.outDirName, "glmReports"])) # GLM will be run later and those reports will be stored here
    utils_os.createDirectory("/".join([args.parentDir, args.outDirName, "glmModels"]))  # GLM will be run later and those models will be stored here
    utils_os.createDirectory("/".join([args.parentDir, args.outDirName, "ids"]))  # txt files of read ids assigned to circular or linear category
    # just doing read1s 
    if args.sampleId.endswith("1"):
    
        try:
            # populate the ignoreIds, regIds, nonRegIds for this file and also print out regIds and nonRegIds to juncNonGR file
            ignoreIds = {} # ribo and genome aligned 
            regIds = {} # regular junction overlapped and not ribo or genome aligned (aligned to reg-only index)
            nonRegIds = {} # junction overlapped and not ribo or genome aligned or regular-junction aligned
            
            # we treat denovo reads as the equivalent of the junction reads in unalignedMode
            if args.unalignedMode:
                useJuncStr = "denovo"
            else:
                useJuncStr = "junction"
                
            selectCandidateIds()
            
            if args.verbose:
                print "ignoreIds (aligned to genome):", str(len(ignoreIds))
                print "linearIds:", str(len(regIds))
                print "scrambledIds:", str(len(nonRegIds))
                
            juncReads = {}  # base read id: juncReadObj
            junctions = {}  # junction id: list of juncObjs mapping to this junction
                
            if not args.unalignedMode:
                # make a pass through the sam file for read 1 to regular junctions to populate dict with juncReadObjs (key is portion of read id shared by Rd1 and Rd2)
                #     and at the same time create entry in junction dictionary and add this read obj to the junction list of reads
                parseSam("".join(["/".join([args.parentDir, "orig", "reg", args.sampleId]), "_reg_output.sam"]), "r1rj")
            
            # make a pass through the sam for read 1 to all junctions file to populate dict with juncReadObjs (key is portion of read id shared by Rd1 and Rd2)
            #     and at the same time create entry in junction dictionary and add this read obj to the junction list of reads
            parseSam("_".join(["/".join([args.parentDir, "orig", useJuncStr, args.sampleId]), useJuncStr, "output.sam"]), "r1j")
            
            # get mate alignment data if this is from paired end sequencing
            if not args.singleEnd:
                if args.unalignedMode:
                    alignedMateId = args.sampleId[10:-1] + "2"  # trim off unaligned_ from start of id, change 1 to 2
                else:
                    alignedMateId = args.sampleId[:-1] + "2"
                
                # make a pass through Rd2 to genome sam file to update mateGenomic in each juncReadObj
                parseSam("".join(["/".join([args.parentDir, "orig", "genome", alignedMateId]), "_genome_output.sam"]), "gMate")
                
                # make a pass through Rd2 to junction sam file to update mateRegJunction in each juncReadObj
                parseSam("".join(["/".join([args.parentDir, "orig", "reg", alignedMateId]), "_reg_output.sam"]), "rjMate")
                
                # make a pass through Rd2 to junction sam file to update mateJunction in each juncReadObj
                parseSam("".join(["/".join([args.parentDir, "orig", "junction", alignedMateId]), "_junction_output.sam"]), "jMate")
                
                # make a pass through Rd2 to de novo sam file to update mateDenovo in each juncReadObj
                if args.unalignedMode:
                    parseSam("".join(["/".join([args.parentDir, "orig", "denovo", args.sampleId[:-1]]), "2_denovo_output.sam"]), "dMate")
                
            if args.verbose:
                print "sample id:", str(args.sampleId)
                print "single end?", str(args.singleEnd)
                print "num junctions with aligned reads:", str(len(junctions))
                print "num reads:", str(len(juncReads))

            updateReads()  # assign reads to categories based on alignment data (circular, linear, decoy, etc)
            
            # print out some alignment statistics after reads have all been assigned to a category
            if args.verbose:
                numCirc = numDecoy = numLinear = numAnomaly = numUnmapped = numMultimapped = numUnknown = 0
                for j in junctions:
                    numCirc = numCirc + len(junctions[j].circularReads)
                    numDecoy = numDecoy + len(junctions[j].decoyReads)
                    numLinear = numLinear + len(junctions[j].linearReads)
                    numAnomaly = numAnomaly + len(junctions[j].anomalyReads)
                    numUnmapped = numUnmapped + len(junctions[j].unmappedReads)
                    numMultimapped = numMultimapped + len(junctions[j].multimappedReads)
                    numUnknown = numUnknown + len(junctions[j].unknownReads)
                print "number of reads kept in each category after update (multimapped & unknown should be 0, anomaly & decoy should be 0 for SE data):"
                print "circ: ", str(numCirc), ", decoy: ", str(numDecoy), ", linear: ", str(numLinear), ", anomaly: ", str(numAnomaly), ", unmapped: ", str(numUnmapped), ", multimapped: ", str(numMultimapped), ", unknown: ", str(numUnknown)
                print str(sum([numCirc, numDecoy, numLinear, numAnomaly, numUnmapped, numMultimapped, numUnknown]))
            
            if args.aScore1 and args.aScore2:
                globalCutOff = (args.aScore1, args.aScore2)
                globalDecoyMMrate = None
                if args.verbose:
                    print "globalCutOff specified:", globalCutOff
            else:
                # use decoy distribution
                globalCutOff = None
                globalDecoyMMrate = getDecoyMismatchRate()
                
                # if there were no decoys, use the default seqErrorRate
                if globalDecoyMMrate == None:
                    globalDecoyMMrate = args.seqErrorRate
                
                if args.verbose:
                    print "using decoy rate:", globalDecoyMMrate
                    
            
            reportCircularReads(globalCutOff)  # output reports
            reportAllReadIds2(globalCutOff)   # output ids used in the reports for further manual investigation if desired + GLM uses these category assignments
            
        except Exception as e:
            print "Exception"
            print e
        