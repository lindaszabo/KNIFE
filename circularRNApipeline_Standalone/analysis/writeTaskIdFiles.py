# Generate a tab-delimited file containing 1 line per fastq file with path to file and
# the sample name identifying that file. This is used by the task array jobs for all alignment
# and analysis steps. Also create all output directories for alignment files now.

# usage: python writeTaskIdFiles.py -r /srv/gsfs0/projects/salzman/Ovarian_test_libraries_2013
#                                   -a /srv/gsfs0/projects/salzman/Linda/alignments
#                                   -d OvarianArrayJob

import argparse
import os
import utils_os
import re

if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--readDir', help='directory containing fastq files', required=True)
    parser.add_argument('-a', '--alignmentParDir', help='parent directory where directory for alignments for this data set will be created', default="/srv/gsfs0/projects/salzman/Linda/alignments")
    parser.add_argument('-t', '--taskDir', help='name of directory under alignmentParDir to output task file', default="taskIdFiles")
    parser.add_argument('-d', '--dataSet', required=True,
                        help='name of directory under alignmentParDir that all alignment files will be written to')
    parser.add_argument('-v', '--verbose', help='print extra debugging info', action='store_true')
    parser.add_argument('-u', '--unalignedMode', help='pass this flag if we are aligning the unaligned reads from previous pipeline run', action='store_true')
    args = parser.parse_args()
    
    try:
        ### create output directories if they don't exist, and the alignment subdirectories
        utils_os.createDirectory("/".join([args.alignmentParDir, args.taskDir]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet]))

        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "denovo"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "still_unaligned"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "genome"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "junction"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ribo"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "reg"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ids"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ids", "denovo"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ids", "genome"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ids", "junction"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ids", "ribo"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ids", "reg"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ids", "juncNonGR"]))
        utils_os.createDirectory("/".join([args.alignmentParDir, args.dataSet, "orig", "ids", "denovoNonGR"]))
        #### open file for output ####
        if args.unalignedMode:
            # we are creating a list of the unaligned fastq file names
            out_handle = open("/".join([args.alignmentParDir, args.taskDir, args.dataSet]) + "_unaligned.txt", "wb")
        else:
            out_handle = open("/".join([args.alignmentParDir, args.taskDir, args.dataSet]) + ".txt", "wb")
        
        # we identify reads of interest as fastq files that have _R1, _R2, _1, or _2 in the name.
        # anything after the R1 is discarded in the file name for all subsequent output.
        patt_file = re.compile("(.+_R?[1|2]).*?\.f(ast)*q(\.gz)*$")
        
        # for each file, write a line with its required task job info
        for f in os.listdir(args.readDir):
            match = patt_file.search(f)
            if match:
                out_handle.write("/".join([args.readDir, f]))  # full path to fastq file for bowtie
                out_handle.write("\t")
                out_handle.write(match.group(1))  # file id we will be using
                out_handle.write("\n")
            
        out_handle.close()
    except Exception as e:
        print "Exception: ", e