import argparse
import shutil
import os
import utils_os

if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--alignmentParDir', help='same alignmentParDir used in original run of findCircularRNA', required=True)
    parser.add_argument('-d', '--dataSet', required=True,
                        help='name of directory that original files were output to under alignmentParDir')
    parser.add_argument('-t', '--taskDir', help='name of directory under alignmentParDir to output task file', default="taskIdFiles")
    parser.add_argument('-m', '--mode', help='use swap if you want to move into the swapped directory, restore if you want to move back to original location',
                        choices=['swap', 'restore'], default="swap")
    parser.add_argument('-v', '--verbose', help='print extra debugging info', action='store_true')
    args = parser.parse_args()
    
    try:
        swappedDataSet = args.dataSet+"Swapped"
        origSam = "/".join([args.alignmentParDir, args.dataSet, "orig"])
        swappedSam = "/".join([args.alignmentParDir, swappedDataSet, "orig"])
        
        if args.mode == "swap":
            ### create output directories if they don't exist, and the alignment subdirectories
            ### these are all that are created by writeTaskIdFiles.py for original run, all other directories are output as usual during the run in analysis mode
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet]))
    
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "denovo"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "still_unaligned"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "genome"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "junction"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ribo"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "reg"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ids"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ids", "denovo"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ids", "genome"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ids", "junction"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ids", "ribo"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ids", "reg"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ids", "juncNonGR"]))
            utils_os.createDirectory("/".join([args.alignmentParDir, swappedDataSet, "orig", "ids", "denovoNonGR"]))
        
            ### just copy the taskDir file, it is going to have the sample ids for the original order files but I don't think that causes a problem since we aren't using the original read files anymore
            shutil.copyfile("/".join([args.alignmentParDir, args.taskDir, args.dataSet]) + ".txt", "/".join([args.alignmentParDir, args.taskDir, swappedDataSet]) + ".txt")
            
            fromSamDir = origSam
            toSamDir = swappedSam
        else:
            fromSamDir = swappedSam
            toSamDir = origSam
            
        ### mv the sam files, renaming them in the process
        for subdir in next(os.walk(fromSamDir))[1]:  # just get the top level subdirectories
            if subdir in ["genome", "junction", "reg", "ribo", "denovo"]:  # if it is one of the alignment output directories
                for f in os.listdir("/".join([fromSamDir, subdir])):
                    newf = f  # avoid exceptions if a file in the directory doesn't contain a read identifier
                    if "_R1" in f:
                        newf = f.replace("_R1", "_R2")
                    elif "_1" in f:
                        newf = f.replace("_1", "_2")
                    elif "_R2" in f:
                        newf = f.replace("_R2", "_R1")
                    elif "_2" in f:
                        newf = f.replace("_2", "_1")
                    
                    #print "/".join([fromSamDir, subdir, f]), "is moving to", "/".join([toSamDir, subdir, newf])
                    shutil.move("/".join([fromSamDir, subdir, f]), "/".join([toSamDir, subdir, newf]))
    except Exception as e:
        print "Exception: ", e