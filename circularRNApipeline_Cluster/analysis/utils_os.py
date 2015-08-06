import os
import errno
import re

# create directory if it does not exist
def createDirectory(dirPath):
    try:
        if not os.path.isdir(dirPath):
            os.makedirs(dirPath)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
        
# parse a string of interest out of a filename
def getFileId(patt, grp, filename):
    match = patt.search(filename)
    if match:
        return match.group(grp)
    else:
        print "filename not structured correctly: " + str(filename)

# good way to grep, particularly for large files
def grep(pattern, file_obj, include_line_nums=False):
    grepper = re.compile(pattern)
    for line_num, line in enumerate(file_obj):
        if grepper.search(line):
            if include_line_nums:
                yield (line_num, line)
            else:
                yield line


