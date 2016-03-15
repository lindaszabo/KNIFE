#!/bin/sh

# This is just a wrapper that passes on the call to combineSwappedReadsNaive.py 

# store commandline params and look up job info in txt file
CLUSTER_TYPE=$1
DIR_A=$2
DIR_B=$3
READ_TYPE=$4

source ../depends.sh ${CLUSTER_TYPE} # load python module if necessary
 
python combineSwappedReadsNaive.py -a ${DIR_A} -b ${DIR_B} -q ${READ_TYPE} -v

  
