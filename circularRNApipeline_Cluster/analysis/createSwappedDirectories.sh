#!/bin/sh

# This is just a wrapper that passes on the call to createSwappedDirectories.py 

# store commandline params and look up job info in txt file
CLUSTER_TYPE=$1
ALIGN_PARDIR=$2
DATASET_NAME=$3
MODE=$4

source ../depends.sh ${CLUSTER_TYPE} # load python module if necessary

python createSwappedDirectories.py -a ${ALIGN_PARDIR} -d ${DATASET_NAME} -m ${MODE} -v
 

  
  
