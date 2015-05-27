#!/bin/sh

# ./completeRun.sh /srv/gsfs0/projects/salzman/Linda/DEBUG_MISEQ/trimmed complete /srv/gsfs0/projects/salzman/Linda/alignments NUGENpipeline 10 circReads 55

## parse parameters because we need some for each of the steps. 
READ_DIR=${1}
READ_STYLE=${2}
ALIGN_PARDIR=${3}
DATASET_NAME=${4}
OVERLAP=${5}
if [ $# -ge 6 ]
then
  MODE=${6}
else
  MODE=sam
fi

if [ $# -ge 7 ]
then
  REPORTDIR_NAME=${7}
else
  REPORTDIR_NAME=circReads
fi

# ntrim for denovo
if [ $# -ge 8 ]
then
  NTRIM=${8}
else
  NTRIM=50
fi

# should denovo contains only circles (1 means only circles)
if [ $# -ge 9 ]
then
  DENOVOCIRC=${9}
else
  DENOVOCIRC=1
fi

JUNCTION_DIR_SUFFIX=${10}
RD1_THRESH=${11}
RD2_THRESH=${12}

JUNCTION_MIDPOINT=${13}
## end parse parameters

# initial alignment
echo "./findCircularRNA.sh \"$@\""
./findCircularRNA.sh "$@"  

# run GLM to output reports  
if [[ $MODE != *skipGLM* ]]
then
  cd analysis
  echo "./parseForAnalysis.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE}"
  ./parseForAnalysis.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} # get info we need from sam file for linear junctions and circular junctions
  echo "./predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${REPORTDIR_NAME}"
  ./predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${REPORTDIR_NAME}
  cd ..
fi

cd qualityStats
# generate samples alignment statistics, also generates the unaligned reads

echo "./qualityStatsAll.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${REPORTDIR_NAME} ${NTRIM} ${JUNCTION_DIR_SUFFIX}"
./qualityStatsAll.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${REPORTDIR_NAME} ${NTRIM} ${JUNCTION_DIR_SUFFIX}
cd ..

 
# run Julia's analysis to generate new fasta file to use for unaligned run below
if [[ $MODE != *skipDenovo* ]]
then
  
  cd denovo_scripts
  echo "perl process_directory_unaligned.pl ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${NTRIM} ${DENOVOCIRC}"
  perl process_directory_unaligned.pl ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${NTRIM} ${DENOVOCIRC}
  cd ..
  
  # run unaligned mode
  # not passing junction midpoint here because we're not changing it at this point anyway and since it has been 
  echo "./findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_unaligned ${REPORTDIR_NAME} ${NTRIM} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}"
  ./findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_unaligned ${REPORTDIR_NAME} ${NTRIM} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}
  
fi



