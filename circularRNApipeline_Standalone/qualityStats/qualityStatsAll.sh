#!/bin/bash

# Sets up required output directory structure for storing unaligned reads and alignment statistics reports.
# This wrapper manages calling qualityStatsSingleSample on up to 6 samples at a time (to prevent overloading
# computer CPU) and consolidating all individual report files into reports for the entire dataset
# under sampleStats.

PAR_DIR=$1
DATASET_NAME=$2
if [ $# -ge 3 ]
then
  REPORTDIR_NAME=$3
else
  REPORTDIR_NAME=circReads
fi

NTRIM=$4
JUNCTION_ID_SUFFIX=$5

TASK_DATA_FILE="${PAR_DIR}/taskIdFiles/${DATASET_NAME}.txt"
NUM_FILES=`cat $TASK_DATA_FILE | wc -l`

# make the output directories
OUT_DIR=${PAR_DIR}/${DATASET_NAME}/sampleStats
mkdir -p $OUT_DIR
mkdir -p ${PAR_DIR}/${DATASET_NAME}/orig/unaligned
mkdir -p ${PAR_DIR}/${DATASET_NAME}/orig/unaligned/forDenovoIndex

# for each sample, generate sample stats in the output directory
for (( i=1; i<=NUM_FILES; i++ ))
do
  READ_FILE=`awk 'FNR == '${i}' {print $1}' $TASK_DATA_FILE`
  SAMPLE_ID=`awk 'FNR == '${i}' {print $2}' $TASK_DATA_FILE`
  READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2
  

  sh qualityStatsSingleSample.sh $READ_FILE $SAMPLE_ID $PAR_DIR $DATASET_NAME $REPORTDIR_NAME ${OUT_DIR} ${NTRIM} ${JUNCTION_ID_SUFFIX} &
  echo "Launched qualityStatsSingleSample.sh into the background "`date`
  run_count=`ps -ealf | grep qualityStatsSingleSample.sh | grep ${USER} | grep -v grep | wc -l`
  while [ "$run_count" -gt 5 ]
  do
    sleep 1
    run_count=`ps -ealf | grep qualityStatsSingleSample.sh | grep ${USER} | grep -v grep | wc -l`
  done
done

# Wait until finished
echo 'Done looping over qualityStatsSingleSample, awaiting the last few qualityStatsSingleSample.sh to finish '`date`
while [ "$run_count" -gt 0 ]
do
  sleep 1
  run_count=`ps -ealf | grep qualityStatsSingleSample.sh | grep ${USER} | grep -v grep | wc -l`
done

# cat all of those files then delete the original separate files
sh qualityStatsCat.sh ${OUT_DIR}

