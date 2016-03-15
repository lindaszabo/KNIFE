#!/bin/sh

## This is just a wrapper that passes on the call to predictJunctions_tableData.r to run the
# GLM analysis. It is not called in unaligned mode.

ALIGN_PARDIR=$1
DATASET_NAME=$2
MODE=$3
REPORTDIR_NAME=$4

MODEL_OUTDIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/glmModels
REPORT_OUTDIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/glmReports

# get info about samples
TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt
NUM_FILES=`cat $TASK_DATA_FILE | wc -l`

for (( i=1; i<=NUM_FILES; i++ ))
do
  SAMPLE_ID=`awk 'FNR == '${i}' {print $2}' $TASK_DATA_FILE`
  READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2
  
  # only looking at read1s for now
  if [ "$READ_NUM" -eq 1 ]
  then
    CLASS_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}_${APPENDED}_output.txt
    MODEL_OUT=${MODEL_OUTDIR}/${SAMPLE_ID}_${APPENDED}_glm.RData
    LINEAR_JUNC_OUT=${REPORT_OUTDIR}/${SAMPLE_ID}_${APPENDED}_linearJuncProbs.txt
    CIRC_JUNC_OUT=${REPORT_OUTDIR}/${SAMPLE_ID}_${APPENDED}_circJuncProbs.txt
    echo "./predictJunctions_tableData.r ${CLASS_FILE} ${MODEL_OUT} ${LINEAR_JUNC_OUT} ${CIRC_JUNC_OUT}"
    ./predictJunctions_tableData.r ${CLASS_FILE} ${MODEL_OUT} ${LINEAR_JUNC_OUT} ${CIRC_JUNC_OUT}
  fi
done
