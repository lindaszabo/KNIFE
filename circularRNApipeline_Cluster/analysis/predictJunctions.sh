#!/bin/sh

## This is just a wrapper that passes on the call to predictJunctions_tableData.r to run the
# GLM analysis. It is not called in unaligned mode.

CLUSTER_TYPE=$1
ALIGN_PARDIR=$2
DATASET_NAME=$3
MODE=$4
REPORTDIR_NAME=$5

MODEL_OUTDIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/glmModels
REPORT_OUTDIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/glmReports

# get info about samples
TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt
NUM_FILES=`cat $TASK_DATA_FILE | wc -l`

source ../sampleInfo.sh ${CLUSTER_TYPE} # get sample-specific variables from TASK_DATA_FILE
source ../depends.sh ${CLUSTER_TYPE}  # load R
ALIGNED_SAMPLE_ID=`echo ${SAMPLE_ID} | sed 's/unaligned_//'` # if we are in unaligned mode, we will need to work with both aligned and unaligned data
READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2
  
# only looking at read1s for now
if [ "$READ_NUM" -eq 1 ]
then 
  CLASS_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}__output.txt
  MODEL_OUT=${MODEL_OUTDIR}/${SAMPLE_ID}_${APPENDED}_glm.RData
  LINEAR_JUNC_OUT=${REPORT_OUTDIR}/${SAMPLE_ID}_${APPENDED}_linearJuncProbs.txt
  CIRC_JUNC_OUT=${REPORT_OUTDIR}/${SAMPLE_ID}_${APPENDED}_circJuncProbs.txt
    
  ./predictJunctions_tableData.r ${CLASS_FILE} ${MODEL_OUT} ${LINEAR_JUNC_OUT} ${CIRC_JUNC_OUT}
fi
