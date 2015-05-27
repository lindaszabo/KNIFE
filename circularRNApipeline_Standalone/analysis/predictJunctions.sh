#!/bin/sh

ALIGN_PARDIR=$1
DATASET_NAME=$2
MODE=$3
REPORTDIR_NAME=$4

MODEL_OUTDIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/glmModels
REPORT_OUTDIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/glmReports

# get info about samples
if [[ $MODE = *unaligned* ]]
then
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}_unaligned.txt
  APPENDED=denovo
else
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt
fi
NUM_FILES=`cat $TASK_DATA_FILE | wc -l`

for (( i=1; i<=NUM_FILES; i++ ))
do
  SAMPLE_ID=`awk 'FNR == '${i}' {print $2}' $TASK_DATA_FILE`
  READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2
  
  # only looking at read1s for now
  if [ "$READ_NUM" -eq 1 ]
  then
    LINEAR_SCORE_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/orig/temp/${SAMPLE_ID}_reg_output.txt
    CIRC_SCORE_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/orig/temp/${SAMPLE_ID}_junction_output.txt
    CLASS_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}_${APPENDED}_output.txt
    MODEL_OUT=${MODEL_OUTDIR}/${SAMPLE_ID}_${APPENDED}_glm.RData
    LINEAR_JUNC_OUT=${REPORT_OUTDIR}/${SAMPLE_ID}_${APPENDED}_linearJuncProbs.txt
    CIRC_JUNC_OUT=${REPORT_OUTDIR}/${SAMPLE_ID}_${APPENDED}_circJuncProbs.txt
    echo "./predictJunctions_tableData.r ${LINEAR_SCORE_FILE} ${CIRC_SCORE_FILE} ${CLASS_FILE} ${MODEL_OUT} ${LINEAR_JUNC_OUT} ${CIRC_JUNC_OUT}"
    ./predictJunctions_tableData.r ${LINEAR_SCORE_FILE} ${CIRC_SCORE_FILE} ${CLASS_FILE} ${MODEL_OUT} ${LINEAR_JUNC_OUT} ${CIRC_JUNC_OUT}
  fi
done
