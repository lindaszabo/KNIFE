#!/bin/sh

## This is just a wrapper that passes on the call to predictJunctions_tableData.r to run the
# GLM analysis. While there is code here for running the GLM on de novo candidates (mode unaligned),
# the GLM did not work well for the de novo candidates. It is possible using the model developed
# during the annotation-dependent run could be used to predict on the de novo candidates, but
# we just used the naive method for the paper instead of exploring this further. So this code does
# not actually get called by findCircularRNA.sh wrapper in unaligned mode.

module load R  ## SPECIFIC TO OUR CLUSTER, need to have R loaded into the computing namespace

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

SAMPLE_ID=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $2}' "${TASK_DATA_FILE}"`
ALIGNED_SAMPLE_ID=`echo ${SAMPLE_ID} | sed 's/unaligned_//'` # if we are in unaligned mode, we will need to work with both aligned and unaligned data
READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2
  
# only looking at read1s for now
if [ "$READ_NUM" -eq 1 ]
then 
  if [[ $MODE = *unaligned* ]]
  then
    # need to combine the linear ids from the original run with the denovo run since linears come from original and circ/dup come from denovo
    grep linear ${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${ALIGNED_SAMPLE_ID}__output.txt > ${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}_output_linear.txt  # get the original linear ids 
    cat ${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}_${APPENDED}_output.txt ${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}_output_linear.txt > ${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}_output_combined.txt  # add those linear to all de novo junctions
    CLASS_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}_output_combined.txt
    rm ${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}_output_linear.txt  # remove the temp file created for catting linear and junction results
  else
    CLASS_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}__output.txt
  fi
  LINEAR_SCORE_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/orig/temp/${ALIGNED_SAMPLE_ID}_reg_output.txt
  CIRC_SCORE_FILE=${ALIGN_PARDIR}/${DATASET_NAME}/orig/temp/${SAMPLE_ID}_junction_output.txt
  MODEL_OUT=${MODEL_OUTDIR}/${SAMPLE_ID}_${APPENDED}_glm.RData
  LINEAR_JUNC_OUT=${REPORT_OUTDIR}/${SAMPLE_ID}_${APPENDED}_linearJuncProbs.txt
  CIRC_JUNC_OUT=${REPORT_OUTDIR}/${SAMPLE_ID}_${APPENDED}_circJuncProbs.txt
    
  ../../../analysis/predictJunctions_tableData.r ${LINEAR_SCORE_FILE} ${CIRC_SCORE_FILE} ${CLASS_FILE} ${MODEL_OUT} ${LINEAR_JUNC_OUT} ${CIRC_JUNC_OUT}
fi
