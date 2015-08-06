#!/bin/sh

## Parse out the relevant fields from the sam file for the alignments we are interested in
# to avoid having to do file parsing in R. This creates text files in orig/temp output
# directory that are used by predictJunctions_tableData.r for doing GLM analysis.
CLUSTER_TYPE=$1
ALIGN_PARDIR=$2
DATASET_NAME=$3
MODE=$4

TEMP_OUT_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/orig/temp  # for output files

# get info about samples
if [[ $MODE = *unaligned* ]]
then
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}_unaligned.txt
  JUNC_TYPE=denovo  # where the sam files are under orig
else
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt
  JUNC_TYPE=junction
fi


source ./sampleInfo.sh ${CLUSTER_TYPE}  # get sample-specific variables from TASK_DATA_FILE
ALIGNED_SAMPLE_ID=`echo ${SAMPLE_ID} | sed 's/unaligned_//'` # if we are in unaligned mode, we will need to work with both aligned and unaligned data
READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2


# only looking at read1s for now
if [ "$READ_NUM" -eq 1 ]
then

  JUNC_SAM=${ALIGN_PARDIR}/${DATASET_NAME}/orig/${JUNC_TYPE}/${SAMPLE_ID}_${JUNC_TYPE}_output.sam
  JUNC_OUT=${TEMP_OUT_DIR}/${SAMPLE_ID}_junction_output.txt
  
  # in unaligned mode there are secondary alignments in the sam file, but we only want to consider the primary alignments
  if [[ $MODE = *unaligned* ]]
  then
    head -n 2 ${JUNC_SAM} > ${JUNC_SAM}.primary  # adding in the 2 comment lines so we can still use the sed '1,2d' needed for original file
    awk '$2 == 0 || $2 == 16' ${JUNC_SAM} >> ${JUNC_SAM}.primary
    USE_JUNC_SAM=${JUNC_SAM}.primary
    NUM_PRIMARY_JUNC=`wc -l ${USE_JUNC_SAM}`
  else
    USE_JUNC_SAM=${JUNC_SAM}
  fi

  REG_SAM=${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${ALIGNED_SAMPLE_ID}_reg_output.sam
  REG_OUT=${TEMP_OUT_DIR}/${ALIGNED_SAMPLE_ID}_reg_output.txt

  # parse junctions sam file. remove first 2 comment lines, parse those with a next best reported,
  sed '1,2d' ${USE_JUNC_SAM} | grep "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $14 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' > ${JUNC_OUT}
  sed '1,2d' ${USE_JUNC_SAM} | grep -v "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $13 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' >> ${JUNC_OUT}
  
  if [[ $MODE = *unaligned* ]]
  then
    rm ${USE_JUNC_SAM}
  fi

  # parse reg sam file
  sed '1,2d' ${REG_SAM} | grep "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $14 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' > ${REG_OUT}
  sed '1,2d' ${REG_SAM} | grep -v "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $13 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' >> ${REG_OUT}
fi


