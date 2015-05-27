#!/bin/sh

# prints out a file that will be used by R code to train GLM
#  read id, position of alignment, mapping quality, alignment score, number of Ns, length, junction id
ALIGN_PARDIR=$1
DATASET_NAME=$2
MODE=$3

TEMP_OUT_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/orig/temp # to store the output files with data parsed from SAM
mkdir -p ${TEMP_OUT_DIR}

# get info about samples
if [[ $MODE = *unaligned* ]]
then
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}_unaligned.txt
  JUNC_TYPE=denovo  # where the sam files are under orig
else
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt
  JUNC_TYPE=junction
fi
NUM_FILES=`cat $TASK_DATA_FILE | wc -l`

for (( i=1; i<=NUM_FILES; i++ ))
do
  SAMPLE_ID=`awk 'FNR == '${i}' {print $2}' $TASK_DATA_FILE`
  ALIGNED_SAMPLE_ID=`echo ${SAMPLE_ID} | sed 's/unaligned_//'` # if we are in unaligned mode, we will need to work with both aligned and unaligned data
  READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2
  echo "READ_NUM ${READ_NUM}"
  # only looking at read1s for now
  if [ "$READ_NUM" -eq 1 ]
  then
  
    JUNC_SAM=${ALIGN_PARDIR}/${DATASET_NAME}/orig/${JUNC_TYPE}/${SAMPLE_ID}_${JUNC_TYPE}_output.sam
    JUNC_OUT=${TEMP_OUT_DIR}/${SAMPLE_ID}_junction_output.txt
  
    REG_SAM=${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${ALIGNED_SAMPLE_ID}_reg_output.sam
    REG_OUT=${TEMP_OUT_DIR}/${SAMPLE_ID}_reg_output.txt
    
    echo "JUNC_SAM ${JUNC_SAM}"
    # parse junctions sam file. remove first 2 comment lines, parse those with a next best reported, 
    sed '1,2d' ${JUNC_SAM} | grep "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $14 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' > ${JUNC_OUT}
    sed '1,2d' ${JUNC_SAM} | grep -v "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $13 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' >> ${JUNC_OUT}
    
    echo "REG_SAM ${REG_SAM}"
    # parse reg sam file
    sed '1,2d' ${REG_SAM} | grep "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $14 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' > ${REG_OUT}
    sed '1,2d' ${REG_SAM} | grep -v "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $13 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' >> ${REG_OUT}
  fi
done



