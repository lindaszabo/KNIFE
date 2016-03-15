#!/bin/sh

# This is just a wrapper that passes on the call to filterFDR.py after converting bam files for
# the current sample to sam if necessary, and then removes the temporary sam files after completion.

# store commandline params and look up job info in txt file
CLUSTER_TYPE=$1
MODE=$2
ALIGN_PARDIR=$3
DATASET_NAME=$4
OUTDIR_NAME=$5
READ_TYPE=$6
OVERLAP=$7

if [[ $MODE = *unaligned* ]]
then
  TASK_DATA_FILE="${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}_unaligned.txt"
else
  TASK_DATA_FILE="${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt"
fi

source ../depends.sh ${CLUSTER_TYPE} # load python module if necessary
source ../sampleInfo.sh ${CLUSTER_TYPE} # get sample-specific variables from TASK_DATA_FILE
SAMPLE_MATE_ID=`echo ${SAMPLE_ID%?}2`  # will be full sample id with trailing 1 replaced by 2
READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2

ALIGNED_SAMPLE_ID=`echo ${SAMPLE_ID} | sed 's/unaligned_//'` # if we are in unaligned mode, we will need to work with both aligned and unaligned data
ALIGNED_SAMPLE_MATE_ID=`echo ${SAMPLE_MATE_ID} | sed 's/unaligned_//'`

# Rd1-centric analysis 
if [ "$READ_NUM" -eq 1 ]
then
  
  # convert sample bam files to sam if necessary. just checking for a single sam file, not all, assume if 1 exists they all will
  f=`find ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction -type f -name ${ALIGNED_SAMPLE_ID}_junction_output.sam`

  # if paired end, we will need the mate too. If mate doesn't exist, assume we are working with single-end reads
  f2=`find ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction -type f -name ${ALIGNED_SAMPLE_MATE_ID}_junction_output.*`

  if [ ! -f "$f" ]
  then
    DO_CONVERT=true
    
    # if read1 file was stored as bam, assume read2 is too if it exists
    if [ -f "$f2" ]
    then
      DO_CONVERT_2=true
    fi
  fi

  if [[ $MODE = *unaligned* ]]
  then
    if [ "$DO_CONVERT" = true ]
    then
      # need Rd1 unaligned bam file converted to sam for analysis
      samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/denovo/${SAMPLE_ID}_denovo_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/denovo/${SAMPLE_ID}_denovo_output.sam
    fi

    if [ "$DO_CONVERT_2" = true ]
    then
      # need Rd2 too
      samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/denovo/${SAMPLE_MATE_ID}_denovo_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/denovo/${SAMPLE_MATE_ID}_denovo_output.sam
      
      # need Rd2 genome, reg, and junction
      samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/genome/${ALIGNED_SAMPLE_MATE_ID}_genome_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/genome/${ALIGNED_SAMPLE_MATE_ID}_genome_output.sam
      samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${ALIGNED_SAMPLE_MATE_ID}_junction_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${ALIGNED_SAMPLE_MATE_ID}_junction_output.sam
      samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${ALIGNED_SAMPLE_MATE_ID}_reg_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${ALIGNED_SAMPLE_MATE_ID}_reg_output.sam
    fi
  else
  
    if [ "$DO_CONVERT" = true ]
    then
      # need Rd1 regular and junction bam files converted to sam for analysis 
      samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${SAMPLE_ID}_junction_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${SAMPLE_ID}_junction_output.sam
      samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${SAMPLE_ID}_reg_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${SAMPLE_ID}_reg_output.sam
    
      if [ "$DO_CONVERT_2" = true ]
      then
        # need Rd2 genome, reg, and junction
        samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/genome/${SAMPLE_MATE_ID}_genome_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/genome/${SAMPLE_MATE_ID}_genome_output.sam
        samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${SAMPLE_MATE_ID}_junction_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${SAMPLE_MATE_ID}_junction_output.sam
        samtools view ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${SAMPLE_MATE_ID}_reg_output.bam > ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${SAMPLE_MATE_ID}_reg_output.sam
      fi
    fi
  fi
    
  # set up alternate python arguments based on parameters passed
  if [ $# -eq 8 ]
  then
    OPT_ARGS=`echo "-j ${8}"`  # id suffix was passed
  elif [ $# -eq 9 ]
  then
    OPT_ARGS=`echo "-a1 ${8} -a2 ${9}"`
  elif [ $# -eq 10 ]
  then
    OPT_ARGS=`echo "-j ${8} -a1 ${9} -a2 ${10}"`
  fi
  
  # add param for single reads if we didn't find a read2 file
  if [ ! -f "$f2" ]
  then
    OPT_ARGS=`echo ${OPT_ARGS} -se`  
  fi
  
  # add param for unaligned mode
  # to keep consistent as much as possible with StandAlone version of code,
  # search for unalign instead of unaligned to handle Standalone code where R2analysis_unalign indicates run in unaligned mode on pre-existing files
  if [[ $MODE = *unalign* ]]
  then
    OPT_ARGS=`echo ${OPT_ARGS} -u`
  fi
  
  python filterFDR.py -p ${ALIGN_PARDIR}/${DATASET_NAME} -s ${SAMPLE_ID} -o ${OUTDIR_NAME} -q ${READ_TYPE} -oh ${OVERLAP} -v ${OPT_ARGS}
  
  if [[ $MODE = *unaligned* ]]
  then
    if [ "$DO_CONVERT" = true ]
    then
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/denovo/${SAMPLE_ID}_denovo_output.sam
    fi
    
    if [ "$DO_CONVERT_2" = true ]
    then
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/denovo/${SAMPLE_MATE_ID}_denovo_output.sam
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/genome/${ALIGNED_SAMPLE_MATE_ID}_genome_output.sam
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${ALIGNED_SAMPLE_MATE_ID}_junction_output.sam
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${ALIGNED_SAMPLE_MATE_ID}_reg_output.sam 
    fi
  else
    # delete those sam files when we're done with them
    if [ "$DO_CONVERT" = true ]
    then 
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${SAMPLE_ID}_junction_output.sam
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${SAMPLE_ID}_reg_output.sam
    fi
    if [ "$DO_CONVERT_2" = true ]
    then
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/genome/${SAMPLE_MATE_ID}_genome_output.sam
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/junction/${SAMPLE_MATE_ID}_junction_output.sam
      rm ${ALIGN_PARDIR}/${DATASET_NAME}/orig/reg/${SAMPLE_MATE_ID}_reg_output.sam
    fi
  fi
fi
