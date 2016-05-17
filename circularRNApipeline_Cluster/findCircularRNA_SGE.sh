#!/bin/sh

## The main wrapper script to generate reports containing read count and statistical score
# for each circular or linear junction detected in the data. Here, the task array feature
# of Sun Grid Engine is used to run tasks in parallel on individual files where possible,
# waiting to proceed until all files successfully processed at a specific step when necessary.
# These memory / time requests work for our cluster, you may need to edit them to meet
# restrictions implemented on your cluster.
# The main steps are:
#  1) align to fastq files to genome, linear junction, circular junction, and ribosomal indices using Bowtie2 
#  2) parse sam files and analyze alignments using "naive" statistical method
#  3) run GLM to report posterior probability and p-value per junction
#  4) generate alignment statistics per sample (useful for evaluating efficiency of ribo depletion, global circular:linear ratio, etc)
#     and identify unaligned reads 
#  5) identify candidate de novo junctions from unaligned reads and generate a Bowtie2 de novo index from these candidates
#  6) align unaligned reads to de novo index and analyze using "naive" statistical method (requires 2nd call to findCircularRNA.sh using "unaligned" mode)

# Usage (see README for descriptions of parameters):
#    sh findCircularRNA_SGE.sh read_directory read_id_style alignment_parent_directory dataset_name junction_overlap [mode] [report_directory_name] [ntrim] [denovoCircMode] [junction_id_suffix] 
#    And then run again with same parameters as before but append "_unaligned" to the mode parameter

# very basic error checking
if [ $# -eq 11 ]
then
  echo "Read 1 threshold was specified without Read 2 threshold. Please include Read 2 threshold or check your arguments if you did not intend to specify Read 1 threshold ${11}"
  exit 2
fi

CODE_DIR=`pwd`
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

if [ $# -ge 13 ]
then
  JUNCTION_MIDPOINT=${13}
else
  JUNCTION_MIDPOINT=150
fi

# change resource allocations based on job size
if [[ "$MODE" = *large* ]]
then
  JUNC_VMEM="80G"
  GENOME_VMEM="80G"
  REG_VMEM="80G"
  RIBO_VMEM="80G"
  ALIGN_MAX_RT="23:0:0"
  PREPROCESS_VMEM="80G"
  FILTER_VMEM="80G"
  PREPROCESS_MAX_RT="23:0:0"
  FILTER_MAX_RT="23:0:0"
  PFA_VMEM="80G"
  QUALITY_VMEM="80G"
  QUALITY_MAX_RT="48:0:0"
else
  JUNC_VMEM="20G"
  GENOME_VMEM="6G"
  REG_VMEM="20G"
  RIBO_VMEM="5G"
  ALIGN_MAX_RT="12:0:0"
  PREPROCESS_VMEM="15G"
  FILTER_VMEM="30G"
  PREPROCESS_MAX_RT="4:0:0"
  FILTER_MAX_RT="12:0:0"
  PFA_VMEM="30G"
  QUALITY_VMEM="20G"
  QUALITY_MAX_RT="48:0:0"
fi

if [[ "$MODE" = *bam* ]]
then
  ALIGN_MAX_RT="23:0:0"
  PREPROCESS_MAX_RT="23:0:0"
  FILTER_MAX_RT="23:0:0"
fi

# a few special flags for the unaligned mode
if [[ $MODE = *unaligned* ]]
then
  TASK_FILE_READ_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/orig/unaligned
  UFLAG="-u"
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}_unaligned.txt
else
  TASK_FILE_READ_DIR=${READ_DIR}
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt
fi  

# set up info about this dataset, create directory structures
if [[ $MODE != *analysis* ]]
then
  # set up
  python analysis/writeTaskIdFiles.py -r ${TASK_FILE_READ_DIR} -a ${ALIGN_PARDIR} -d ${DATASET_NAME} ${UFLAG}
  
  # select correct prefix name to use for bowtie index files
  if [[ $MODE = *grch38* ]]
  then
    bt_prefix="grch38"
  elif [[ $MODE = *mouse* ]]
  then
    bt_prefix="mm10"
  elif [[ $MODE = *rat* ]]
  then
    bt_prefix="rn5"
  elif [[ $MODE = *fly* ]]
  then
    bt_prefix="dm3"
  else
    bt_prefix="hg19"
  fi
fi

# break out of alignment-only loop since this is needed for both alignment and analysis
NUM_FILES=`cat $TASK_DATA_FILE | wc -l`

# create directory and subdirectories for log output
LOG_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/logs
mkdir -p ${LOG_DIR}
mkdir -p ${LOG_DIR}/align
mkdir -p ${LOG_DIR}/analysis
mkdir -p ${LOG_DIR}/glm
mkdir -p ${LOG_DIR}/sampleStats
mkdir -p ${LOG_DIR}/denovo_script_out
mkdir -p ${LOG_DIR}/denovo_index

# alignments
if [[ $MODE != *analysis* ]]
then
  if [[ $MODE = *unaligned* ]]
  then
    qsub -t 1-${NUM_FILES}:1 -N DeNovoAlign${DATASET_NAME}${DENOVOCIRC} -l h_vmem=6G -l h_rt=${ALIGN_MAX_RT} -wd ${CODE_DIR}/index -o ${LOG_DIR}/align/ -e ${LOG_DIR}/align/ analysis/align.sh SGE $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE denovo_${DATASET_NAME}_${DENOVOCIRC} denovo denovo_${DATASET_NAME}_onlycircles${DENOVOCIRC}.fa
    depend_str="-hold_jid_ad "DeNovoAlign${DATASET_NAME}${DENOVOCIRC}
  else
    qsub -t 1-${NUM_FILES}:1 -N GenomeAlign${DATASET_NAME} -l h_vmem=${GENOME_VMEM} -l h_rt=${ALIGN_MAX_RT} -wd ${CODE_DIR}/index -o ${LOG_DIR}/align/ -e ${LOG_DIR}/align/ analysis/align.sh SGE $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_genome genome ${bt_prefix}_genome.fa
    qsub -t 1-${NUM_FILES}:1 -N JunctionAlign${DATASET_NAME} -l h_vmem=${JUNC_VMEM} -l h_rt=${ALIGN_MAX_RT} -wd ${CODE_DIR}/index -o ${LOG_DIR}/align/ -e ${LOG_DIR}/align/ analysis/align.sh SGE $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_junctions_scrambled junction ${bt_prefix}_junctions_scrambled.fa
    qsub -t 1-${NUM_FILES}:1 -N RiboAlign${DATASET_NAME} -l h_vmem=${RIBO_VMEM} -l h_rt=${ALIGN_MAX_RT} -wd ${CODE_DIR}/index -o ${LOG_DIR}/align/ -e ${LOG_DIR}/align/ analysis/align.sh SGE $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_ribosomal ribo ${bt_prefix}_ribosomal.fa
    qsub -t 1-${NUM_FILES}:1 -N RegAlign${DATASET_NAME} -l h_vmem=${REG_VMEM} -l h_rt=${ALIGN_MAX_RT} -wd ${CODE_DIR}/index -o ${LOG_DIR}/align/ -e ${LOG_DIR}/align/ analysis/align.sh SGE $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_junctions_reg reg ${bt_prefix}_junctions_reg.fa
    depend_str="-hold_jid_ad "GenomeAlign${DATASET_NAME},JunctionAlign${DATASET_NAME},RiboAlign${DATASET_NAME},RegAlign${DATASET_NAME}
  fi
fi

if [[ $MODE = *R2analysis* ]]
then
  # set up directories and swap alignment file names to run with R2 as R1 
  qsub -N CSD${DATASET_NAME} -l h_vmem=10G -l h_rt=5:0:0 -wd ${CODE_DIR}/analysis -o ${LOG_DIR}/analysis/${DATASET_NAME}CSD.out -e ${LOG_DIR}/analysis/${DATASET_NAME}CSD.err analysis/createSwappedDirectories.sh SGE $ALIGN_PARDIR $DATASET_NAME swap 
  
  UNSWAPPED_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}
  UNSWAPPED_DATASET_NAME=${DATASET_NAME}
  DATASET_NAME=${DATASET_NAME}Swapped  # from here on out, the analysis mode will operate on the swapped directory
  SWAPPED_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}
  
  # if we ran previously in unaligned mode, we will want to analyze that output too. But no errors thrown if you didn't run in unaligned mode previously
  DENOVO_TASK_DATA_FILE=${DATASET_NAME}_unaligned.txt  # this is the to-be-created-below file
  f=`find ${ALIGN_PARDIR}/taskIdFiles -type f -name ${UNSWAPPED_DATASET_NAME}_unaligned.txt`
  if [ ! -f "$f" ]
  then
    NUM_DENOVO_FILES=0
  else
    # make a copy of the unaligned task data file with the correct name for the swapped run
    cp ${ALIGN_PARDIR}/taskIdFiles/${UNSWAPPED_DATASET_NAME}_unaligned.txt ${ALIGN_PARDIR}/taskIdFiles/${DENOVO_TASK_DATA_FILE}  
    NUM_DENOVO_FILES=`cat ${ALIGN_PARDIR}/taskIdFiles/${DENOVO_TASK_DATA_FILE} | wc -l`
  fi
  depend_str="-hold_jid "CSD${UNSWAPPED_DATASET_NAME}
fi

# preprocessing
qsub -t 1-${NUM_FILES}:1 -N ${DATASET_NAME}Preprocess -l h_vmem=${PREPROCESS_VMEM} -l h_rt=${PREPROCESS_MAX_RT} -wd ${CODE_DIR} -o ${LOG_DIR}/analysis/ -e ${LOG_DIR}/analysis/ ${depend_str} analysis/preprocessAlignedReads.sh SGE ${TASK_DATA_FILE} ${ALIGN_PARDIR} ${DATASET_NAME} ${JUNCTION_MIDPOINT} ${OVERLAP} ${JUNCTION_DIR_SUFFIX}

# analysis
qsub -t 1-${NUM_FILES}:1 -N ${DATASET_NAME}Analysis -l h_vmem=${FILTER_VMEM} -l h_rt=${FILTER_MAX_RT} -wd ${CODE_DIR}/analysis -o ${LOG_DIR}/analysis/ -e ${LOG_DIR}/analysis/ -hold_jid_ad ${DATASET_NAME}Preprocess analysis/filterFDR.sh SGE ${MODE} ${ALIGN_PARDIR} ${DATASET_NAME} ${REPORTDIR_NAME} ${READ_STYLE} ${OVERLAP} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH}

# if this is R2analysis mode then we need to re-call preprocessing and analysis with the unaligned task id file
if [[ $MODE = *R2analysis* ]]
then
  if [ ${NUM_DENOVO_FILES} -gt 0 ]
  then
    qsub -t 1-${NUM_DENOVO_FILES} -N ${DATASET_NAME}UnalPreprocess -l h_vmem=${PREPROCESS_VMEM} -l h_rt=${PREPROCESS_MAX_RT} -wd ${CODE_DIR} -o ${LOG_DIR}/analysis/ -e ${LOG_DIR}/analysis/ ${depend_str} analysis/preprocessAlignedReads.sh SGE ${ALIGN_PARDIR}/taskIdFiles/${DENOVO_TASK_DATA_FILE} ${ALIGN_PARDIR} ${DATASET_NAME} ${JUNCTION_MIDPOINT} ${OVERLAP} ${JUNCTION_DIR_SUFFIX} 
    qsub -t 1-${NUM_DENOVO_FILES} -N ${DATASET_NAME}UnalAnalysis -l h_vmem=${FILTER_VMEM} -l h_rt=${FILTER_MAX_RT} -wd ${CODE_DIR}/analysis -o ${LOG_DIR}/analysis/ -e ${LOG_DIR}/analysis/ -hold_jid_ad ${DATASET_NAME}UnalPreprocess analysis/filterFDR.sh SGE ${MODE}_unaligned ${ALIGN_PARDIR} ${DATASET_NAME} ${REPORTDIR_NAME} ${READ_STYLE} ${OVERLAP} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} 
  fi
fi

# quality stats, glm, and denovo only if this is not unaligned mode. if this is R2analysis, then quality stats and denovo junction creation do not get run 
if [[ ${MODE} != *unaligned* ]]
then
  
  if [[ ${MODE} != *skipGLM* ]]
  then
    qsub -t 1-${NUM_FILES}:1 -N PFA${DATASET_NAME} -l h_vmem=${PFA_VMEM} -l h_rt=${PREPROCESS_MAX_RT} -wd ${CODE_DIR}/analysis -o ${LOG_DIR}/glm/ -e ${LOG_DIR}/glm/ -hold_jid_ad ${DATASET_NAME}Analysis analysis/predictJunctions.sh SGE ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${REPORTDIR_NAME}
  fi
  
  if [[ ${MODE} != *R2analysis* ]]
  then
    # make the output directories
    OUT_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/sampleStats
    mkdir -p ${OUT_DIR}
    mkdir -p ${ALIGN_PARDIR}/${DATASET_NAME}/orig/unaligned
    mkdir -p ${ALIGN_PARDIR}/${DATASET_NAME}/orig/unaligned/forDenovoIndex
    
    # for each sample, generate sample stats in the output directory
    qsub -t 1-${NUM_FILES}:1 -N Quality${DATASET_NAME} -l h_vmem=${QUALITY_VMEM} -l h_rt=${QUALITY_MAX_RT} -wd ${CODE_DIR}/qualityStats -o ${LOG_DIR}/sampleStats/ -e ${LOG_DIR}/sampleStats/ -hold_jid_ad ${DATASET_NAME}Analysis qualityStats/qualityStatsSingleSample.sh SGE $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $REPORTDIR_NAME ${OUT_DIR} ${NTRIM} ${JUNCTION_DIR_SUFFIX} 
    
    # cat all of those files then delete the original separate files
    qsub -N Cat${DATASET_NAME} -l h_vmem=5G -l h_rt=1:0:0 -wd ${CODE_DIR}/qualityStats -o ${LOG_DIR}/sampleStats/ -e ${LOG_DIR}/sampleStats/ -hold_jid Quality${DATASET_NAME} qualityStats/qualityStatsCat.sh ${OUT_DIR}

    if [[ ${MODE} != *skipDenovo* ]]    # create the index for this dataset
    then
      qsub -N Denovo_${DENOVOCIRC}_${DATASET_NAME} -l h_vmem=55G -l h_rt=23:0:0 -wd ${CODE_DIR}/denovo_scripts -o ${LOG_DIR}/denovo_script_out/ -e ${LOG_DIR}/denovo_script_out/ -hold_jid Quality${DATASET_NAME} denovo_scripts/denovo_scripts.batch ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${NTRIM} ${DENOVOCIRC} SGE
    fi
  fi
fi

if [[ ${MODE} = *R2analysis* ]]
then
  if [ ${NUM_DENOVO_FILES} -gt 0 ]
  then
    # this will be added to the dependency for the annotation-dependent analysis
    naive_depend_str=","${DATASET_NAME}UnalAnalysis
  fi
  # combine the swapped reads output and initial read output into single output file
  qsub -N CSR${DATASET_NAME} -l h_vmem=25G -l h_rt=23:0:0 -wd ${CODE_DIR}/analysis -o ${LOG_DIR}/analysis/${DATASET_NAME}CSR.out -e ${LOG_DIR}/analysis/${DATASET_NAME}CSR.err -hold_jid PFA${DATASET_NAME} analysis/combineSwappedReadsGLM.sh SGE ${UNSWAPPED_DIR} ${SWAPPED_DIR} ${READ_STYLE}
  qsub -N CSN${DATASET_NAME} -l h_vmem=25G -l h_rt=23:0:0 -wd ${CODE_DIR}/analysis -o ${LOG_DIR}/analysis/${DATASET_NAME}CSN.osut -e ${LOG_DIR}/analysis/${DATASET_NAME}CSN.err -hold_jid ${DATASET_NAME}Analysis${naive_depend_str} analysis/combineSwappedReadsNaive.sh SGE ${UNSWAPPED_DIR} ${SWAPPED_DIR} ${READ_STYLE}
  # swap back all of the alignment files to their original directories. this will not run if either of the above commands failed, will need to run manually in that case
  qsub -N CSD2${DATASET_NAME} -l h_vmem=10G -l h_rt=5:0:0 -wd ${CODE_DIR}/analysis -o ${LOG_DIR}/analysis/${DATASET_NAME}CSD2.out -e ${LOG_DIR}/analysis/${DATASET_NAME}CSD2.err -hold_jid CSR${DATASET_NAME},CSN${DATASET_NAME} analysis/createSwappedDirectories.sh SGE ${ALIGN_PARDIR} ${UNSWAPPED_DATASET_NAME} restore
fi

