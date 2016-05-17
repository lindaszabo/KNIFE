#!/bin/sh

## Calls Bowtie2 to align a single fastq file to a single index, constructing commandline arguments
# for Bowtie based on the index being aligned to. SAM files are output by default, but if user
# specified "bam" or "bam_sort" in the mode then bam files are output.

# store commandline params and look up job info in txt file
CLUSTER_TYPE=$1
TASK_DATA_FILE=$2
ALIGN_PARDIR=$3
DATASET_NAME=$4
MODE=$5
INDEX=$6  # name of the bowtie index passed to bowtie2 call
ALIGN_OUTDIR=$7  # name of directory under orig to output alignment files, also used in output file names
FASTA=$8 # name of fasta file under index to use

source ../sampleInfo.sh ${CLUSTER_TYPE} # get sample-specific variables from TASK_DATA_FILE
source ../depends.sh ${CLUSTER_TYPE} # get sample-specific variables from TASK_DATA_FILE

OUTFILE_BASE=${ALIGN_PARDIR}/${DATASET_NAME}/orig/${ALIGN_OUTDIR}/${SAMPLE_ID}_${ALIGN_OUTDIR}_output

NFLAGS="--threads 4 --no-sq --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50 --no-unal" # flags common to all alignments

if [ "$ALIGN_OUTDIR" = "junction" -o "$ALIGN_OUTDIR" = "reg" -o "$ALIGN_OUTDIR" = "denovo" ]
then
  # no N penalty or cap on N
  NFLAGS=${NFLAGS}" --n-ceil L,0,1 --np 1"
fi

# For unaligned reads we report all reads. Only the first will be used for further analysis in the pipeline though,
# the additional alignments are present in the sam file only for manual exploration.
if [[ "$MODE" = *unaligned* ]]
then
  NFLAGS=${NFLAGS}" -a --un ${ALIGN_PARDIR}/${DATASET_NAME}/orig/still_unaligned/${SAMPLE_ID}_${ALIGN_OUTDIR}_output.fq"
  INDEX=${ALIGN_PARDIR}/${DATASET_NAME}/logs/denovo_index/${INDEX}
  FASTA=${ALIGN_PARDIR}/${DATASET_NAME}/logs/denovo_index/${FASTA}
fi

if [[ "$MODE" = *phred64* ]]
then
  NFLAGS=${NFLAGS}" --phred64"
fi

bowtie2 -x ${INDEX} -U ${READ_FILE} -S ${OUTFILE_BASE}.sam ${NFLAGS}
  
# save output as bam file of aligned reads and remove sam file
if [[ "$MODE" = *bam* ]]
then
  samtools view -bT ${FASTA} ${OUTFILE_BASE}.sam > ${OUTFILE_BASE}_unsorted.bam
  
  # sorting takes a long time, so only sort if user requested
  if [ "$MODE" = "*bam_sort*" ]
  then
    samtools sort ${OUTFILE_BASE}_unsorted.bam ${OUTFILE_BASE}

    # remove intermediates. Although in theory I should have been able to pipe output and not
    # save intermediate files, the samtools calls would fail due to prematurely terminated
    # files so I just delete the sam and unsorted bam when done creating the sorted bam
    rm ${OUTFILE_BASE}_unsorted.bam
  else
    # rename the unsorted so we have consistent files names for downstream processing
    mv ${OUTFILE_BASE}_unsorted.bam ${OUTFILE_BASE}.bam
  fi
  
  rm ${OUTFILE_BASE}.sam
fi
