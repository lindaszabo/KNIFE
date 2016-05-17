#!/bin/sh

# store commandline params and look up job info in txt file
READ_FILE=$1
SAMPLE_ID=$2
ALIGN_PARDIR=$3
DATASET_NAME=$4
MODE=$5
INDEX=$6  # name of the bowtie index passed to bowtie2 call
ALIGN_OUTDIR=$7  # name of directory under orig to output alignment files, also used in output file names
FASTA=$8 # name of fasta file under index to use

OUTFILE_BASE=${ALIGN_PARDIR}/${DATASET_NAME}/orig/${ALIGN_OUTDIR}/${SAMPLE_ID}_${ALIGN_OUTDIR}_output

NFLAGS="--threads 4 --no-sq --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50 --no-unal" # flags common to all alignments

if [ "$ALIGN_OUTDIR" = "junction" -o "$ALIGN_OUTDIR" = "reg" -o "$ALIGN_OUTDIR" = "denovo" ]
then
  # no N penalty or cap on N
  NFLAGS=${NFLAGS}" --n-ceil L,0,1 --np 1 -p 8"
fi

# for unaligned reads we report all reads. Only the first will be used for further analysis in the pipeline though.
# also output unaligned reads only for the unaligned mode
if [[ "$MODE" = *unaligned* ]]
then
  NFLAGS=${NFLAGS}" -a --un ${ALIGN_PARDIR}/${DATASET_NAME}/orig/still_unaligned/${SAMPLE_ID}_${ALIGN_OUTDIR}_output.fq"
fi

if [[ "$MODE" = *phred64* ]]
then
  NFLAGS=${NFLAGS}" --phred64"
fi

bowtie2 -x ${INDEX} -U ${READ_FILE} -S ${OUTFILE_BASE}.sam ${NFLAGS}
  
# save output as bam file of aligned reads and remove sam file
if [[ "$MODE" = *bam* ]]
then
  echo "outputting bam file: ${OUTFILE_BASE}_unsorted.bam"
  samtools view -bT ${FASTA} ${OUTFILE_BASE}.sam > ${OUTFILE_BASE}_unsorted.bam
  
  # sorting takes a long time, so only sort if user requested
  if [ "$MODE" = "*bam_sort*" ]
  then
    echo "sorting bam file: ${OUTFILE_BASE}_unsorted.bam"
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

