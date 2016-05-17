#!/bin/sh

# Create text files with ids of Ribo-aligned, genome-aligned, 
# junction-overlapped, and reg junction-overlapped reads.
# Assumes the sam alignment files are written to junction, ribo, reg, and genome directories.

# it will first look for sam file for this sample, if none found then it will try to convert bam to sam

# store commandline params and look up job info in txt file
CLUSTER_TYPE=$1
TASK_DATA_FILE=$2
PAR_DIR=$3  
DATASET_NAME=$4 
MIDPOINT=$5
OVERHANG=$6

# declare some variables for the directory structure so I find & put files in correct location
ORIG_DIRNAME=orig
DENOVO_DIRNAME=denovo
RIBO_DIRNAME=ribo
JUNC_DIRNAME=junction
GENOME_DIRNAME=genome
REG_DIRNAME=reg
ID_DIRNAME=ids
JUNC_NONGR_DIRNAME=juncNonGR
DENOVO_NONGR_DIRNAME=denovoNonGR

ALIGN_DIR=$PAR_DIR/$DATASET_NAME/$ORIG_DIRNAME
DENOVO_ALIGN_DIR=$ALIGN_DIR/$DENOVO_DIRNAME
RIBO_ALIGN_DIR=$ALIGN_DIR/$RIBO_DIRNAME
JUNC_ALIGN_DIR=$ALIGN_DIR/$JUNC_DIRNAME
GENOME_ALIGN_DIR=$ALIGN_DIR/$GENOME_DIRNAME
REG_ALIGN_DIR=$ALIGN_DIR/$REG_DIRNAME

ID_DIR=$ALIGN_DIR/$ID_DIRNAME
RIBO_ID_DIR=$ID_DIR/$RIBO_DIRNAME
GENOME_ID_DIR=$ID_DIR/$GENOME_DIRNAME

# if a directory suffix was passed, we want to output the junction overlaps to a different
# directory than the default
if [ $# -eq 7 ]
then
  DENOVO_ID_DIR=$ID_DIR/$DENOVO_DIRNAME$7
  mkdir -p $DENOVO_ID_DIR
  JUNC_ID_DIR=$ID_DIR/$JUNC_DIRNAME$7
  mkdir -p $JUNC_ID_DIR
  REG_ID_DIR=$ID_DIR/$REG_DIRNAME$7
  mkdir -p $REG_ID_DIR
  mkdir -p $ID_DIR/$JUNC_NONGR_DIRNAME$7
else
  DENOVO_ID_DIR=$ID_DIR/$DENOVO_DIRNAME
  JUNC_ID_DIR=$ID_DIR/$JUNC_DIRNAME
  REG_ID_DIR=$ID_DIR/$REG_DIRNAME
fi

source ./sampleInfo.sh ${CLUSTER_TYPE} # get sample-specific variables from TASK_DATA_FILE

# just need to to the unaligned reads if sample id starts with unaligned_
if [[ "$SAMPLE_ID" = unaligned_* ]]
then
  # get read ids that overlap the denovo junctions
  f=`find ${DENOVO_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${DENOVO_DIRNAME}_output.sam`
  if [ ! -f "$f" ]
  then
    DO_CONVERT=true
    # first convert to sam file for reading
    f=`find ${DENOVO_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${DENOVO_DIRNAME}_output.bam`
    if [ -f "$f" ]
    then
      samtools view ${f} > ${DENOVO_ALIGN_DIR}/${SAMPLE_ID}_${DENOVO_DIRNAME}_output.sam
      f=`find ${DENOVO_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${DENOVO_DIRNAME}_output.sam`
    fi
  fi

  if [ -f "$f" ]
  then
    filename=$(basename "$f")
    filename="${filename%.*}"
    outfile="${DENOVO_ID_DIR}/${filename}.txt"
    awk '$1 !~ /^@/ && $4 >= '$MIDPOINT'-length($10)+'$OVERHANG'+1 && $4 <= '$MIDPOINT'-'$OVERHANG'+1 {print $1 "\t" $2 "\t" $3}' $f > $outfile

    if [ "$DO_CONVERT" = true ]
    then
      # and then remove the temp file
      rm ${f}
    fi
  fi
else
  # this is a regular first-time run, preprocess all of the genome, junction, ribo
  # get read ids that overlap a junction
  f=`find ${JUNC_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${JUNC_DIRNAME}_output.sam`
  if [ ! -f "$f" ]
  then
    DO_CONVERT=true
    f=`find ${JUNC_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${JUNC_DIRNAME}_output.bam`
    # first convert to sam file for reading
    samtools view ${f} > ${JUNC_ALIGN_DIR}/${SAMPLE_ID}_${JUNC_DIRNAME}_output.sam
    f=`find ${JUNC_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${JUNC_DIRNAME}_output.sam`
  fi
  filename=$(basename "$f")
  filename="${filename%.*}"
  outfile="${JUNC_ID_DIR}/${filename}.txt"
  awk '$1 !~ /^@/ && $4 >= '$MIDPOINT'-length($10)+'$OVERHANG'+1 && $4 <= '$MIDPOINT'-'$OVERHANG'+1 {print $1 "\t" $2 "\t" $3}' $f > $outfile

  if [ "$DO_CONVERT" = true ]
  then 
    # and then remove the temp file
    rm ${f}
  fi

  # get read ids that overlap a regular junction
  f=`find ${REG_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${REG_DIRNAME}_output.sam`
  if [ ! -f "$f" ]
  then
    f=`find ${REG_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${REG_DIRNAME}_output.bam`
    # first convert to sam file for reading
    samtools view ${f} > ${REG_ALIGN_DIR}/${SAMPLE_ID}_${REG_DIRNAME}_output.sam
    f=`find ${REG_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${REG_DIRNAME}_output.sam`
  fi
  filename=$(basename "$f")
  filename="${filename%.*}"
  outfile="${REG_ID_DIR}/${filename}.txt"
  awk '$1 !~ /^@/ && $4 >= '$MIDPOINT'-length($10)+'$OVERHANG'+1 && $4 <= '$MIDPOINT'-'$OVERHANG'+1 {print $1 "\t" $2 "\t" $3}' $f > $outfile

  if [ "$DO_CONVERT" = true ]
  then 
    # and then remove the temp file
    rm ${f}
  fi

  # get read ids that aligned to the ribosome
  f=`find ${RIBO_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${RIBO_DIRNAME}_output.sam`
  if [ ! -f "$f" ]
  then
    f=`find ${RIBO_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${RIBO_DIRNAME}_output.bam`
    # first convert to sam file for reading
    samtools view ${f} > ${RIBO_ALIGN_DIR}/${SAMPLE_ID}_${RIBO_DIRNAME}_output.sam
    f=`find ${RIBO_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${RIBO_DIRNAME}_output.sam`
  fi
  filename=$(basename "$f")
  filename="${filename%.*}"
  outfile="${RIBO_ID_DIR}/${filename}.txt"
  awk '$1 !~ /^@/ && $2 != 4 {print $1  "\t" $2 "\t" $3}' $f > $outfile

  if [ "$DO_CONVERT" = true ]
  then 
    # and then remove the temp file
    rm ${f}
  fi

  # get read ids that aligned to the genome
  f=`find ${GENOME_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${GENOME_DIRNAME}_output.sam`
  if [ ! -f "$f" ]
  then
    f=`find ${GENOME_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${GENOME_DIRNAME}_output.bam`
    # first convert to sam file for reading
    samtools view ${f} > ${GENOME_ALIGN_DIR}/${SAMPLE_ID}_${GENOME_DIRNAME}_output.sam
    f=`find ${GENOME_ALIGN_DIR} -type f -name ${SAMPLE_ID}_${GENOME_DIRNAME}_output.sam`
  fi
  filename=$(basename "$f")
  filename="${filename%.*}"
  outfile="${GENOME_ID_DIR}/${filename}.txt"
  awk '$1 !~ /^@/ && $2 != 4 {print $1  "\t" $2 "\t" $3 "\t" $4}' $f > $outfile

  if [ "$DO_CONVERT" = true ]
  then 
    # and then remove the temp file
    rm ${f}
  fi
fi
