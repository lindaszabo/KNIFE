#!/bin/bash

# Generates reports on alignment statistics under sampleStats directory.

# create text files with total reads, total genomic reads, total ribosomal reads, reads per ribosomal subunit, reads per junctionbucket, HBB reads
CLUSTER_TYPE=$1
TASK_DATA_FILE=$2
PAR_DIR=$3
DATASET_NAME=$4
REPORTDIR_NAME=$5
OUTDIR_NAME=$6
NTRIM=$7

if [ $# -eq 8 ]
then
  JUNC_DIRNAME=junction$8
  OUTF=SampleAlignmentStats$8.txt
  OUT_READS=SampleCircularStats$8.txt
else
  JUNC_DIRNAME=junction
  OUTF=SampleAlignmentStats.txt
  OUT_READS=SampleCircularStats.txt
fi

ORIG_DIRNAME=orig
RIBO_DIRNAME=ribo
GENOME_DIRNAME=genome
REG_DIRNAME=reg
ID_DIRNAME=ids

SAM_DIR=$PAR_DIR/$DATASET_NAME/$ORIG_DIRNAME
ID_DIR=$SAM_DIR/$ID_DIRNAME
RIBO_ID_DIR=$ID_DIR/$RIBO_DIRNAME
GENOME_ID_DIR=$ID_DIR/$GENOME_DIRNAME
JUNC_ID_DIR=$ID_DIR/$JUNC_DIRNAME
REG_ID_DIR=$ID_DIR/$REG_DIRNAME

source ../sampleInfo.sh ${CLUSTER_TYPE}  # get sample-specific variables from TASK_DATA_FILE
source ../depends.sh ${CLUSTER_TYPE}  # load python module if necessary

# get number of reads by counting lines in fastq file and dividing by 4 (handles gzipped files too)
if file --mime-type $READ_FILE | grep -q gzip$
then
  TOTAL_READ_LINES=`zcat $READ_FILE | wc -l`
else
  TOTAL_READ_LINES=`wc -l $READ_FILE | awk '{print $1}'`
fi
TOTAL_READS=`echo "scale=0; $TOTAL_READ_LINES/4" | bc -l`
  
GENOME_FILE=${ID_DIR}/${GENOME_DIRNAME}/${SAMPLE_ID}_genome_output.txt  
TOTAL_GENOME=`wc -l $GENOME_FILE | awk '{print $1}'`
PERCENT_GENOME=`echo "scale=4; $TOTAL_GENOME/$TOTAL_READS" | bc -l`
GENOME_POS=`awk '$2 == 0 {print $2}' $GENOME_FILE | wc -l`
GENOME_POS=`echo "scale=2; $GENOME_POS/$TOTAL_GENOME" | bc -l`
GENOME_NEG=`awk '$2 == 16 {print $2}' $GENOME_FILE | wc -l`
GENOME_NEG=`echo "scale=2; $GENOME_NEG/$TOTAL_GENOME" | bc -l`
  
# count reads mapping with offset within HBB CDS (chr11 617350-621570) from ucsc hg19 knowngenes file
HBB_READS=`awk '$3 == "chr11" && $4 >= 617350 && $4 <= 621570 {print $1}' $GENOME_FILE | wc -l`
  
JUNC_FILE=${ID_DIR}/${JUNC_DIRNAME}/${SAMPLE_ID}_junction_output.txt
TOTAL_JUNC=`wc -l $JUNC_FILE | awk '{print $1}'`
PERCENT_JUNC=`echo "scale=4; $TOTAL_JUNC/$TOTAL_READS" | bc -l`
# since - strand junctions are actually stored as reverse complement in index, neg alignments to - strand genes are actually pos alignments
JUNC_POS=`awk '($3 ~ /+$/ && $2 == 0) || ($3 ~ /-$/ && $2 == 16) {print $2}' $JUNC_FILE | wc -l`
JUNC_POS=`echo "scale=2; $JUNC_POS/$TOTAL_JUNC" | bc -l`
JUNC_NEG=`awk '($3 ~ /+$/ && $2 == 16) || ($3 ~ /-$/ && $2 == 0) {print $2}' $JUNC_FILE | wc -l`
JUNC_NEG=`echo "scale=2; $JUNC_NEG/$TOTAL_JUNC" | bc -l`
  
REG_FILE=${ID_DIR}/${REG_DIRNAME}/${SAMPLE_ID}_reg_output.txt
TOTAL_REG=`wc -l $REG_FILE | awk '{print $1}'`
PERCENT_REG=`echo "scale=4; $TOTAL_REG/$TOTAL_READS" | bc -l`
REG_POS=`awk '($3 ~ /+$/ && $2 == 0) || ($3 ~ /-$/ && $2 == 16) {print $2}' $REG_FILE | wc -l`
REG_POS=`echo "scale=2; $REG_POS/$TOTAL_REG" | bc -l`
REG_NEG=`awk '($3 ~ /+$/ && $2 == 16) || ($3 ~ /-$/ && $2 == 0) {print $2}' $REG_FILE | wc -l`
REG_NEG=`echo "scale=2; $REG_NEG/$TOTAL_REG" | bc -l`
  
RIBO_FILE=${ID_DIR}/${RIBO_DIRNAME}/${SAMPLE_ID}_ribo_output.txt
TOTAL_RIBO=`wc -l $RIBO_FILE | awk '{print $1}'`
PERCENT_RIBO=`echo "scale=4; $TOTAL_RIBO/$TOTAL_READS" | bc -l`
RIBO_POS=`awk '$2 == 0 {print $2}' $RIBO_FILE | wc -l`
RIBO_POS=`echo "scale=2; $RIBO_POS/$TOTAL_RIBO" | bc -l`
RIBO_NEG=`awk '$2 == 16 {print $2}' $RIBO_FILE | wc -l`
RIBO_NEG=`echo "scale=2; $RIBO_NEG/$TOTAL_RIBO" | bc -l`
  
# parse out the reads per ribo subunit
RIBO_28S=`grep "gi|224514641:112587-119176" $RIBO_FILE | wc -l`
RIBO_18S=`grep "gi|224514641:109078-110946" $RIBO_FILE | wc -l`
RIBO_58S=`grep "gi|224514641:112025-112180" $RIBO_FILE | wc -l`
RIBO_5SDNA=`grep "gi|23898|emb|X12811.1|" $RIBO_FILE | wc -l`
RIBO_5SrRNA=`grep "hg19_5SrRna" $RIBO_FILE | wc -l`
RIBO_28S=`echo "scale=4; $RIBO_28S/$TOTAL_RIBO" | bc -l`
RIBO_18S=`echo "scale=4; $RIBO_18S/$TOTAL_RIBO" | bc -l`
RIBO_58S=`echo "scale=4; $RIBO_58S/$TOTAL_RIBO" | bc -l`
RIBO_5SDNA=`echo "scale=4; $RIBO_5SDNA/$TOTAL_RIBO" | bc -l`
RIBO_5SrRNA=`echo "scale=4; $RIBO_5SrRNA/$TOTAL_RIBO" | bc -l`
  
# get count of reads that did not align anywhere
# first convert bam to sam if necessary, just check one file and assume all sam files exist if one does
f=`find ${SAM_DIR}/${GENOME_DIRNAME} -type f -name ${SAMPLE_ID}_${GENOME_DIRNAME}_output.sam`
if [ ! -f "$f" ]
then
  DO_CONVERT=true
  echo "Did not find sam files, need to convert bam to sam for further processing"
  f=`find ${SAM_DIR}/${GENOME_DIRNAME} -type f -name ${SAMPLE_ID}_${GENOME_DIRNAME}_output.bam`
  # first convert to sam file for reading
  samtools view ${f} > ${SAM_DIR}/${GENOME_DIRNAME}/${SAMPLE_ID}_${GENOME_DIRNAME}_output.sam
  samtools view ${SAM_DIR}/${JUNC_DIRNAME}/${SAMPLE_ID}_${JUNC_DIRNAME}_output.bam > ${SAM_DIR}/${JUNC_DIRNAME}/${SAMPLE_ID}_${JUNC_DIRNAME}_output.sam
  samtools view ${SAM_DIR}/${REG_DIRNAME}/${SAMPLE_ID}_${REG_DIRNAME}_output.bam > ${SAM_DIR}/${REG_DIRNAME}/${SAMPLE_ID}_${REG_DIRNAME}_output.sam
  samtools view ${SAM_DIR}/${RIBO_DIRNAME}/${SAMPLE_ID}_${RIBO_DIRNAME}_output.bam > ${SAM_DIR}/${RIBO_DIRNAME}/${SAMPLE_ID}_${RIBO_DIRNAME}_output.sam
fi
  

UNALIGNED=`python getUnalignedReadCount.py -r ${READ_FILE} -n ${SAMPLE_ID} -a ${SAM_DIR} -t ${NTRIM}`
PERCENT_UNALIGNED=`echo "scale=4; $UNALIGNED/$TOTAL_READS" | bc -l`
  
if [ "$DO_CONVERT" = true ]
then 
  # and then remove the temp file
  rm ${SAM_DIR}/${GENOME_DIRNAME}/${SAMPLE_ID}_${GENOME_DIRNAME}_output.sam
  rm ${SAM_DIR}/${JUNC_DIRNAME}/${SAMPLE_ID}_${JUNC_DIRNAME}_output.sam
  rm ${SAM_DIR}/${REG_DIRNAME}/${SAMPLE_ID}_${REG_DIRNAME}_output.sam
  rm ${SAM_DIR}/${RIBO_DIRNAME}/${SAMPLE_ID}_${RIBO_DIRNAME}_output.sam
fi

outfile="${OUTDIR_NAME}/${SAMPLE_ID}${OUTF}"
echo -e "$SAMPLE_ID\t$TOTAL_READS\t$UNALIGNED ($PERCENT_UNALIGNED)\t$TOTAL_GENOME ($PERCENT_GENOME)\t$GENOME_POS, $GENOME_NEG\t$TOTAL_REG ($PERCENT_REG)\t$REG_POS, $REG_NEG\t$TOTAL_JUNC ($PERCENT_JUNC)\t$JUNC_POS, $JUNC_NEG\t$TOTAL_RIBO ($PERCENT_RIBO)\t$RIBO_POS, $RIBO_NEG\t$RIBO_28S\t$RIBO_18S\t$RIBO_58S\t$RIBO_5SDNA\t$RIBO_5SrRNA\t$HBB_READS" > $outfile


# generate the read assignment stats for R1 only
READ_NUM=`echo ${SAMPLE_ID:(-1)}` # will be a 1 or 2
if [ "$READ_NUM" -eq 1 ]
then
  REPORT_DIR=${PAR_DIR}/${DATASET_NAME}/${REPORTDIR_NAME}/ids/${SAMPLE_ID}__output.txt
  REPORT_FILE=`find $REPORT_DIR -type f -name "${SAMPLE_ID}__output.txt"` 

  if [ -f "${REPORT_FILE}" ]
  then
    CIRC_READS=`awk '$2 == "circStrong" {print $4}' ${REPORT_FILE} | wc -l`
    CIRC_ARTIFACT=`awk '$2 == "circArtifact" {print $4}' ${REPORT_FILE} | wc -l`
    DECOY_READS=`awk '$2 == "decoy" {print $4}' ${REPORT_FILE} | wc -l`
    LINEAR_READS=`awk '$2 == "linearStrong" {print $4}' ${REPORT_FILE} | wc -l`
    LINEAR_ARTIFACT=`awk '$2 == "linearArtifact" {print $4}' ${REPORT_FILE} | wc -l`
    ANOM_READS=`awk '$2 == "anomaly" {print $4}' ${REPORT_FILE} | wc -l`
    UN_READS=`awk '$2 == "unmapped" {print $4}' ${REPORT_FILE} | wc -l`
    NON_GR=`expr $CIRC_READS + $CIRC_ARTIFACT + $DECOY_READS + $LINEAR_READS + $LINEAR_ARTIFACT + $ANOM_READS + $UN_READS`
    CIRC_FRACTION=`echo "scale=6; $CIRC_READS/$NON_GR" | bc -l`
    LINEAR_FRACTION=`echo "scale=6; $LINEAR_READS/$NON_GR" | bc -l`
    CIRC_OVER_LINEAR=`echo "scale=6; $CIRC_FRACTION/$LINEAR_FRACTION" | bc -l`

    outreadsfile="${OUTDIR_NAME}/${SAMPLE_ID}${OUT_READS}"
    echo -e "$SAMPLE_ID\t$CIRC_READS\t$CIRC_ARTIFACT\t$DECOY_READS\t$LINEAR_READS\t$LINEAR_ARTIFACT\t$ANOM_READS\t$UN_READS\t$NON_GR\t$CIRC_FRACTION\t$LINEAR_FRACTION\t$CIRC_OVER_LINEAR" >> $outreadsfile
  fi
fi
