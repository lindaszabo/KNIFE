#!/bin/sh

# This assumes that createExonDB.py has already been called, creating a directory OUT_DIR that contains exons, genes, and records directories.
# This code now parses that exon info and creates junction fasta files for scrambled and linear junctions and bowtie2 indices and places these
# into the pipeline code so they can be used by specifying ${INDEX_MODE_ID} in the MODE parameter when calling findCircularRNA.sh.

# pipelineDirectory: path to the circularRNApipeline directory. Created files will be placed inside the index directory.
# outputDirectory: the output directory specified in the call to createExonDB.sh for this species
# fileIdentifier: should be short String to help distinguish it from other genomes. This will be included in all of the fasta file names and bowtie index file names
# windowSize: size of sliding window to include exon pairs in junction database, default 1000000

# example usage: ./createJunctionIndex.sh /home/linda/circularRNApipeline
#                                        /home/linda/index/pombe
#                                        ASM294v2_23_test
 
PIPELINE_DIR=$1
OUT_DIR=$2
FILE_ID=$3
if [ $# -ge 4 ]
then
  WINDOW=${4}  
else
  WINDOW=1000000  
fi

python makeJunctionsAndWriteFasta.py -w ${WINDOW} -e ${OUT_DIR}/exons -r ${OUT_DIR}/records -f ${OUT_DIR}/fastas -v

cat ${OUT_DIR}/fastas/*.fa > ${FILE_ID}.fa  # combine into single file

# split the junctions into files containing only reg, only rev, and only dup junctions
python limitFasta.py -s ${FILE_ID}.fa -o ${OUT_DIR}/fastas/ -t reg -p _junctions_reg
python limitFasta.py -s ${FILE_ID}.fa -o ${OUT_DIR}/fastas/ -t dup -p _junctions_dup
python limitFasta.py -s ${FILE_ID}.fa -o ${OUT_DIR}/fastas/ -t rev -p _junctions_rev

# combine rev and dup junctions into scrambled fasta file and put in the pipeline index
cat ${OUT_DIR}/fastas/${FILE_ID}_junctions_rev.fa ${OUT_DIR}/fastas/${FILE_ID}_junctions_dup.fa > ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled.fa

# put the regular junction fasta file into the pipeline index directory
mv ${OUT_DIR}/fastas/${FILE_ID}_junctions_reg.fa ${PIPELINE_DIR}/index

# remove the temp fastas created along the way 
rm ${FILE_ID}.fa
rm -r ${OUT_DIR}/fastas/

# create scrambled junction bowtie2 index and put into index directory
echo "bowtie2-build ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled.fa ${FILE_ID}_junctions_scrambled"
bowtie2-build ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled.fa ${FILE_ID}_junctions_scrambled

echo "mv ${FILE_ID}_junctions_scrambled.* ${PIPELINE_DIR}/index"
mv ${FILE_ID}_junctions_scrambled.* ${PIPELINE_DIR}/index

# create linear junction bowtie2 index and put into index directory 
echo "bowtie2-build ${PIPELINE_DIR}/index/${FILE_ID}_junctions_reg.fa ${FILE_ID}_junctions_reg"
bowtie2-build ${PIPELINE_DIR}/index/${FILE_ID}_junctions_reg.fa ${FILE_ID}_junctions_reg

echo "mv ${FILE_ID}_junctions_reg.* ${PIPELINE_DIR}/index"
mv ${FILE_ID}_junctions_reg.* ${PIPELINE_DIR}/index

