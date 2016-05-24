#!/bin/sh

# This assumes that makeExonDB.py has already been called, creating a directory OUT_DIR that contains exons, genes, and records directories.
# This code now parses that exon info and creates junction fasta files for scrambled and linear junctions and bowtie2 indices and places these
# into the pipeline code so they can be used by specifying ${INDEX_MODE_ID} in the MODE parameter when calling findCircularRNA.sh.

# pipelineDirectory: path to the circularRNApipeline directory. Created files will be placed inside the index directory.
# outputDirectory: the output directory specified in the call to createExonDB.sh for this species
# fileIdentifier: should be short String to help distinguish it from other genomes. This will be included in all of the fasta file names and bowtie index file names
# primaryGeneName: field in gtf used to assign gene name, default: gene_name
# secondaryGeneName: field in gtf used to assign gene name if primary does not exist, default: gene_id
# windowSize: size of sliding window to include exon pairs in junction database, default: 1000000

# example usage: ./createJunctionIndex.sh /home/linda/circularRNApipeline
#                                        /home/linda/index/pombe
#                                        ASM294v2_23_test
#                                        1000000
#                                        gene_name
#                                        gene_id
 
PIPELINE_DIR=$1
OUT_DIR=$2
FILE_ID=$3

if [ $# -ge 4 ]
then
  WINDOW=${4}  
else
  WINDOW=1000000  
fi

if [ $# -ge 5 ]
then
  GENE_NAME_1=${5}
else
  GENE_NAME_1=gene_name
fi

if [ $# -ge 6 ]
then
  GENE_NAME_2=${6}
else
  GENE_NAME_2=gene_id
fi

python makeJunctionsAndWriteFasta.py -w ${WINDOW} -e ${OUT_DIR}/exons -r ${OUT_DIR}/records -f ${OUT_DIR}/fastas -n1 ${GENE_NAME_1} -n2 ${GENE_NAME_2} -v

# combine into single file, using xargs method to avoid argument list too long error in bash

echo "find ${OUT_DIR}/fastas/ -size 0 -delete"
find ${OUT_DIR}/fastas/ -size 0 -delete

echo "cat ${OUT_DIR}/fastas/*.fa > ${OUT_DIR}/${FILE_ID}.fa"
cat ${OUT_DIR}/fastas/*.fa > ${OUT_DIR}/${FILE_ID}.fa  # combine into single file

#ls ${OUT_DIR}/fastas | xargs -n 32 -P 8 cat >> ${OUT_DIR}/${FILE_ID}.fa  

# split the junctions into files containing only reg, only rev, and only dup junctions
python limitFasta.py -s ${OUT_DIR}/${FILE_ID}.fa -o ${OUT_DIR}/fastas/ -t reg -p _junctions_reg
python limitFasta.py -s ${OUT_DIR}/${FILE_ID}.fa -o ${OUT_DIR}/fastas/ -t dup -p _junctions_dup
python limitFasta.py -s ${OUT_DIR}/${FILE_ID}.fa -o ${OUT_DIR}/fastas/ -t rev -p _junctions_rev

# combine rev and dup junctions into scrambled fasta file and put in the pipeline index
echo "cat ${OUT_DIR}/fastas/${FILE_ID}_junctions_rev.fa ${OUT_DIR}/fastas/${FILE_ID}_junctions_dup.fa > ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled.fa"
cat ${OUT_DIR}/fastas/${FILE_ID}_junctions_rev.fa ${OUT_DIR}/fastas/${FILE_ID}_junctions_dup.fa > ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled.fa

# put the regular junction fasta file into the pipeline index directory
mv ${OUT_DIR}/fastas/${FILE_ID}_junctions_reg.fa ${PIPELINE_DIR}/index

# remove the temp fastas created along the way 
rm ${OUT_DIR}/${FILE_ID}.fa
rm -r ${OUT_DIR}/fastas/

# create scrambled junction bowtie2 index in index directory
echo "bowtie2-build ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled.fa ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled"
bowtie2-build ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled.fa ${PIPELINE_DIR}/index/${FILE_ID}_junctions_scrambled

# create linear junction bowtie2 index in index directory 
echo "bowtie2-build ${PIPELINE_DIR}/index/${FILE_ID}_junctions_reg.fa ${PIPELINE_DIR}/index/${FILE_ID}_junctions_reg"
bowtie2-build ${PIPELINE_DIR}/index/${FILE_ID}_junctions_reg.fa ${PIPELINE_DIR}/index/${FILE_ID}_junctions_reg

