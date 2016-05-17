#!/bin/bash

# Called after qualityStatsSingleSample.sh has completed to generate an inidividual report
# per sample. Consolidates all into 2 reports for the entire dataset and deletes the individual
# report files.

STATS_DIR=$1

ALIGN_FILE=${STATS_DIR}/SampleAlignStats.txt
CIRC_FILE=${STATS_DIR}/SampleCircStats.txt

echo -e "ID\tREADS\tUNALIGNED\tGENOME\tG_strand +,-\tREG\tReg_strand +,-\tJUNC\tJ_strand +,-\tRIBO\tR_strand +,-\t28S\t18S\t5.8S\t5SDNA\t5SrRNA\tHBB" > ${ALIGN_FILE}
cat ${STATS_DIR}/*AlignmentStats.txt >> ${ALIGN_FILE}

echo -e "ID\tCIRC_STRONG\tCIRC_ARTIFACT\tDECOY\tLINEAR_STRONG\tLINEAR_ARTIFACT\tANOMALY\tUNMAPPED\tTOTAL\tCIRC_FRACTION\tLINEAR_FRACTION\tCIRC/LINEAR" > ${CIRC_FILE}
cat ${STATS_DIR}/*CircularStats.txt >> ${CIRC_FILE}

rm ${STATS_DIR}/*AlignmentStats.txt
rm ${STATS_DIR}/*CircularStats.txt




