READ_FILE=`awk 'FNR == '${SGE_TASK_ID}' {print $1}' $TASK_DATA_FILE`
SAMPLE_ID=`awk 'FNR == '${SGE_TASK_ID}' {print $2}' $TASK_DATA_FILE`