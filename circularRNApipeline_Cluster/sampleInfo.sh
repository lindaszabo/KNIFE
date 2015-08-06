
if [ "$CLUSTER_TYPE" = "SLURM" ]
then
  TASK_VAR=${SLURM_ARRAY_TASK_ID}
elif [ "$CLUSTER_TYPE" = "SGE" ]
then
  TASK_VAR=${SGE_TASK_ID}
fi

READ_FILE=`awk 'FNR == '${TASK_VAR}' {print $1}' $TASK_DATA_FILE`
SAMPLE_ID=`awk 'FNR == '${TASK_VAR}' {print $2}' $TASK_DATA_FILE`