
#if [ "$CLUSTER_TYPE" = "SLURM" ]
#then
#  module load R  # must be version with data.table installed
#  module load python/2.7.5  # must be 2.7+, will not work with python 3
#elif [ "$CLUSTER_TYPE" = "SGE" ]
#then
#  module load r/3.2.0  # must be version with data.table installed
#  module load python/2.7  # must be 2.7+, will not work with python 3
#  module load bowtie/2.2.4 # must be 2.2.2 or higher
#fi
