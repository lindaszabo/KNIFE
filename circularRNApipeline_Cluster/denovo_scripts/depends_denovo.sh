#if [ "$CLUSTER_TYPE" = "SLURM" ]
#then
#  myBowtie=`which bowtie`
#  echo "my bowtie: $myBowtie"
#elif [ "$CLUSTER_TYPE" = "SGE" ]
#then
#  module load bowtie/1.1.1  # denovo uses bowtie1 instead of bowtie2 for selecting junctions to be included in denovo index
#  module load bowtie/2.2.4 # must be 2.2.2 or higher
#fi
