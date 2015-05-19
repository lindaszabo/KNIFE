# KNIFE
Known and Novel IsoForm Explorer. Statistically based splicing detection for circular and linear isoforms

# Overview
This statistical algorithm increases the sensitivity and specificy of circularRNA detection by quantifying circular and linear RNA splicing events at both annotated and un-annotated exon boundaries. A version of the code that runs on a single Linux machine is provided in circularRNApipeline_Standalone. circularRNApipeline_SLURM contains code that provides the same results, but includes wrapper scripts to run on a cluster using the SLURM scheduler. This is the functional code we ran on our cluster, but each scheduler deployment is configured slightly differently. Therefore, this is not meant to be used as is like the standalone version of the code is, but rather to be used as a template that you can modify to work with your scheduler. See README files within each of these directories for detailed usage instructions.

As for all RNA-Seq analysis tools, our algorithm may perform poorly on low quality data. Therefore, we strongly recommend using our algorithm only on reads that have been processed by a tool such as cutAdapt to trim poor quality ends and adapters as predictions will be adversely affected by low-quality reads.

# Software Requirements
Bowtie 2.2.2; Bowtie 0.12.7; samtools; perl; python (tested with 2.7.5) and the numpy and scipy libraries; R (tested with 3.0.2) and the data.table package

# Cite
Szabo L, Morey R, Palpant NJ, Wang PL, Afari N, Jiang C, Parast MM, Murry CE, Laurent LC, Salzman J. Tissue-specific induction of circular RNA during human fetal development revealed by statistically based splicing detection.
