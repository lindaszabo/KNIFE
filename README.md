# KNIFE
Known and Novel IsoForm Explorer. Statistically based splicing detection for circular and linear isoforms

# Overview
This statistical algorithm increases the sensitivity and specificy of circularRNA detection by quantifying circular and linear RNA splicing events at both annotated and un-annotated exon boundaries. A version of the code that runs on a single Linux machine is provided in circularRNApipeline_Standalone. circularRNApipeline_SLURM contains code that provides the same results, but includes wrapper scripts to run on a cluster using the SLURM scheduler. This is the functional code we ran on our cluster, but each scheduler deployment is configured slightly differently. Therefore, this is not meant to be used as is like the standalone version of the code is, but rather to be used as a template that you can modify to work with your scheduler. See README files within each of these directories for detailed usage instructions.

As for all RNA-Seq analysis tools, our algorithm may perform poorly on low quality data. Therefore, we strongly recommend using our algorithm only on reads that have been processed by a tool such as cutAdapt to trim poor quality ends and adapters as predictions will be adversely affected by low-quality reads.

# Using Available Genomes
Annotated junction indices are available for Human (hg19), Mouse (mm10), Rat (rn5) and Drosophila (dm3). For convenience, we have also packaged up all of the the transcriptome, genome, and ribosomal index files and fasta files for each of these genomes named as required for use with our scripts. The Bowtie2 tar must be downloaded and uncompressed (use tar zxvf filename.tar.gz) into the circularRNApipeline/index directory. The gtf file must be downloaded and uncompressed (use gunzip filename.gtf.gz) into the circularRNApipeline/denovo_scripts directory. Files are available here: https://mega.nz/#F!59UQxJ4a!GLaLflWtyOrNTEdyBKqLEQ

You will also need to obtain Bowtie1 genome index files for the species (available from iGenomes at http://support.illumina.com/sequencing/sequencing_software/igenome.html), rename the files so that everything before the first "." is the same as the Bowtie2 index files we provided for the species, and place these index files in the  circularRNApipeline/denovo_scripts/index directory. 

Code and instructions for creating a new index are provided in createJunctionIndex.

# Software Requirements
Bowtie 2.2.2; Bowtie 0.12.7; samtools; perl; python (tested with 2.7.5) and the numpy and scipy libraries; R (tested with 3.0.2) and the data.table package

# Cite
Szabo L, Morey R, Palpant NJ, Wang PL, Afari N, Jiang C, Parast MM, Murry CE, Laurent LC, Salzman J. Tissue-specific induction of circular RNA during human fetal development revealed by statistically based splicing detection.
