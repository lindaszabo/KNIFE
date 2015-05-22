# KNIFE
Known and Novel IsoForm Explorer. Statistically based splicing detection for circular and linear isoforms

# Overview
This statistical algorithm increases the sensitivity and specificy of circularRNA detection from RNA-Seq data by quantifying circular and linear RNA splicing events at both annotated and un-annotated exon boundaries. It can analyze single-end or paired-end reads stored in plain-text fastq or gzipped fastq files. A version of the code that runs on a single Linux machine is provided in circularRNApipeline_Standalone. circularRNApipeline_SLURM contains code that provides the same results, but includes wrapper scripts to run on a cluster using the SLURM scheduler. This is the functional code we ran on our cluster, but each scheduler deployment is configured slightly differently. Therefore, this is not meant to be used as is like the standalone version of the code is, but rather to be used as a template that you can modify to work with your scheduler. Parameters and output files are described below, and README files within each of these code directories contain detailed instructions for exactly how to run the code.

As for all RNA-Seq analysis tools, our algorithm may perform poorly on low quality data. Therefore, we strongly recommend using our algorithm only on reads that have been processed by a tool such as cutAdapt to trim poor quality ends and adapters as predictions will be adversely affected by low-quality reads.

# Software Requirements

- Bowtie 2.2.2 
- Bowtie 0.12.7
- perl 
- python (tested with 2.7.5) and the numpy and scipy libraries
- R (tested with 3.0.2) and the data.table package
- samtools (if you want alignment files output as bam instead of sam files)

# Using Available Genomes
Annotated junction indices are available for Human (hg19), Mouse (mm10), Rat (rn5) and Drosophila (dm3). For convenience, we have also packaged up all of the the transcriptome, genome, and ribosomal index files and fasta files for each of these genomes named as required for use with our scripts. The Bowtie2 tar must be downloaded and uncompressed (use tar zxvf filename.tar.gz) into the circularRNApipeline/index directory. The gtf file must be downloaded and uncompressed (use gunzip filename.gtf.gz) into the circularRNApipeline/denovo_scripts directory. Files are available here: https://mega.nz/#F!59UQxJ4a!GLaLflWtyOrNTEdyBKqLEQ

You will also need to obtain Bowtie1 genome index files for the species (available from iGenomes at http://support.illumina.com/sequencing/sequencing_software/igenome.html), rename the files so that everything before the first "." is the same as the Bowtie2 index files we provided for the species, and place these index files in the  circularRNApipeline/denovo_scripts/index directory. 

# Using Other Genomes
Code and instructions for creating a new index are provided in createJunctionIndex.

# Script Parameters
Running the algorithm is slightly different in the Standalone or SLURM versions (see README files within those directories for exact usage information), but the parameters are the same. There are a total of 10 possible parameters, the first 5 are required and defaults for the rest are described here:

- read_directory: full path to directory containing fastq files for alignment. Paired-end reads (PE) must have read1 and read2 in separate files. The file names for a given sample must start with an identical string that identifies the sample and then have either _1, _R1, _2, or _R2 identifying them as read1 or read2. Any part of the file name after this will be ignored. Reads
may be in gzipped files or plain text fastq files. For simplicity, single end read files must follow naming requirements for read1 shown below (must have _1 or _R1 in file name). Valid naming examples are:

  ```
  SAMPLENAME_1.fastq, SAMPLENAME_2.fastq
  SAMPLENAME_1.fq.gz, SAMPLENAME_2.fq.gz
  SAMPLENAME_R1.fq, SAMPLENAME_R2.fq
  ```

- read_id_style: complete|appended (use complete for single end)

  ```
  complete: read ids in the fastq files for read 1 and read 2 are identical for mates
  appended: the last character of the read id is different in the 2 read files for mates
  ```
  
- alignment_parent_directory: full path to directory where the dataset analysis output will be stored. This directory must already exist, and a directory named dataset_name (see below) will be created under this directory for all output files.

- dataset_name: string identifier for this dataset. A folder of this name will be created under alignment_parent_directory (see above) and all output for this run will be stored in this directory.

- junction_overlap: minimum number of bases in the read which must be on each side of the junction to consider that the read is overlapping the junction. Values that have empirically worked well are 8 for paired-end (PE) reads of length < 70, 13 for longer
PE, and 10 for single-end (SE) reads of length < 70, 15 for longer SE reads.

- mode: the mode that the script will be run in. This affects which portions of the code are executed, output file types, etc. Portions of the mode string can be specified in any order. For example, mouse_phred64 is the same as phred64_mouse.

  ```
  sam (default): output alignments as sam files
  bam: output alignments as bam files. takes longer to run, but reduces disk space required to store 
       alignments to about 1/6 of the sam files.
  bam_sort: output alignments as sorted bam files. takes a lot longer to run, and reduces disk space 
            to about 2/3 the size of the unsorted bam files.
  rat|mouse|fly: need to append this string to the mode for the respective organism.
  analysis: alignment sam or bam files have already been generated, so just go straight to analysis 
            instead of running Bowtie2
  phred64: need to append this string to the mode to indicate that quality scores are in phred64, 
           else phred33 is assumed
  ```
  
- report_directory_name: Default circReads. This is the name of the directory created under alignment_parent_directory/dataset_name where the report files listing circular or linear RNA junction read counts and statistical scores will be output.

- ntrim: Default 50. For creating de novo index, what is the number of bases to trim from each end of the read. Suggested value is 2n/3 where n is the length of the raw (untrimmed) reads in the dataset.

- denovoCircMode: 0|1, default 1. Should we report only circular junctions in the denovo analysis? changing to 0 means linear denovo junctions are included too. NOTE: denovoCircMode=0 and 1 for the same dataset cannot be run at the same time as there are some dataset-specific files that will get overwritten. If you want to run denovoCircMode=0 and 1 on the same dataset, use a different output directory for the 2 unaligned runs so your files do not get overwritten. If you also want to save the sam files for both unaligned mode runs, you will need to copy orig/denovo and orig/still_unaligned before the 2nd run because these files will be overwritten.

- junction_id_suffix: Default none. During the analysis step, the ids of all scrambled junction-overlapping reads are output to alignment_parent_directory/dataset_name/orig/ids/junction and the ids of all canonical-order junction-overlapping reads are output to alignment_parent_directory/dataset_name/orig/ids/reg. Subsequent steps then read these files to report on reads aligning to each junction. If you want to try requiring different amounts of junction overlap in order to compare results, you can specify a string be added to the id directory name so that the original results are not overwritten. If you use this option, also make sure that you use a different report_directory_name so that your reports for this run are output to a different directory instead of overwriting those in the default directory.


# Output


# Cite
Szabo L, Morey R, Palpant NJ, Wang PL, Afari N, Jiang C, Parast MM, Murry CE, Laurent LC, Salzman J. Tissue-specific induction of circular RNA during human fetal development revealed by statistically based splicing detection.
