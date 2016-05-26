# KNIFE
Known and Novel IsoForm Explorer. Statistically based splicing detection for circular and linear isoforms from RNA-Seq data.

Note: as of release v1.3, circular and linear junctions are quantified based on both R1 and R2 junctional alignments. Prior versions considered only R1 junctional alignments. If both R1 and R2 aligned to the same junction as can be the case for small circles, that read is only counted once in the reported totals.

Note: It has been brought to our attention that there was a discrepancy in the SRA annotations in the NCBI short read archive for the HeLa samples we analyzed for Figure 3 in our manuscript and their true identity. Please see misAnnotationUpdate.zip available in Downloads under release v1.2 for updated analysis.  

# Overview
This statistical algorithm increases the sensitivity and specificy of circularRNA detection from RNA-Seq data by quantifying circular and linear RNA splicing events at both annotated and un-annotated exon boundaries. It can analyze single-end or paired-end reads stored in plain-text fastq or gzipped fastq files. A version of the code that runs on a single Linux machine is provided in circularRNApipeline_Standalone. circularRNApipeline_Cluster contains code that provides the same results, but includes wrapper scripts to run on a cluster using either a SLURM or Sun Grid Engine scheduler. Parameters and output files are described below, and README files within each of these code directories contain detailed instructions for exactly how to run the code. 

As for all RNA-Seq analysis tools, our algorithm may perform poorly on low quality data. Therefore, we strongly recommend using our algorithm only on reads that have been processed by a tool such as cutAdapt to trim poor quality ends and adapters as predictions will be adversely affected by low-quality reads.

# Installation

See INSTALL file for instructions on installing KNIFE.

# Software Requirements

- Bowtie2 2.2.1 or higher (tested with 2.2.1)
- Bowtie 0.12.7 or higher (tested with 0.12.7)
- perl 
- python (tested with 2.7.5) and the numpy and scipy libraries
- R (tested with 3.0.2) and the data.table package
- samtools (if you want alignment files output as bam instead of sam files)

# Using Available Genomes
Annotated junction indices are available for Human (hg19), Mouse (mm10), Rat (rn5) and Drosophila (dm3). We have also packaged up all of the the transcriptome, genome, and ribosomal index, fasta, and gtf files for each of these genomes named as required for use with our scripts. See INSTALL file for directions on installing these genomes.

# Using Other Genomes
Code and instructions for creating a new index are provided in createJunctionIndex.

# Script Parameters
Running the algorithm is slightly different in the Standalone or Cluster implementaions (see README files within those directories for exact usage information), but the parameters are the same. There are a total of 10 possible parameters, the first 5 are required and defaults for the rest are described here:

- read_directory: absolute path to directory containing fastq files for alignment. Paired-end reads (PE) must have read1 and read2 in separate files. The file names for a given sample must start with an identical string that identifies the sample and then have either _1, _R1, _2, or _R2 identifying them as read1 or read2. Any part of the file name after this will be ignored. Reads may be in gzipped files or plain text fastq files, but the file extension must be .fq, .fq.gz, .fastq, or .fastq.gz. Read file names must not contain the string 'unaligned' as this will interfere with logic around identifying de novo junctions. For simplicity, single end read files must follow naming requirements for read1 shown below (must have _1 or _R1 in file name). Valid naming examples are:

  ```
  SAMPLENAME_1.fastq, SAMPLENAME_2.fastq
  SAMPLENAME_1.fq.gz, SAMPLENAME_2.fq.gz
  SAMPLENAME_R1.fq, SAMPLENAME_R2.fq
  ```

- read_id_style: complete|appended (use complete for single end). 
```
Note that the read id is only what is up to the first space on the sequence identifier line. 
Everything after that is an optional description and is ignored. 
So for example, if the sequence identifier lines in your 2 files are: 
 @SRR445016.1 UNC15-SN850_0140:1:1101:1247:2164/1
 @SRR445016.1 UNC15-SN850_0140:1:1101:1247:2164/2
 then the read id is SRR445016.1 and the style "complete" should be used
 
  complete: read ids in the fastq files for read 1 and read 2 are identical for mates
  appended: the last character of the read id is different in the 2 read files for mates 
            (ex: @27847_787_327_212/1, @27847_787_327_212/2)
  ```
  
- alignment_parent_directory: absolute path to directory where the dataset analysis output and log files will be stored. This directory must already exist, and a directory named dataset_name (see below) will be created under this directory for all output files.

- dataset_name: string identifier for this dataset. A folder of this name will be created under alignment_parent_directory (see above) and all output for this run will be stored in this directory.

- junction_overlap: minimum number of bases in the read which must be on each side of the junction to consider that the read is overlapping the junction. Values that have empirically worked well are 8 for paired-end (PE) reads of length < 70, 13 for longer PE, and 10 for single-end (SE) reads of length < 70, 15 for longer SE reads.

- mode: the mode that the script will be run in. This affects which portions of the code are executed, output file types, etc. Portions of the mode string can be specified in any order. For example, mouse_phred64 is the same as phred64_mouse.

  ```
  sam (default): output alignments as sam files
  bam: output alignments as bam files. takes longer to run, but reduces disk space required to store 
       alignments to about 1/6 of the sam files.
  bam_sort: output alignments as sorted bam files. takes a lot longer to run, and reduces disk space 
            to about 2/3 the size of the unsorted bam files.
  rat|mouse|fly: need to append this string to the mode for the respective organism, otherwise will default to human.
  analysis: alignment sam or bam files have already been generated, so just go straight to analysis step 
            instead of running Bowtie2
  phred64: need to append this string to the mode to indicate that quality scores are in phred64, 
           else phred33 is assumed
  ```
  
- report_directory_name: Default circReads. This is the name of the directory created under alignment_parent_directory/dataset_name where the report files listing circular or linear RNA junction read counts and statistical scores will be output.

- ntrim: Default 50. When creating the de novo index, this is the number of bases to trim from each end of an unaligned read to try to align independently to the genome. Suggested value is 2n/3 where n is the length of the raw (untrimmed) reads in the dataset.

- denovoCircMode: 0|1, default 1. Should we report only circular junctions in the denovo analysis? Changing to 0 means linear denovo junctions are included too. NOTE: denovoCircMode=0 and 1 for the same dataset cannot be run at the same time as there are some dataset-specific files that will get overwritten. If you want to run denovoCircMode=0 and 1 on the same dataset, use a different output directory for the 2 unaligned runs so your files do not get overwritten. If you also want to save the sam files for both unaligned mode runs, you will need to copy orig/denovo and orig/still_unaligned before the 2nd run because these files will be overwritten.

- junction_id_suffix: Default none. During the analysis step, the ids of all scrambled junction-overlapping reads are output to alignment_parent_directory/dataset_name/orig/ids/junction and the ids of all canonical-order junction-overlapping reads are output to alignment_parent_directory/dataset_name/orig/ids/reg. Subsequent steps then read these files to report on reads aligning to each junction. If you want to try requiring different amounts of junction overlap in order to compare results, you can specify a string be added to the id directory name so that the original results are not overwritten. If you use this option, also make sure that you use a different report_directory_name so that your reports for this run are output to a different directory instead of overwriting those in the default directory.


# Output
All output files can be found under [alignment_parent_directory]/[dataset_name] as specified when the script is called. Within this directory the following subdirectories will be created:

1. circReads (or the directory will be named as specified by the [report_directory_name] parameter). The primary output files you will be interested in looking at are in the following subdirectories. 
    1. reports: read count and p-value per junction using naive method. 2 files created per sample, 1 for annotated junctions (linear and circular) and the other for de novo junctions. For single end reads and de novo junctions from either single end or paired end data, these are the output files of interest as GLM reports are for annotated junctions using paired end data only. You will want to select a threshold on the p-value for which of these junctions are considered true positive circles. For the publication, we considered all junctions with a p-value of 0.9 or higher and a decoy/circ read count ratio of 0.1 or lower.
      - junction: chr|gene1_symbol:splice_position|gene2_symbol:splice_position|junction_type|strand
        - junction types are reg (linear), rev (circle formed from 2 or more exons), or dup (circle formed from single exon)
      - linear: number of reads where read1 aligned to this linear junction and read2 was consistent with presumed splice event, or just number of aligned reads to this linear junction for SE reads
      - anomaly: number of reads where read2 was inconsistent with read1 alignment to this linear junction
      - unmapped: number of reads where read1 aligned to this junction and read2 did not map to any index
      - multimapped: N/A
      - circ: number of reads where read1 aligned to this circular junction and read2 was consistent with presumed splice event, or just number of aligned reads to this circular junction for SE reads
      - decoy: number of reads where read2 was inconsistent with read1 alignment to this circular junction
      - pvalue: naive method p-value for this junction based on all aligned reads (higher = more likely true positive). You will want to select a threshold on the p-value for which of these junctions are considered true positive circles. 
      - scores: (read1, read2) Bowtie2 alignment scores for each read aligning to this junction, or scores at each 10th percentile for junctions with more than 10 reads
    2. glmReports: read count and posterior probability per junction using GLM (only for PE reads, annotation-dependent junctions). 2 files created per sample, 1 with circular splice junctions and the other with linear splice junctions. You will want to select a threshold on the posterior probability for which of these junctions are considered true positive circles. For the publication, we considered all junctions with a posterior probability of 0.9 or higher.
      - junction: chr|gene1_symbol:splice_position|gene2_symbol:splice_position|junction_type|strand
        - junction types are reg (linear), rev (circle formed from 2 or more exons), or dup (circle formed from single exon)
      - numReads: number of reads where read1 aligned to this junction and read2 was consistent with presumed splice event
      - p_predicted: posterior probability that the junction is a true junction (higher = more likely true positive). 
      - p_value: p-value for the posterior probability to control for the effect of total junctional counts on posterior probability (linear reports only)
    3. glmModels: RData files containing the model used to generate the glmReports 
    4.  ids: alignment and category assignment per read
    5.  combinedReports: aggregated read counts per junction considering cases where R1 or R2 aligned to the junction. If R1 and R2 both aligned to the same junction, they are only counted once. 
    
    circleJuncProbs.txt and linearJuncProbs.txt are the summaries of the GLM reports. 
      - The first 3 columns contain values from glmReports from the original run
      - swapped_count: number of reads where R2 aligned to this junction and R1 was consistent with presumed splice event
      - swapped_posterior: posterior probability that the junction is a true junction based only on R2 junctional alignments
      - total_reads: number of reads where R1 or R2 aligned to this junction. This may not be the sum of orig_count and swapped_count because R1 and R2 may have aligned to the same junction and should not be double-counted
    
    naiveunaligned.txt is the summary of alignments to the denovo index. 
      - orig_circOrLinear: number of reads where read1 aligned to this junction and read2 was consistent with presumed splice 
      - orig_decoyOrAnom: number of reads where read2 was inconsistent with read1 alignment to this junction
      - orig_unmapped: number of reads where read1 aligned to this junction and read2 did not map to any index
      - orig_pval: naive method p-value for this junction based on all R1s aligned to the junction (higher = more likely true positive).
      - swapped_circOrLinear: number of reads where read2 aligned to this junction and read1 was consistent with presumed splice 
      - swapped_decoyOrAnom: number of reads where read1 was inconsistent with read2 alignment to this junction
      - swapped_unmapped: number of reads where read2 aligned to this junction and read1 did not map to any index
      - swapped_pval: naive method p-value for this junction based on all R2s aligned to the junction
      - total_reads: number of reads where R1 or R2 aligned to this junction. This may not be the sum of orig_circOrLinear and swapped_circOrLinear because R1 and R2 may have aligned to the same junction and should not be double-counted
2. sampleStats: Contains 2 txt files with high-level alignment statistics per sample (read1 and read2 reported separately).
   * SampleAlignStats.txt: useful for evaluating how well the library prep worked, for example ribosomal depletion. Number of reads are reported, with fraction of total reads listed in ()  
     - READS: number of reads in original fastq file
     - UNMAPPED: number of reads that did not align to any of the junction, genome, transcriptome, or ribosomal indices
     - GENOME: number of reads aligning to the genome
     - G_STRAND: percentage of GENOME reads aligning to forward strand and percentage aligning to reverse strand
     - JUNC: number of reads aligning to the scrambled or linear junction index and overlapping the junction by required amount
     - J_STRAND: percentage of JUNC reads aligning to forward strand and percentage aligning to reverse strand
     - RIBO: number of reads aligning to the ribosomal index
     - R_STRAND: percentage of RIBO reads aligning to forward strand and percentage aligning to reverse strand
     - 28S, 18S, 5.8S, 5SDNA, 5SrRNA: percentage of RIBO aligning to each of these ribosomal subunits (for human samples only)
     - HBB: number of reads aligning to HBB genomic location (per hg19 annotation)
   * SampleCircStats.txt: useful for comparing circular and linear ratios per sample
      - CIRC_STRONG: number of reads that aligned to a circular junction that has a p-value >= 0.9 using the naive method (very high confidence of true circle)
     - CIRC_ARTIFACT: number of reads that aligned to a circular junction that has a p-value < 0.9 using the naive method 
     - DECOY: number of reads where read2 did not align within the circle defined by read1 alignment
     - LINEAR_STRONG: number of reads that aligned to a linear junction that has a p-value >= 0.9 using the naive method (very high confidence of true linear splicing)
     - LINEAR_ARTIFACT: number of reads that aligned to a linear junction that has a p-value < 0.9 using the naive method
     - ANOMALY: number of reads where read2 did not support a linear transcript that includes the read1 junction alignment
     - UNMAPPED: number of reads where read1 aligned to a linear or scrambled junction but read2 did not map to any index
     - TOTAL: sum of all previous columns, represents total number of reads mapped to junction but not to genome or ribosome
     - CIRC_FRACTION: CIRC_STRONG / TOTAL
     - LINEAR_FRACTION: LINEAR_STRONG / TOTAL
     - CIRC / LINEAR: CIRC_FRACTION / LINEAR_FRACTION
3. orig: contains all sam/bam files output and information used to assign reads to categories. In general there is no reason to dig into these files since the results, including the ids of reads that aligned to each junction, are output in report files under circReads as described above, but sometimes it is useful to dig back through if you want to trace what happened to a particular read.
  1. genome: sam/bam files containing Bowtie2 alignments to the genome index
  2. junction: sam/bam files containing Bowtie2 alignments to the scrambled junction index
  3. reg: sam/bam files containing Bowtie2 alignments to the linear junction index
  4. ribo: sam/bam files containing Bowtie2 alignments to the ribosomal index
  5. unaligned: fastq and fasta files for all reads that did not align to any index
    1. forDenovoIndex: fastq files containing subset of the unaligned reads that are long enough to be used for creating the denovo junction index
  6. denovo: sam/bam files containing Bowtie2 alignments to the de novo junction index
  7. still_unaligned: fastq files containing the subset of the unaligned reads that did not align to the denovo index either 
  8. ids: text files containing the ids of reads that aligned to each index, location of alignment, and any other relevant data from the sam/bam files used in subsequent analysis. The reads reported in the junction and reg subdirectories are only those that overlapped the junction by user-specified amount. In juncNonGR and denovoNonGR, the reported read ids are the subset of reads that overlapped a junction and did not align to the genome or ribosomal index.
4. logs (Cluster version only, for Standalone stderr and stdout need to be redirected to log file)
  1. align: logs from all bowtie2 alignment calls
  2. analysis: logs from processing the Bowtie2 output to generate the naive-method and glm reports
  3. denovo_index: fasta and bowtie2 index files generated for the de novo index 
  4. denovo_script_out: debugging output generated during creation of de novo index.
  5. sampleStats: logs from generating sample statistics files

If you have generated the combinedReports directory above (this occurs automatically in the circularRNApipeline_Standalone code and requires an additional call in the circularRNApipeline_Cluster code - see circularRNApipeline_Cluster README for details), an additional output directory is created as [alignment_parent_directory]/[dataset_name]Swapped. The output under the "circReads" subdirectories here contains equivalent information as described above for read1, but generated from treating R2 as R1. The output detailed above for the "orig" subdirectories and "sampleStats" is not present because Bowtie2 sam outputs from the initial run are re-used for this step.

# Selecting High Quality Circular/Linear Junctions
True circular and linear junctions are a subset of the reported junctions. All identified junctions are reported along with a posterior probability (GLM reports) or p-value (naive reports) which is used to flag false positives. Instead of selecting a threshold on read counts as you would for other circular RNA algorithms, you will select a threshold on these statistical scores and filter the reports. For the publication, we generally considered junctions with a score of 0.9 or higher to be true positive circular RNA (see publication for details). You can filter the results on read count as well if you want to limit your analysis to highly expressed RNA.

# Verifying Installation
Pre-trimmed fastq files for 1 human sample are provided with the release of v1.1. See the README under testData for instructions to run KNIFE on this sample. Expected results are provided in testData/testOutput.

# Contact
This code was developed and is maintained by Linda Szabo, lszabo at stanford.edu. Please understand that I am continuing to develop new algorithms, but I will do my best to respond in a timely manner. 

# Cite
Szabo L, Morey R, Palpant NJ, Wang PL, Afari N, Jiang C, Parast MM, Murry CE, Laurent LC, Salzman J. Statistically based splicing detection reveals neural enrichment and tissue-specific induction of circular RNA during human fetal development. Genome Biology. 2015, 16:126.  http://www.ncbi.nlm.nih.gov/pubmed/26076956

