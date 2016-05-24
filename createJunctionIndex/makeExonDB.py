# Given a gtf genome annotation file and a genome fasta file, create and pickle
# dicts of exon objects. 1 file is created with a dict for each chromosome, and
# each dict contains 2 lists, 1 for forward strand exons on this chromosome
# and one for reverse strand exons. These can then be used to create junction objects
# and junction fasta files.

# Output directories will be created if they do not exist, and default values
# are provided. The genes, exons, and records directories will be created as
# subdirectories of the output directory.

# usage: python makeExonDB3.py -f genome.fasta -a genome.gtf
#                              [-o outDir] [-e exonDir] [-r recordDir] [-g geneDir] [-v]


### TODO:
### 1. try catching the deprecation warning from Biopython
### 2. allow single file option (as currently implemented) or to loop through all files in a directory
### 4. streamline code
### 5. add logging instead of printing to stdout?
### 6. create script to combine all genes so I can easily find data for any gene in 1 structure

import argparse
from Bio import SeqIO
from BCBio import GFF
import cPickle as pickle
import utils_os
import sys

# not currently using genes, but figured it may be useful to have these objects
# around for future analysis if necessary. Each gene contains a chromosome id
# and a SeqFeature object which contains all of the exons for the gene. Sequence
# data is not stored in this object, just positional info for each exon.
class gene:
    
    def __init__(self, chr, obj):
        self.chr = chr  # string chromosome id  
        self.feature = obj # SeqFeature for entire gene, contains exons as subfeatures
    
    def __str__(self):
        msg = "chr: " + str(self.chr)
        msg += " feature: " + str(self.feature)
        return msg

# tries getting primary gene name, if that doesn't exist tries getting secondary gene name.
# if neither primary or secondary gene name fields exist, an exception will be thrown
# which is caught and handled in the calling function
def getGeneName(use_feature):
    if args.name1 in use_feature.qualifiers:
        use_name = args.name1
    else: 
        use_name = args.name2
        
    return use_feature.qualifiers[use_name][0]
    
# the workhorse function.
# param chrSeqRecord: a SeqRecord for 1 chromosome
def parseChrFeatures(chrSeqRecord):
    chrExonsByDirection = [] # each element will be a (chrID, strand, dict with (exon start, stop):SeqFeature, 2 per chromosome (1 per strand)
    exonsPos = {}  # all positive strand exons
    exonsNeg = {}  # all neg strand exons
    geneExons = {}  # gene_name: RecFeature. don't actually use these right now, but may be useful at some point
    
    
    for ftr in chrSeqRecord.features: # this is 1 entire gene
        try:
            sftr = None # genes with a single exon are at the top level ftr
            for sftr in ftr.sub_features:
                # add exon to the exon dict for the current chr
                # keeping in separate dicts for now so we are only looking at 1 direction with 100kb
                if sftr.strand == 1:
                    exonsPos[(int(sftr.location.start), int(sftr.location.end))] = sftr
                else:
                    exonsNeg[(int(sftr.location.start), int(sftr.location.end))] = sftr
    
            # gene names are tied to the subrecords, and chromosome to SeqRecords, so need to track this info with the feature object 
            if sftr:
                geneExons[getGeneName(sftr)] = gene(chrSeqRecord.id, ftr)
            else:
                # this is a single exon gene so gene name is in the feature instead of subfeatures
                if args.verbose:
                    print "single exon gene " + str(getGeneName(ftr))
                geneExons[getGeneName(ftr)] = gene(chrSeqRecord.id, ftr)
                # and need to add the feature to the appropriate exonDict
                if ftr.strand == 1:
                    exonsPos[(int(ftr.location.start), int(ftr.location.end))] = ftr
                else:
                    exonsNeg[(int(ftr.location.start), int(ftr.location.end))] = ftr
        except Exception as e:
            print "Exception"
            print e
            print "error:", sys.exc_info()[0]
            print "parsing features for", ftr
            
    
    if len(exonsPos) > 0:
        chrExonsByDirection.append((chrSeqRecord.id, 1, exonsPos))  # chrSeqRecord.id is chr1 for example, 1 indicates pos strand
    if len(exonsNeg) > 0:
        chrExonsByDirection.append((chrSeqRecord.id, -1, exonsNeg))  # -1 indicates neg strand
        
    ########### save these objects for future use #############
    pk_geneExons = open(geneOutFullPath + '/geneExons_' + chrSeqRecord.id + '.pkl', 'wb')
    pickle.dump(geneExons, pk_geneExons)
    pk_geneExons.close()
    
    if args.verbose:
        print "num genes: " + str(len(geneExons))
    
    pk_exonsByDirection = open(exonOutFullPath + '/exonsByStrand_' + chrSeqRecord.id + '.pkl', 'wb')
    pickle.dump(chrExonsByDirection, pk_exonsByDirection)
    pk_exonsByDirection.close()
    
    if args.verbose:
        for elem in chrExonsByDirection:
            print "chromosome:" + str(elem[0])
            print "strand: " + str(elem[1])
            print "number of exons: " + str(len(elem[2]))
            
    pk_rec = open(recOutFullPath + '/rec_' + chrSeqRecord.id + '.pkl', 'wb')
    pickle.dump(chrSeqRecord, pk_rec)
    pk_rec.close()

if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastaFile', required=True, help='path to fasta file with chromosome sequences')
    parser.add_argument('-a', '--annotationFile', required=True, help='path to gff file with exon annotations')
    parser.add_argument('-o', '--outDir', help='directory to output files, will be created if it does not exist', default='output')
    parser.add_argument('-n1', '--name1', help='name of field in gtf to use for gene names', default='gene_name')
    parser.add_argument('-n2', '--name2', help='name of field in gtf to use for gene names if n1 does not exist', default='gene_id')
    parser.add_argument('-e', '--exonOutDir', help='directory to output exon pickle files, will be created within outDir if it does not exist', default='exons')
    parser.add_argument('-g', '--geneOutDir', help='directory to output gene pickle files, will be created within outDir if it does not exist', default='genes')
    parser.add_argument('-r', '--recOutDir', help='directory to output seq record pickle files, will be created within outDir if it does not exist', default='records')
    parser.add_argument('-v', '--verbose', help='print info about data obtained', action='store_true')
    args = parser.parse_args()
    
    exonOutFullPath = '/'.join([args.outDir, args.exonOutDir])
    geneOutFullPath = '/'.join([args.outDir, args.geneOutDir])
    recOutFullPath = '/'.join([args.outDir, args.recOutDir])
    
    # create output directories if necessary
    utils_os.createDirectory(args.outDir)
    utils_os.createDirectory(exonOutFullPath)
    utils_os.createDirectory(geneOutFullPath)
    utils_os.createDirectory(recOutFullPath)
    
    # read in the sequences 
    f_handle = open(args.fastaFile, "rU")
    f_dict = SeqIO.to_dict(SeqIO.parse(f_handle, "fasta"))
    f_handle.close()
    
    # only want exons for now. unfortunately can't limit by strand so have to do that after the fact
    # we can speed things up a bit by limiting to just the chromosomes we have sequences for
    limit_info = dict(
        gff_id = f_dict.keys(),
        gff_type = ["exon"])
    
    if args.verbose:
        print "data we are searching for in gff: " + str(limit_info)
    
    ########### read in the annotations, adding annotations to the sequences 
    a_handle = open(args.annotationFile, "rU")
    if args.verbose:
            print "opening annotation file", args.annotationFile
    for rec in GFF.parse(a_handle, base_dict=f_dict, limit_info=limit_info): # each rec is a SeqRecord (for 1 chromosome)
        if args.verbose:
            print "####### starting new chromosome: " + str(rec.id)
        # populate the data structures with genes and exons from this chromosome and pickle the objects
        parseChrFeatures(rec)
    a_handle.close()
    
    