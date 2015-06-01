from Bio import SeqIO

EMPTY_SEQUENCE = "NOSEQUENCE" # to indicate junction sequence contains all Ns so do not include in fasta output

class junction: 
    def __init__(self, rec, chr, exon1, exon2):
        self.chr = chr  # string chromosome id  
        self.exons = [exon1, exon2] # SeqFeatures
        self.type = self.setType(exon1, exon2) # reg, rev, dup
        self.strand = exon1.strand # 1 or -1 
        self.seq = self.setSeq(rec, exon1, exon2) # 150 from each exon
        
    def __str__(self):
        return str(self.seq) + "\n"
    
    def setType(self, exon1, exon2):
        if exon1.location == exon2.location:
            return "dup"
        # exon1 comes before exon2
        elif exon1.location.start < exon2.location.start:
            if exon1.strand == 1:
                return "reg"
            else:
                return "rev"
        # exon2 comes before exon1
        else:
            if exon1.strand == 1:
                return "rev"
            else:
                return "reg"
            
    def setSeq(self, rec, exon1, exon2):
        leftSide = exon1.extract(rec.seq)
        rightSide = exon2.extract(rec.seq)
        
        # SeqRecord features automatically give you the reverse complement for - strand genes
        # but I want the + strand to simplify analysis downstream
        if self.strand == -1:
            rSide = leftSide.reverse_complement()
            lSide = rightSide.reverse_complement()
            rightSide = rSide
            leftSide = lSide
            
        padding = "N" * 150  # excessive padding so we know we always have 150 bases
        
        leftSide = padding + leftSide
        rightSide = rightSide + padding
        
        # take 150 from each side
        completeJuncSeq = leftSide[-150:] + rightSide[:150]
        
        # for less well-annotated genomes, sometimes we get an entire junction of Ns which is not useful
        if completeJuncSeq != padding*2:
            return completeJuncSeq
        
        return EMPTY_SEQUENCE
    
    # the SeqRecord features switch the start and end apparently for genes on the minus strand
    # So I need to switch them back since I want to print the coordinates of the junction
    def printHeader(self):
        msg = ">" + str(self.chr)
        msg += "|" + str(self.exons[0].qualifiers["gene_id"][0])
        if self.strand == 1:
            msg += ":" + str(self.exons[0].location.end)
        else:
            msg += ":" + str(self.exons[0].location.start)
        msg += "|" + str(self.exons[1].qualifiers["gene_id"][0])
        if self.strand == 1:
            msg += ":" + str(self.exons[1].location.start)
        else:
            msg += ":" + str(self.exons[1].location.end)
        msg += "|" + str(self.type)
        if self.strand == 1:
            msg += "|+"
        else:
            msg += "|-"
        msg += "\n"
        return msg

# new junction class for whole-genome fusion index
class junction_fusion: 
    def __init__(self, chr1, rec1, exon1, chr2, rec2, exon2):
        self.chr1 = chr1  # string chromosome id
        self.chr2 = chr2  # string chromosome id 
        self.exons = [exon1, exon2] # SeqFeatures
        self.strand1 = exon1.strand # 1 or -1
        self.strand2 = exon2.strand # 1 or -1 
        self.seq = self.setSeq(rec1, rec2) # 150 from each exon
        self.type = self.setType() # reg, rev, dup, fusion, strandcross
        
    def __str__(self):
        return str(self.seq) + "\n"
    
    def setType(self):
        # across chromosomes
        if self.chr1 != self.chr2:
            return "fusion"
        
        # same chromosome but different strands
        if self.strand1 != self.strand2:
            return "strandcross"
        
        # on same chromosome, figure out dup, rev, reg
        if (self.exons[0].location.start == self.exons[1].location.start
            and self.exons[0].location.end == self.exons[1].location.end):
            return "dup"
        
        # exon1 comes before exon2
        if self.exons[0].location.start < self.exons[1].location.start:
            if self.exons[0].strand == 1:
                return "reg"
            else:
                return "rev"
        # exon2 comes before exon1
        if self.exons[0].strand == 1:
            return "rev"
        else:
            return "reg"
            
    def setSeq(self, rec1, rec2):
        leftSide = self.exons[0].extract(rec1.seq)
        rightSide = self.exons[1].extract(rec2.seq)
        
        padding = "N" * 150  # excessive padding so we know we always have 150 bases
        
        leftSide = padding + leftSide
        rightSide = rightSide + padding
        
        # take 150 from each side
        return leftSide[-150:] + rightSide[:150]

    def printHeader(self):
        msg = ">" + str(self.chr1)
        msg += ":" + str(self.exons[0].qualifiers["gene_id"][0])
        msg += ":" + str(self.exons[0].location.end)
        if self.strand1 == 1:
            msg += ":+"
        else:
            msg += ":-"
        msg += "|" + str(self.chr2)
        msg += ":" + str(self.exons[1].qualifiers["gene_id"][0])
        msg += ":" + str(self.exons[1].location.start)
        if self.strand2 == 1:
            msg += ":+"
        else:
            msg += ":-"
        msg += "|" + str(self.type)
        
        msg += "\n"
        return msg

