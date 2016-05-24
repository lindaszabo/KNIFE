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
    
    # tries getting primary gene name, if that doesn't exist tries getting secondary gene name.
    # if neither primary or secondary gene name fields exist, an exception will be thrown
    # which is caught and handled in the calling function
    def getGeneName(self, use_feature, use_n1, use_n2):
        if use_n1 in use_feature.qualifiers:
            use_name = use_n1
        else: 
            use_name = use_n2
            
        return use_feature.qualifiers[use_name][0]
    
    # the SeqRecord features switch the start and end apparently for genes on the minus strand
    # So I need to switch them back since I want to print the coordinates of the junction
    # n1 is primary name field to use, n2 is backup name field to use 
    def printHeader(self, n1, n2):
        msg = ">" + str(self.chr)
        msg += "|" + str(self.getGeneName(self.exons[0], n1, n2))
        if self.strand == 1:
            msg += ":" + str(self.exons[0].location.end)
        else:
            msg += ":" + str(self.exons[0].location.start + 1)
        msg += "|" + str(self.getGeneName(self.exons[1], n1, n2))
        if self.strand == 1:
            msg += ":" + str(self.exons[1].location.start + 1)
        else:
            msg += ":" + str(self.exons[1].location.end )
        msg += "|" + str(self.type)
        if self.strand == 1:
            msg += "|+"
        else:
            msg += "|-"
        msg += "\n"
        return msg

