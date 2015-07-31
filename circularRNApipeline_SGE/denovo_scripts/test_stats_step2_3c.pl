$tempOutDir=$ARGV[0];  # output directory 
$screen=$ARGV[1]; # unique identifier for dataset

####################### COMPUTES DISGRIBUTION OF TESTSTATSO 
use List::Util qw(sum);
 use List::Util qw(first max maxstr min minstr reduce shuffle sum);
############# script should also store original reads and offsets of each alignment side 
open (F,$tempOutDir."hamming_".$screen);
open (O,">>".$tempOutDir."lindaOUTPUTEST".$screen);
# DBR4KXP1:210:D1PKYACXX:6:1101:8595:7784 2:N:0:AGTCAA17239454600094546000chr14-98ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZGGAAAATGGAAGGAAGTGAAGATTGACCCAAATATGTTTGCAAATTTCAGACAAAGGGAATCAAAGTTGTGGGAAAATGGZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZGGAAAATGGAAGGAAGTGAAGATTGACCCAAATATGTTTGCAAATTTCAGACAAAGGGAATCAAAGTTGTGGGAAAATGGAAGGAAGTGAAGATTGACZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
## hamming_ is an intermediate file with sequences padded with Zs so that all sequences with ends from the same bin pair should align "perfectly" if they are a contig. It is used to compute base by base hamming distance

## this loop iterates through the hamming file and orders all sequences from the same bin into a hash indexed by bin id.

while ($l=<F>){
($id,$start1,$start2,$bin1,$bin2,$chr,$strand,$length,$seq)=split(/\t/,$l);
$binid=$bin1."_".$bin2."_".$chr."_".$strand; # completely determines binid
$diff=$start1-$start2;

if (!(defined($info{$binid}))){
$info{$binid}=$diff;
$seqhash{$binid}=$seq;
$oseqhash{$binid}=$seq;
$readids{$binid}=$id;
$startlist{$binid}=$start1."\t".$start2;
$startlist1{$binid}=$start1;
}
if (defined($info{$binid})){
## put ids in according to startlinst1
@allstarts=split(/\t/,$startlist1{$binid});
$nc=scalar(@allstarts1);

### THESE DON"T HAVE ORDERING DEPENDENCE OR USE
$startlist{$binid} = $startlist{$binid}."\t". $start1."\t".$start2;
$info{$binid}= $info{$binid}."\t" .$diff;
$seqhash{$binid}=$seqhash{$binid}."\t".$seq;
$readids{$binid}=$id."\t".$readids{$binid};

### THESE DO HAVE ORDERING DEPENDENCE OR USE AND ONLY RECORD "record values" highest or lowest sequences
if ($start1>$allstarts[($nc-1)]){

$startlist1{$binid} = $startlist1{$binid}."\t". $start1;
$oseqhash{$binid}=$seqhash{$binid}."\t".$seq;
    print "$binid\t$start1 bigger thatn $allstarts[($nc -1)] list is $startlist1{$binid}\n";
    print " list is $oseqhash{$binid}\t $seq\n";
}
if ($start1<$allstarts[0]){
$startlist1{$binid} = $start1."\t".$startlist1{$binid};
$oseqhash{$binid}=$seq."\t".$seqhash{$binid};
    print "$binid\t$start1 smaller thatn $allstarts[0] list is $startlist1{$binid}\n";
 print "$seq and  list is $oseqhash{$binid}\n";
}
## now startlist starts with smallest start and ends with largest start; dittofor sequence
}

}

foreach $k (keys %info){
############ this loop addresses variance in alignment, alternatively, we can compute selfconsistency metricks with the zs
############ this can also make the largest contig if we want;
    $mypsum=0;

    @mys= split(/\t/,$seqhash{$k}); ## all sequences
$myn=scalar(@mys); # myn is the number of reads, duplicates counted
@mystarts=split(/\t/,$startlist{$k}); # two reads per sequence
my $numberdiff= numunique(@mystarts); # numberdiff could be between 2 and 2*numtot below
$numtot=scalar(@mystarts)/2; ## number of unique (start,stops), cannot exceed read totals.

## numtot is the number of distinct starting positions divided by 2

if ( $numberdiff> 0){

$counter=0;
$mysa=split(//,@mys);
## mys are all of the sequences
my @trans;
for ($counter=0; $counter<=length($mys[0]); $counter++){
## loop through each nt of the first sequence TO BE conservative, assuming we are only computing HD between sets of reads that overlap with the first sequence
    $trans[$counter]="";
foreach $sk (@mys){
   chomp($sk);

  $lsk=length($sk);
   $next=substr($sk,$counter,1);
   $trans[$counter]=$trans[$counter].$next;
## SO, at end of loop, trans[$counter] is a string that represents the vertical alignments of all sequences, and we desire to compute the consistency of these sequences. 
}

## these variables compute the number of instances of each string
@ca=    ($trans[$counter] =~ /A/g);
@cc=    ($trans[$counter] =~ /C/g);
@cg=    ($trans[$counter] =~ /G/g);
@ct=    ($trans[$counter] =~ /T/g);
@ns="";
    push @ns,scalar(@ca);
    push @ns,scalar(@cc);
    push @ns,scalar(@cg);
    push @ns,scalar(@ct);

## IGNORES Zs
## ns now has the number of instances of each letter in the alpabet

  shift @ns; ## removes the dummy initializer
   @nss= sort {$a <=> $b} @ns;
## SORTS the values so that (4,4,4,6) or (0,0,0,10) might be examples, the former is "random", the latter might represent something real.
## p is lack of fit, then print out p vector
$p=NA;
if (sum(@nss)>0){ ## ie there is something to count
$p = 1 - ( $nss[3] / sum(@nss) ) ; ## 1-  6/(4+4+4+6) or 10/10, ie sum the discrepancy at this position 

$mypsum=$p+$mypsum;
}
}

#>7__+_chr8_128903200_128902900_psum_0_n_7_ndiff_4_0

$newk=$k."_psum_$mypsum"."_nseq_$myn"."_numberdiff_".$numberdiff."_";

###########################

    @stuff=split(/_/,$k);
#    >7__+_chr8_128903200_128902900_psum_0_n_7_ndiff_4_0
#50_16500_chrM_-_psum_1.48934341258113_nseq_69_numberdiff_23
($bin1,$bin2,$chrom,$strand,$psd,$psum,$nsd,$myn,$mtnd,$numberdiff)=split(/_/,$newk);
 $header=">".$numtot."__".$strand."_".$chrom."_".$bin1."_".$bin2."_psum_".$mypsum."_nseq_".$myn."_numberdiff_".$numberdiff."_";

print ">"."_".$newk."\nINFO\t@stuff\t@nums\n$readids{$k}\n";
    @myseqs=split(/\t/,$seqhash{$k});
    @mystarts=split(/\t/,$startlist1{$k});
    $myn=scalar(@mystarts);
    $mys=scalar(@myseqs);

    $nletters=length($myseqs[0]);
    print "LETTERLENGTH $nletters\n";
print O "$header\n";
    for ($j=0;$j<$nletters;$j++){
	$s1=$myseqs[0];
	$s2=$myseqs[($mys-1)];
	$next=minstr(substr($s1,$j,1),substr($s2,$j,1));
##### since sequences are ordered according to start and end, $s1 and $s2 will correspond to the extrema. minstr will choose the nt corresponding to ACGT over Z.
	print O "$next";     }
# new line character INCLDUED
}
}

#}
sub mean {

$d=scalar(@s);
$sum=sum(@_)/@_;

}
sub numunique { 
    @mysort = sort {$a<=>$b} @_;
    $n=scalar(@mysort);
    $c=0;
    for ($i=1;$i<=$n;$i++){
	if ($mysort[$i] == $mysort[$i-1]){
	    $c=$c+1;
	}
    }
    return $n-$c;
}



sub variance {

@nums=@_;
$n=scalar(@_);
$mym=sum(@nums)/($n);
## ssq is x^2/n - 2x*n*mean +mean2

$mysd=  sum(@numsq) - $n* $mym*$mym;
    return $mysd;
}
