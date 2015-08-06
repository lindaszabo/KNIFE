use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Util qw(sum);
 use POSIX;

################################################################ JUST ANNOTATION
         $binsize=50;   
$boundary=0; # boundary is the nt increase on the bin which I have set to 0 but should be taken into account when assigning decoy status
open (FL,$ARGV[0]); ## fasta file from processing unaligned reads
open (F,$ARGV[1]);  # open the gtf file
while ($l=<F>){
    chomp($l);
    @s=split(/\t/,$l);
    @genestuff=split(/"/,$s[8]);
	#chr1   unknown exon    12613   12721   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS14844";                      
    $length = ceil(abs($s[3]-$s[4]) / 50);
    for ( $i =-50;$i< ($length+50);$i++){
        $place=$i* 50 + 50 *ceil($s[3] / 50) ;
        $myg{ $s[0]."_".$place} = $genestuff[1];
    }
}
close(F);
################################################################ END ANNOTATION

#>13__-_chr2_189862100_189862450_psum_0.153846153846154_n_13_ndiff_9_2663_0.222222222222222_COL3A1_output_Fetal_Stomach_408_GTGAAA_L006_R1.fq
#>4__-_chr5_150269950_65466500_psum_9_nseq_4_numberdiff_4_

while ($l=<FL>){
@s=split(/_/,$l);
$strand=$s[2];
if ($s[6] eq 'psum'){

    $psum=ceil(1000*$s[7])/1000;
    $tot=$s[9];
    $ndiff=$s[11];
$chr=$s[3];
$bin1=$s[4];
$bin2=$s[5];

$usediff=0;
    $id1=$chr."_".$bin1;
    $id2=$chr."_".$bin2;
############# this way we get the two ids...
	$textgene="UNAN";
    if (defined($myg{$id1})){
	$textgene=$myg{$id1};	
    }if (defined($myg{$id2})){
	$textgene=$myg{$id2};	
    }

$seq=<FL>;
chomp($seq);
$seq =~ tr/Z//d;
$seqlength=length($seq);
$gene=$s[14]."-".$psum."-".$tot."-".$ndiff."_seqlength-".$seqlength."-".$textgene;
$n=300-length($seq);
$n1=ceil($n/2);
$n2=300-$n1-length($seq);
$z5= "N"x$n1;
$z3= "N"x$n2;

#chr12|ACTR6:100599522|ACTR6:100598717|rev|+
$low=min($s[4],$s[5])-$boundary; ## TO BE SAFE
$high=max($s[4],$s[5])+$boundary;

$chr=$s[3];
#$hashid=$chr."|".$gene.":".$low."|".$gene.":".$high."|rev|+\n";
$info=$chr."|".$gene.":".$low."|".$gene.":".$high."|rev|".$strand."\n";

    if ($tot>1){
 $myid{$info}="$z5$seq$z3\n";
 if ($strand eq '-'){

     $seq =~ tr/ACGT/TGCA/;
     $seq = reverse $seq;

 $myid{$info}="$z5$seq$z3\n";
}
    }
## assume that the midpoint of the sequence is in the middle of the read and pad with a variety of Ns. could also blat to retreive sequence

}
}
foreach $k (keys %myid){
print ">".$k.$myid{$k};
}
