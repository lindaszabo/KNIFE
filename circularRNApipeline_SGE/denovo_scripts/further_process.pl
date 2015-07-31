### no longer used ###

#$screen = $ARGV[0]; ## this is the input prefix

use POSIX;
#NC_000002.12   Gnomon  exon    40425999        40426149        .       -       .       ID=id115316;Parent=rna9868;Dbxref=GeneID:6546,Genbank:XM_006712086.1,HGNC:11068,HPRD:01659,             
#MIM:182305;gbkey=mRNA;gene=SLC8A1;product=solute carrier family 8 %28sodium%2Fcalcium exchanger%29%2C member 1%2C transcript variant X10;transcript_id=XM_006712086.1                          
#-bash-4.1$ grep SLC8A1  Linda/GRCh38/ref_GRCh38_top_level.gff3|grep exon                                                                                                                    $binsize=50;   

open (F,$ARGV[1]);  # open the gtf file
while ($l=<F>){
#    for ($h=0;$h<10000;$h++){
#	$l=<F>;
    chomp($l);
    @s=split(/\t/,$l);

    @genestuff=split(/"/,$s[8]);
#chr1   unknown exon    12613   12721   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS14844";                                           
    $length = ceil(abs($s[3]-$s[4]) / 50);
    for ( $i =-50;$i< ($length+50);$i++){
        $place=$i* 50 + 50 *ceil($s[3] / 50) ;
        $myg{ $s[0]."_".$place} = $genestuff[1];
#print "putting in". $s[0]."_".$place."\n";
    }
#print "DONE w gene\n";
}

close(F);
$usediff=0;
$wd1="/home/julia/denovo_scripts/testing/";

opendir(my $dh, $wd1) || die "can't opendir $wd: $!";
@filesa = grep { /output_Fetal/ && -f "$wd1/$_" } readdir($dh);
@files = grep { /R1/ } @filesa;
#@files = grep { /output_Fetal/ && -f "$wd1/$_" } readdir($dh);
#@files = grep { /output_HAP1/ && -f "$wd1/$_" } readdir($dh);
#@files = grep { /output_subset_ENC/ && -f "$wd1/$_" } readdir($dh);
#@files = grep { /output_lane/ && -f "$wd1/$_" } readdir($dh);
#@files = grep { /output_140701_HAVERS_0463/ && -f "$wd1/$_" } readdir($dh);
closedir $dh;

foreach $f (@files){

open (INF,$f);
while ($l=<INF>){
    chomp($l);
    @s=split(/_/,$l);
    $seq=<INF>;
#>4__-_chr15_40734300_40734350_psum_0.5_n_4_ndiff_4_0 
    $psum=$s[7]; 
    $ndiff=$s[11]; 
    $chr=$s[3];
    $p01=$s[4];
    $p02=$s[5];
#print "positions are $p01 and $p02 and line $l and f $f\n";

    $p1= 50 * ceil ( $p01/50);
    $p2= 50 * ceil ( $p02/50);
    $chr=$s[3];
    $id1=$chr."_".$p1;
    $id2=$chr."_".$p2;
############# this way we get the two ids...
    $id=$id2;
#print "$id1 and $id2 are the ids\n";

    if (($ndiff> 2)){# was an old filter: && (abs($p01-$p02)>0)){
    if (defined($myg{$id1})){
	$id=$id1;
    }
    $n0=substr($s[0],1,length($s[0]));
    $correctedp=$psum* ($n0/$ndiff);
    $next=$l."_".$correctedp."_".$myg{$id1}."_".$f."\n".$seq;

    if (defined($mycount{$myg{$id}})){
    $mycount{$myg{$id}}=	$mycount{$myg{$id}}+1;
    $mylist{$myg{$id}}=	$mylist{$myg{$id}}.$next;    
    }
if (!(defined($mycount{$myg{$id}}))){
    $mycount{$myg{$id}}=	0;
    $mylist{$myg{$id}}=	$next;    
}
}
}
}

foreach $k ( sort {$a<=>$b} keys %mycount){
#    if (!($k eq "")){

    print "------------COUNT is $mycount{$k} and key is $k end\n";
    print $mylist{$k}."\n";

#}
}
