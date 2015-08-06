
use POSIX;

open (F,$ARGV[0]);
open (F2,"");

## IF F2 is specified as a superset of F, it will print out elements of F2 not in F 
## otherwise, prints out elements in F after taking maximum counts over all bins
##  format below is p-num1-num2-seqlength-(number of non-Ns) and then bins 
## NOTE-- for decoy processing, should be rounded up/down to be conservative on the boundary
##### note, the hash lids assumes that there are 2 files input, and the second file contains linear + circular, the first file contains circular only.

$flabel="";#linear+circular";
$f2label="linear";

## these are lower bounds for the distance between bins, otherwise, indels will be called for linears
$flower=0;
$f2lower=200;

#>chr1|-0-3-4_seqlength-102-IRF2BP2:234741900|-0-3-4_seqlength-102-IRF2BP2:234742050|rev|+NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTCTTGCAACTCTTTCAAAATAAAGGACAACAGCAACAACAAACAGAAACAAAACCAATCCCACAATCTCCAAGTTCACCTGGACTGTAACTTCTCTTGCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

    while ($l=<F>){

   if (substr($l,0,1) eq '>'){
       @s=split(/-/,$l);
       $seq=<F>;
       ###NEW
       @s=split(/[_,:,|,-]/,$l);
       $prob=$s[2];
       $numr=$s[3];
       $binid=$s[0].$s[8].$s[16];

	   $length=abs($s[8]-$s[16]);

# if you wanted no Ns  $seq =~ tr/N//d;     
       chomp($l);
       $l=$l."\t$flabel\n";
       if ( ($prob<5)){
       $any{$binid}=1;# hashes on 'any' junction 
if ($length>$flower ){

if ( (!(defined($nums{$binid})))){
      $nums{$binid}=$numr;
       $ids{$binid}=$l.$seq;
}       
if ($numr>$nums{$binid}){
     $nums{$binid}=$numr;
       $ids{$binid}=$l.$seq;
   }
}
   }
}
 }
## NOW OPEN F2 and if any isn't defined, proceed to linears
##########################################################

    while ($l=<F2>){ 
  if (substr($l,0,1) eq '>'){
      @s=split(/-/,$l);
       $seq=<F2>;
       ###NEW
       @s=split(/[_,:,|,-]/,$l);
       $prob=$s[2];
       $numr=$s[3];
       $binid=$s[0].$s[8].$s[16];
	   $length=abs($s[8]-$s[16]);

# if you wanted no Ns:  $seq =~ tr/N//d;     
       chomp($l);
       $l=$l."\t$f2label\n";
       if ((!(defined($any{$binid}))) && ($length>$f2lower ) && ($prob<5)){
       
if ( (!(defined($lnums{$binid})))){
      $lnums{$binid}=$numr;
       $lids{$binid}=$l.$seq;
}       
if ($numr>$lnums{$binid}){
     $lnums{$binid}=$numr;
       $lids{$binid}=$l.$seq;
   }
}
   }

}
  

####################### endf2
##########################################################3

  foreach $k (keys %ids){
       if (($nums{$k}>1 &&(!(substr($ids{$k},0,5) eq '>chrM')))){
	   print "$ids{$k}";
   }
}
################# linears
  foreach $k (keys %lids){
       if (($lnums{$k}>1 &&(!(substr($lids{$k},0,5) eq '>chrM')))){
   
	print "$lids{$k}";
   }
}
