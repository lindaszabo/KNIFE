
$orig=$ARGV[0];
$inputname=$ARGV[1];
$output =">". $ARGV[2];

open (I,$orig);
open(F,$inputname);
#
$i=0;

while ($l=<F>){

chomp($l);
@s=split(/_/,$l);
## 
$bin1=$s[1];
$bin2=$s[2];
$chr=$s[3];
$strand=$s[4];

$psum=$s[6];
$n=$s[8];
$ndiff=$s[10];

if ( ($ndiff>2 )){# && ($bin1<$bin2)) {
    $l=<F>;
$l=<F>;

@mids=split(/\t/,$l);
foreach $id (@mids){
$retreive{$id}=$n."_".$sep1."_".$strand."_".$chr."_".$bin1."_".$bin2."_psum_".$psum."_n_".$n."_ndiff_".$ndiff."_";

}

$i=$i+1;
}
}

################ retreive ids
print "STARTING TO RETRIEVE\n";
print "STARTING TO RETRIEVE\n";
print "STARTING TO RETRIEVE\n";

while ($l=<I>){

chomp($l);
$id=substr($l,1,length($l));

if (  defined($retreive{$id})){

print  $l."\n";
$l=<I>;
if (!(defined($output{$retreive{$id}} ))){
$output{$retreive{$id}}= $l;
}
if (defined($output{$retreive{$id}} )){
$output{$retreive{$id}}= $output{$retreive{$id}}."\n".$l;
}
$l=<I>;
print O $l;
$l=<I>;
print O $l;
}
}

$i=0;
open (O, $output);


foreach $k (keys %output){ 
print "$k is key\n";
@s=split(/\n/,$output{$k}); 
print O ">".$k.$i."\n$s[0]\n"; 
$i=$i+1;
$i=$i+1;

foreach $ss (@s){
$i=$i+1;
}
}
