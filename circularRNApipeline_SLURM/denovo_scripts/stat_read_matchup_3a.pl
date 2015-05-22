######### TAKES TWO ARGUMENTS, $fiveprime and $threeprime files;

$fiveprime= $ARGV[0];
$threeprime= $ARGV[1];

$onlycircles= $ARGV[2];
$input = $ARGV[3];
$screen = $ARGV[4];
$tempOutDir = $ARGV[5];

## input is original file 
############# ONLY FOR VIRAL, can just use the 3rd field (enumerated from 0, and 2nd)
print "files are $fiveprime \t $threeprime \t$onlycircles\t$input\n";
 use POSIX;
$C=1;
open (F,$fiveprime);
my %ids;

while ($l=<F>){
chomp($l);
@s= split(/\t/,$l);
@mms=split(/,/,$s[7]);
$mm=scalar(@mms);
if ($mm<2){
$ids{$s[0]}=$l;

}
}
close(F);
open (F,$threeprime);

while ($l=<F>){

	$addtohash=0;
chomp($l);
@s= split(/\t/,$l);

@mms=split(/,/,$s[7]);
$mm=scalar(@mms);

if (defined($ids{$s[0]}) && ($mm<2)){
	# retreive the other read
	@s1=split(/\t/,$ids{$s[0]});
#### test

	if (($s[2] eq $s1[2]) && (($s[1] eq $s1[1])) && (abs($s[3]- $s1[3]) > 0 ) ) {
	    ##### 3' end is upstream 
$C=50;

######## print "WAS just for cirlces, now for both linear and circular" 
## NOTE: the equality of s[1] and s1[1] requires same sense of orientation


if ( ($onlycircles==0)|| (($s[1] eq '+')  && ( $s[3]<$s1[3] ) && ($onlycircles==1) ) ){

			### put in a hash
			$addtohash=1;	
			$rounder=$C*ceil($s[3]/$C);	
			$rounder1=$C*ceil($s1[3]/$C);	
		}
#### on the - strand
	
if ( ($onlycircles==0) || (($s[1] eq '-')  && ( $s[3]>$s1[3] ) && ($onlycircles==1 ))){
			$addtohash=1;	#
			$rounder=$C*ceil($s[3]/$C);	
			$rounder1=$C*ceil($s1[3]/$C);	
		}
####### splithashtag is the set of all tags with a particular bound.. 
	if ($addtohash==1){
############s1 is the 5' end, from trimming 3'
		$ids{$s[0]}="RETREIVE";
		$tag=$rounder."\t".$rounder1."\t".$s[2]."\t".$s[1];
		$offsets=$s1[3]."-".$s[3];
		$value=$s[3]-$s1[3];
		$offsetvalue=$s[3]."_".$s1[3];
#### seqvalue stores the difference in read position along with the read positiong
		$seqvalue=">starts".$s[3]."_".$s[1]."_5prime-3prime\n".$s1[4].$s[4];
		if ((defined($splithash{$tag}))){
			$splithash{$tag}= $splithash{$tag}."\t".$value;
			$splitseq{$tag}= $splitseq{$tag}."\n".$seqvalue;
			$retreiveseq{$tag}= $retreiveseq{$tag}."\t".$s[0];
 #print "ADDING and nextu $splitseq{$tag}\n";
			}
		if (!(defined($splithash{$tag}))){
			$splithash{$tag}=$value;
			$splitseq{$tag}= $seqvalue;
			$retreiveseq{$tag}= $s[0];
			}
		$tag= $rounder1."\t".$rounder."\t".$s[2];
####################################### info for retreiving sequence: $tag and
		$id=$s[0];
	    $offset=-$s[3]+$rounder; 
	    $offset1=-$s1[3]+$rounder1; 
		$hdkey{$id}="$id\t$offset\t$offset1\t$tag\t$s[1]";

	}
	}
#### test
}
}

################ INSERTED MODULE TO COMPUTE HAMMING DISTS
open (F,$input); ## NOW hamming has enough information to continue to evaluate joint hamming distance.
open (U, "> ".$tempOutDir."hamming_".$screen);
while ($l=<F>){

chomp($l);
$id=substr($l,1,length($l));

if (  defined($hdkey{$id})){
print O $l."\n";
$l=<F>;
print O $l;
@s=split(/\t/,$hdkey{$id}); ### this gives information on offsets
chomp($l);
$rl=length($l);
$ls=scalar(@s)-1;
$strand = $s[$ls];
if ($strand eq '-'){
$fiveprimepad= 1+$s[1] + (200-$rl);
$threeprimepad=100 - ($s[1]);
}
if ($strand eq '+'){
    $fiveprimepad= 100-$s[1];
    $threeprimepad=1 +$s[1]+(200-$rl) ;
}

$n5="Z"x$fiveprimepad;
$n3="Z"x$threeprimepad;
$new=$n5.$l.$n3."\n";
$lnew=length($new);

$a=substr($new,0,50);
$a1=substr($new,50,50);
$a2=substr($new,100,50);
$a3=substr($new,150,50);

print U $hdkey{$id}."\t".$rl."\t".$new;

## this is the read -- print it out

$l=<F>;
print O $l;
$l=<F>;
print O $l;
#}
}
}
##################################################


open (O,">".$tempOutDir."teststatO_".$screen);
foreach $k (sort keys %splitseq){
@s=split(/\n/, $splitseq{$k});

print O "INFO\t$k\t$splithash{$k}\n$retreiveseq{$k}\n$splitseq{$k}\n";

}

open (F,$input);

open (O, ">".$tempOutDir."retreived_all".$screen);

while ($l=<F>){
print $l;
chomp($l);
$id=substr($l,1,length($l));

if (  defined($ids{$id})){
print O $l."\n";
$l=<F>;
print O $l;

## this is the read -- print it out


$l=<F>;
print O $l;
$l=<F>;
print O $l;
}
}
