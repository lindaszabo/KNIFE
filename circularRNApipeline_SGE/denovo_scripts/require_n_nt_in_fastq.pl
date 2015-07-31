### NOT USED ###

##### the step to generate unaligned reads is currently done on SCG3!! but outputs a set of unaligned readscull_unaligned_reads_step1.pl
## file will have R* and be split this way:
# to generate the 3' and 5' reads

###########################################################################
############## examples of genomes
###########################################################################

$nmin = 85; ## number of bases trimmed, 65 for fetal

###########################################################################
############## examples of directories
###########################################################################

$wd="/home/julia/denovo_scripts/testing/Fetal_new_8_21_unaligned/";
$outdir="/home/julia/denovo_scripts/testing/Fetal_new_8_21_unaligned/min_80nt/";

$wd="/home/julia/denovo_scripts/testing/H1EscUnaligned/";

$outdir="/home/julia/denovo_scripts/testing/H1EscUnaligned/min_85nt/";
$wd="/home/julia/alignments/HAP1Alignments/orig/unaligned/";
$outdir="/home/julia/alignments/HAP1Alignments/orig/unaligned/min_80nt/";

#$wd="/home/linda/alignments/OvarianCancer2014_cutAdapt/orig/unaligned/";
#$outdir="/home/linda/alignments/OvarianCancer2014_cutAdapt/orig/unaligned/min_65nt_fastq/";
######################################################################

opendir(my $dh, $wd) || die "can't opendir $wd: $!";
@filesa = grep { /unaligned/ && -f "$wd/$_" } readdir($dh);
@files = grep { //  } @filesa;

print "FILES are @files\n";

closedir $dh;

foreach $f (@files){
$full_path=$wd.$f;
open (F,$full_path);
open (O,">".$outdir.$f);
while ($id=<F>){
$seq=<F>;
$l1=<F>;
$l2=<F>;
if (length($seq)>$nmin){
print O $id;
print O $seq;
print O $l1;
print O $l2;
}
}
}
