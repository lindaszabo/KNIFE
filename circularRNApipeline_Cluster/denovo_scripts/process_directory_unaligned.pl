### NOTE: DO NOT RUN circonly=1 and circonly=0 for the same dataset at the same time. This will cause overwriting of intermediate files ####

$alignDir = $ARGV[0];
$screen = $ARGV[1];  # screen is a variable used for further processing
$mode = $ARGV[2];  # any string containing mouse, fly, rat, pombe, or crypto will change from default of human genome
$ntrim = $ARGV[3]; # number of nt to trim from each end (2N/3 is good rule of thumb)
$onlycircles = $ARGV[4];  # 0 if want linear and circular, 1 if want only circular in denovo index

# make output directory for this dataset
$tempOutDir=$alignDir."/".$screen."/logs/denovo_script_out/";
$tempResultDir=$alignDir."/".$screen."/logs/denovo_index/";

###########################################################################
############## select reference based on mode
###########################################################################

if (index($mode, "grch38") != -1){
    $reference="index/grch38_genome";
    $gtf="grch38_genes.gtf";
} elsif (index($mode, "mouse") != -1){
    $reference="index/mm10_genome";
    $gtf="mm10_genes.gtf";
} elsif (index($mode, "rat") != -1){
    $reference="index/rn5_genome";
    $gtf="rn5_genes.gtf";
} elsif (index($mode, "fly") != -1){
    $reference="index/dm3_genome";
    $gtf="dm3_genes.gtf";
} elsif (index($mode, "pombe") != -1){
    $reference="index/dm3_genome";
    $gtf="ASM294v2_23_genes.gtf";
} elsif (index($mode, "crypto") != -1){
    $reference="index/cryptococcus_neoformans_grubii_h99_genome";
    $gtf="cryptococcus_neoformans_grubii_h99_genes.gtf";
} elsif (index($mode, "cerevisiae") != -1){
    $reference="index/Scer_genome";
    $gtf="Scer_genes.gtf";
} elsif (index($mode, "mikatae") != -1){
    $reference="index/Smik_genome";
} elsif (index($mode, "bayanus") != -1){
    $reference="index/Sbay_genome";
} elsif (index($mode, "HSV") != -1){
    $reference="index/KOS_genome";
    $gtf="KOS_genes.gtf";
} elsif (index($mode, "capsas") != -1){
    $reference="capsaspora_owczarzaki_atcc_30864_2_genome";
    $gtf="capsaspora_owczarzaki_atcc_30864_2_genes.gtf"; 
} elsif (index($mode, "rosetta") != -1){
    $reference="index/salpingoeca_rosetta_1_genome";
    $gtf="salpingoeca_rosetta_1_genes.gtf";   
} else {
    $reference="index/hg19_genome";
    $gtf="hg19_genes.gtf";
}

###########################################################################
############## requires wd where there are files unaligned_NAME which is propogated to output
###########################################################################

$wd=$alignDir."/".$screen."/orig/unaligned/forDenovoIndex/";

# DO YOU WANT TO USE DIFF CHROM, currently not supported
$diffchrom="";
#$diffchrom="_diffchrom";

opendir(my $dh, $wd) || die "can't opendir $wd: $!";
@filesa = grep { /unaligned/ && -f "$wd/$_" } readdir($dh);
@files = grep { /_R?1/  } @filesa;
print "FILES are @files\n";
closedir $dh;

foreach $f (@files){
$full_path=$wd.$f;
($start,$name)=split(/unaligned_/,$f);

print "start is $start and name is $name";
print "perl align_unaligned_in_pieces_step2.pl $name $full_path $reference $ntrim $tempOutDir";
system ("perl align_unaligned_in_pieces_step2.pl ".$name . " ". $full_path. " ". $reference ." ".$ntrim." ".$tempOutDir);

$threeprimefile=$tempOutDir."3prime_".$ntrim."_".$name.".out";
$fiveprimefile=$tempOutDir."5prime_".$ntrim."_".$name.".out";

print "sample is $sample and COMPLETED ALIGNMENTS files are $threeprimefile and $fiveprimefile and $path$input_file\n";
print "threeprimeinputfile is $threeprimefile\n";
print  " perl stat_read_matchup".$diffchrom."_3a.pl ".$fiveprimefile." ".$threeprimefile." ".$onlycircles." ".$full_path." ".$screen." ".$tempOutDir." > ".$tempOutDir."test".$screen."\n";
system ( " perl  stat_read_matchup".$diffchrom."_3a.pl ".$fiveprimefile." ".$threeprimefile." ".$onlycircles." ".$full_path." ".$screen." ".$tempOutDir." > ".$tempOutDir."test".$screen);

print  "perl test_stats_step2_3c.pl ".$tempOutDir." ".$screen." > ".$tempOutDir."testnew".$screen."\n";
system ( "perl test_stats_step2_3c.pl ".$tempOutDir." ".$screen." > ".$tempOutDir."testnew".$screen);

system ("mv ".$tempOutDir."testnew".$screen." ".$tempOutDir."testnew".$name);

print "perl retreive_represent".$diff.".pl ".$full_path." ".$tempOutDir."testnew".$name." ".$tempOutDir."output_".$name."\n";
system ("perl retreive_represent".$diff.".pl ".$full_path." ".$tempOutDir."testnew".$name." ".$tempOutDir."output_".$name);

}

## rest of commands are run on entire directory:
## renames files properly
system ("mv ".$tempOutDir."lindaOUTPUTEST".$screen." ".$tempOutDir."fasta_for_".$screen);
## generates representatives for each bin and probability, etc.
system ("perl process_representatives_for_lindas_pipeline.pl ".$tempOutDir."fasta_for_".$screen." ".$gtf." > ".$tempOutDir."denovo_".$screen);
## consolidates all of these genes and values to be parsimonious
system ("perl process_max_denovo.pl ".$tempOutDir."denovo_".$screen." > ".$tempResultDir."denovo_".$screen."_onlycircles".$onlycircles.".fa");

print "chdir(".$tempResultDir.") or die\n";
chdir($tempResultDir) or die "$!";
system("bowtie2-build denovo_".$screen."_onlycircles".$onlycircles.".fa denovo_".$screen."_".$onlycircles);

system ("rm ".$tempOutDir."hamming_".$screen);
