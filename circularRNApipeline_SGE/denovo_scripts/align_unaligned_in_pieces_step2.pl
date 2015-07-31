# script takes the unaligned reads and aligns each piece one-by-one

$bowtie = "bowtie ";

$name=$ARGV[0];
$input_file=$ARGV[1];
$reference=$ARGV[2];
$ntrim=$ARGV[3];
$tempOutDir=$ARGV[4];

## FOR EACH INPUT_FILE
$threeprimefile=$tempOutDir."3prime_".$ntrim."_".$name.".out";
$fiveprimefile=$tempOutDir."5prime_".$ntrim."_".$name.".out";

$np="2";

print "OPENING\n input file is $input_file\n";

if (! (-e $threeprimefile)){

# -m 1 means only report unique alignments
# -v 2 means only allow 2 mismatches
# -p is num processors
print  $bowtie ."--chunkmbs 256 -v 2 --un ".$tempOutDir."unaligned_3prime_".$name." --trim5 ".$ntrim." -p $np -m 1 ".$reference." ".$input_file. " > ".$threeprimefile."\n";
print  $bowtie ."--chunkmbs 256 -v 2 --un ".$tempOutDir."unaligned_3prime_".$name." --trim3 ".$ntrim."  -p $np -m 1 ".$reference." ".$input_file. " > ".$fiveprimefile."\n";
system( $bowtie ."--chunkmbs 256 -v 2 --un ".$tempOutDir."unaligned_3prime_".$name." --trim5 ".$ntrim." -p ".$np." -m 1 ".$reference." ".$input_file. " > ".$threeprimefile);
system( $bowtie ."--chunkmbs 256 -v 2 --un ".$tempOutDir."unaligned_5prime_".$name." --trim3 ".$ntrim."  -p ".$np." -m 1 ".$reference." ".$input_file. " > ".$fiveprimefile);
}



