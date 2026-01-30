#!usr/bin/perl
use 5.010;

# Sort files and get coverage information
$sample_names = 'mealworm_sample_names.txt';
#open($in, '<', $sample_names);
#while(<$in>) {
#  chomp;
#  system("samtools sort -@ 28 -m 4G $_.bam > $_.sorted.bam");
#  system("samtools coverage $_.sorted.bam > $_.coverage.txt");
#  system("sort $_.coverage.txt > $_.coverage.sorted.txt");
#}
#close($in);

# Combine coverage files
$sample_names = 'mealworm_sample_names.txt';
open($in, '<', $sample_names);
while(<$in>) {
  chomp;
  if($. == 1) {
    system("cut -f1,3,4 $_.coverage.sorted.txt | sed 's/numreads/$_/' > mealworm_coverage.txt");
  }
  else{
    system("cut -f4 $_.coverage.sorted.txt | sed 's/numreads/$_/' > temp.txt");
    system("paste mealworm_coverage.txt temp.txt > temp2.txt");
    system("mv temp2.txt mealworm_coverage.txt");
    system("rm temp.txt");
  }
}
close($in);

# Get toal read counts
$sample_names = 'mealworm_sample_names.txt';
open($in, '<', $sample_names);
while(<$in>) {
  chomp;
  $file_name = $_ . '.txt';
  open($in2, '<', $file_name);
  while(<$in2>) {
    if($. == 1) {
      @line = split(' ', $_);
      @line[0] = @line[0] * 2;
      push(@read_count, @line[0]);
    }
  }
  close($in2);
}
close($in);

# Get MAG read depth and convert to reads per million
$MAG_names = 'mealworm_MAG_names.txt';
$output = 'mealworm_MAG_abundances_reads.txt';
system("tail -1 mealworm_coverage.txt | sed 's/\#rname/MAG_name/' | sed 's/\tendpos//' | sed 's/Mealworm_MAGs_combined_//g' > mealworm_MAG_abundances_reads.txt");
open($in, '<', $MAG_names);
while(<$in>) {
  chomp;
  $MAG = $_;
  @total_coverage = (0) x 22;
  @total_length = (0) x 22;
  open($in2, '<', 'mealworm_coverage.txt');
  while(<$in2>) {
    if(/$MAG/) {
      @line = split("\t", $_);
      for($m = 1; $m <= 22; $m++) {
        @total_coverage[$m-1] = (@line[$m+1]) + @total_coverage[$m-1];
        @total_length[$m-1] = @total_length[$m-1] + @line[1];
      }
    }
  }
  close($in2);
  @averages = (0) x 22;
  open($out, '>>', $output);
  print $out ($MAG);
  for($m = 1; $m <= scalar(@read_count); $m++) {
    @averages[$m-1] = 1000000 * @total_coverage[$m-1] / @read_count[$m-1];
    print $out ("\t@averages[$m-1]");
  }
  print $out ("\n");
  close($out);
}
