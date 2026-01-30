#!usr/bin/perl
use 5.010;

# Sort files and get coverage information
$sample_names = 'superworm_sample_names.txt';
open($in, '<', $sample_names);
while(<$in>) {
  chomp;
  system("samtools sort -@ 28 -m 4G $_.bam > $_.sorted.bam");
  system("samtools coverage $_.sorted.bam > $_.coverage.txt");
  system("sort $_.coverage.txt > $_.coverage.sorted.txt");
}
close($in);

# Combine coverage files
$sample_names = 'superworm_sample_names.txt';
open($in, '<', $sample_names);
while(<$in>) {
  chomp;
  if($. == 1) {
    system("cut -f1,3,4 $_.coverage.sorted.txt | sed 's/numreads/$_/' > superworm_coverage.txt");
  }
  else{
    system("cut -f4 $_.coverage.sorted.txt | sed 's/numreads/$_/' > temp.txt");
    system("paste superworm_coverage.txt temp.txt > temp2.txt");
    system("mv temp2.txt superworm_coverage.txt");
    system("rm temp.txt");
  }
}
close($in);

# Get toal read counts
$sample_names = 'superworm_sample_names.txt';
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
$MAG_names = 'superworm_MAG_names.txt';
$output = 'superworm_MAG_abundances_depth.txt';
system("head -1 superworm_coverage.txt | sed 's/\#rname/MAG_name/' | sed 's/\tendpos//' | sed 's/Superworm_MAGs_combined_//g' > superworm_MAG_abundances_depth.txt");
open($in, '<', $MAG_names);
while(<$in>) {
  chomp;
  $MAG = $_;
  @total_reads = (0) x 24;
  @total_length = (0) x 24;
  open($in2, '<', 'superworm_coverage.txt');
  while(<$in2>) {
    if(/$MAG/) {
      @line = split("\t", $_);
      for($m = 1; $m <= 24; $m++) {
        @total_reads[$m-1] = (@line[$m+1]) + @total_reads[$m-1];
        @total_length[$m-1] = @total_length[$m-1] + @line[1];
      }
    }
  }
  close($in2);
  @averages = (0) x 24;
  open($out, '>>', $output);
  print $out ($MAG);
  for($m = 1; $m <= scalar(@read_count); $m++) {
    $depth = @total_reads[$m-1] * 150 / @total_length[$m-1];
    @averages[$m-1] = 1000000 * $depth / @read_count[$m-1];
    print $out ("\t@averages[$m-1]");
  }
  print $out ("\n");
  close($out);
}
