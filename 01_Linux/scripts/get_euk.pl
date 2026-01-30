#!/usr/bin/perl
use 5.010;

# Save the MAG name
$input = @ARGV[0];

# Get the input file
$busco_output = '../busco/' . $input . '/short_summary.specific.eukaryota_odb10.' . $input . '.txt';

# Check MAG quality and copy the good ones
open($in, '<', $busco_output);
while(<$in>) {
    if(/single-copy/) {
        @line = split("\t", $_);
        if(@line[1] >= 255 * 0.7) {
            $completeness = 1;
        }
        else {
            $completeness = 0;
        }
    }
    elsif(/duplicated/) {
        @line = split("\t", $_);
        if(@line[1] <= 255 * 0.1) {
            $contamination = 1;
        }
        else {
            $contamination = 0;
        }
    }
}
close($in);
if($completeness == 1 && $contamination == 1) {
    system("cp $input ../final/");
}
