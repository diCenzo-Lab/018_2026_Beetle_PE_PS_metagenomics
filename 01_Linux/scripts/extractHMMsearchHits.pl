#!usr/bin/perl
use 5.010;

# variable for the protein fasta file
$proteinSequences = @ARGV[0];
$hmmsearch_hits = @ARGV[1];
$count = 0;

# Make an array with the names of the proteins to extract
$test = -2;
open($in, '<', $hmmsearch_hits);
while(<$in>) {
	chomp;
	$_ =~ s/\ \ /\ /g;
	$_ =~ s/\ \ /\ /g;
	$_ =~ s/\ \ /\ /g;
	$_ =~ s/\ \ /\ /g;
	$_ =~ s/\ \ /\ /g;
	$_ =~ s/\ /,/g;
	$_ =~ s/\[//g;
	$_ =~ s/\]//g;
	@inputLine = split(',',$_);
	$test++;
	if($test > 0) {
		push(@hits,@inputLine[9]);
	}
}
close($in);

@hitsSorted = sort(@hits);

# Extract the protein hits
for($n = 1; $n < 1000000; $n++) {
	open($aa,'<',$proteinSequences);
	while(<$aa>) {
		if($count < $test) {
			@line = split("\t",$_);
			$i = @hitsSorted[$count];
			if(/$i/) {
				print("@line[0]\n@line[1]");
				$count++;
			}
		}
	}
	if($count == $test) {
		$n = 1000000;
	}
}


