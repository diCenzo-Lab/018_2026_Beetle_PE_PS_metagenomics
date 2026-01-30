#!usr/bin/perl
use File::Find;
use 5.010;
use Cwd;

$pwd = cwd();
$prefix = @ARGV[0];
$parent = @ARGV[1];
$out_dir = @ARGV[2];
$file_count = @ARGV[3];
$count = 0;


find( \&search_all_folder, $parent );

sub search_all_folder {
	chomp $_;
	return if $_ eq '.' or $_ eq '..';
	&read_files ($_) if (-f);
}

sub read_files {
	($filename) = @_;
	open($fh, '<', $filename);
  $file_count++;
	if($file_count < 10) {
		$output = $pwd . '/' . $out_dir . $prefix . '_000' . $file_count . '.fasta';
		say("$filename\t$output");
		$file_count2 = '000' . $file_count;
	}
	elsif($file_count < 100) {
		$output = $pwd . '/' . $out_dir . $prefix . '_00' . $file_count . '.fasta';
		say("$filename\t$output");
		$file_count2 = '00' . $file_count;
	}
	elsif($file_count < 1000) {
		$output = $pwd . '/' . $out_dir . $prefix . '_0' . $file_count . '.fasta';
		say("$filename\t$output");
		$file_count2 = '0' . $file_count;
	}
	else {
		$output = $pwd . '/' . $out_dir . $prefix . '_' . $file_count . '.fasta';
		say("$filename\t$output");
		$file_count2 = $file_count;
	}
  $count = 0;
  open($out, '>', $output);
  while(<$fh>) {
    if(/>/) {
      $count++;
      if($count < 10) {
        $header = '>' . $prefix . '_' . $file_count2 . '_contig_000' . $count;
        say $out ("$header");
      }
      elsif($count < 100) {
        $header = '>' . $prefix . '_' . $file_count2 . '_contig_00' . $count;
        say $out ("$header");
      }
      elsif($count < 1000) {
        $header = '>' . $prefix . '_' . $file_count2 . '_contig_0' . $count;
        say $out ("$header");
      }
      else {
        $header = '>' . $prefix . '_' . $file_count2 . '_contig_' . $count;
        say $out ("$header");
      }
    }
    else {
      print $out ("$_");
    }
  }
	close($fh);
	close($out);
}
