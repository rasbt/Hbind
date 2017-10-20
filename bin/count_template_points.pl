#!perl
#
# Time-stamp: <Paul Sanschagrin -- Wed Jul 18 10:41:11 EDT 2001>
#
# count_template_points.pl: count the number of each type of point and the
#   number of key/nonkey points in a hbind template
# Usage: count_template_points.pl [-s] <template file>

use Getopt::Std;
use File::Basename;

# Check usage
sub usage {
    $this = basename ($0);
    print "Usage: $this [-s] <template file>\n".
        "Options:  -s:  simplified output (A/D/N/H/Total)\n";
    exit 1;
}
if (getopts('s') == false) {
  usage();
}
if ($opt_s) {
  $simpleflag = 1;
}
usage() unless @ARGV == 1;

# set some defaults
%typenames = (A => "Acceptor",
	      D => "Donor",
	      N => "Doneptor",
	      H => "Hydrophobic");
@typelist = (A,D,N,H);
%typecount = (A => 0,
	      D => 0,
	      N => 0,
	      H => 0);

# Parse CL args and open file
$fn = $ARGV[0];
open FILE, $fn or die "Unable to open file: $fn\n";

# Process the file
while (<FILE>) {
  next if ($_ =~ /^\#/); # Skip comments
  $type = substr($_,0,1);
  $typecount{$type}++;
  if (substr($_,1,1) eq "*") {
    $key{$type}++;
    $key++;
  }
  $num++;
}

# Output
if (!$opt_s) {
  print "# Type        K    NK   Total\n";
  foreach $type (@typelist) {
    printf ("%-12s  %3d  %3d  %3d\n",$typenames{$type},$key{$type},
	    $typecount{$type}-$key{$type},$typecount{$type});
  }
  printf ("%-12s  %3d  %3d  %3d\n","Total",$key,$num-$key,$num);
} else {
   foreach $type (@typelist) {
     print "$typecount{$type}/";
   }
   print "$num (A/D/N/H/Total)\n";
}

