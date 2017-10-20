#!perl
#
# Time-stamp: <Paul Sanschagrin -- Tue Jul 17 15:16:14 EDT 2001>
#
# unmark_key_template_points.pl: removes the key point marking
#   from all the points in a template file
# Usage: unmark_key_template_points.pl <template file>
#

# check usage
die "Usage: $0 <template file>\n" unless @ARGV == 1;

# parse CL args and open file
$fn = $ARGV[0];
open FILE, $fn or die "Unable to open file: $fn\n";

# process the file
while (<FILE>) {
  substr($_,1,1) = " ";
  print $_;
}

